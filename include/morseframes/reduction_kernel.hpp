#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <future>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "morseframes/complex_view.hpp"
#include "morseframes/task_executor.hpp"

namespace morseframes {

struct ReductionKernelEvent {
  SimplexId sigma = kInvalidSimplex;
  SimplexId tau = kInvalidSimplex;

  bool is_perforation() const { return tau == kInvalidSimplex; }
};

static_assert(sizeof(ReductionKernelEvent) == 2 * sizeof(SimplexId));
static_assert(std::is_trivially_destructible_v<ReductionKernelEvent>);

enum class ReductionKernelExecutionPolicy {
  Sequential,
  Parallel,
};

struct ReductionKernelExecutionOptions {
  ReductionKernelExecutionPolicy policy =
      ReductionKernelExecutionPolicy::Sequential;
  std::size_t max_workers = 0;
  bool collect_metrics = false;
};

struct ReductionKernelMetrics {
  std::uint64_t facet_nanoseconds = 0;
  std::uint64_t essential_nanoseconds = 0;
  std::uint64_t core_nanoseconds = 0;
  std::uint64_t local_reduction_nanoseconds = 0;
  std::uint64_t aggregation_nanoseconds = 0;
  std::uint64_t merge_nanoseconds = 0;
  std::uint64_t closure_nanoseconds = 0;
  // Initial includes packed preparation. The remaining timers partition
  // sparse preparation work, with dispatch/clock overhead left in closure.
  std::uint64_t closure_initial_nanoseconds = 0;
  std::uint64_t closure_packed_nanoseconds = 0;
  std::uint64_t closure_traversal_nanoseconds = 0;
  std::uint64_t closure_sort_nanoseconds = 0;
  std::uint64_t closure_materialize_nanoseconds = 0;
  std::size_t closure_sparse_cells = 0;
  std::size_t closure_sparse_entries = 0;
  std::size_t closure_boundary_visits = 0;
  std::size_t closure_duplicate_faces = 0;
  std::size_t closure_index_growths = 0;
  std::size_t closure_entry_growths = 0;
  // Per-level elapsed time, including facet dispatch/wait; do not add this to
  // cumulative core/local times, which are nested inside facet execution.
  std::uint64_t facet_execution_nanoseconds = 0;
  std::size_t levels = 0;
  std::size_t kernel_rounds = 0;
  std::size_t facet_kernels = 0;
  std::size_t reductions = 0;
  std::size_t perforations = 0;
  // One batch per parallel facet round, with at most executor_workers tasks.
  std::size_t parallel_batches = 0;
  std::size_t facet_parallel_tasks = 0;
  std::size_t max_parallel_facets = 0;
  std::size_t parallel_level_batches = 0;
  std::size_t max_parallel_levels = 0;
  std::size_t executor_workers = 1;
  std::size_t facet_discovery_parallel_tasks = 0;
  std::size_t essential_parallel_tasks = 0;
  std::size_t aggregation_rounds = 0;
  std::size_t aggregation_parallel_tasks = 0;
  std::size_t facet_discovery_coboundary_visits = 0;
  std::size_t facet_discovery_mask_tests = 0;
  std::size_t incidence_cell_visits = 0;
  std::size_t facet_cell_visits = 0;
  std::size_t local_candidate_visits = 0;
  std::size_t local_coboundary_visits = 0;
  std::size_t local_coboundary_mask_tests = 0;
  std::size_t local_membership_tests = 0;
  std::size_t local_membership_comparisons = 0;
  std::size_t local_large_membership_tests = 0;
  std::size_t local_large_membership_comparisons = 0;
  std::size_t local_sparse_scan_passes = 0;
  std::size_t local_sparse_candidate_visits = 0;
  std::size_t local_removed_candidate_visits = 0;
  std::size_t local_protected_candidate_visits = 0;
  std::size_t inline_cell_overflows = 0;
  std::size_t inline_event_overflows = 0;
};

struct ReductionKernelLevelResult {
  std::vector<ReductionKernelEvent> events;
  ReductionKernelMetrics metrics;
};

// Mutable scratch space for Algorithm 1 on immutable ComplexView data. Facet
// workers only read the round snapshot and return isolated results; active-set
// mutation is confined to the deterministic coordinator merge.
template <class ComplexView>
class ReductionKernelWorkspace {
  static_assert(is_complex_view_v<ComplexView>,
                "ReductionKernelWorkspace requires a Morse complex-view type.");

 private:
  template <typename View, typename = void>
  struct HasSameLevelClosureCache : std::false_type {};

  template <typename View>
  struct HasSameLevelClosureCache<
      View,
      std::void_t<
          decltype(std::declval<const View&>()
                       .has_same_level_closure_cache()),
          decltype(std::declval<const View&>()
                       .same_level_closure_entries()),
          decltype(std::declval<const View&>()
                       .same_level_closure_ranges()),
          decltype(std::declval<const View&>()
                       .same_level_coboundary_entries()),
          decltype(std::declval<const View&>()
                       .same_level_coboundary_ranges())>> : std::true_type {};

  using Clock = std::chrono::steady_clock;
  // A tetrahedron has 15 nonempty faces and admits at most seven local pairs.
  static constexpr std::size_t kInlineCellCapacity = 16;
  static constexpr std::size_t kInlineEventCapacity = 8;
  static constexpr std::size_t kPackedClosureBucketCapacity = 128;
  // Facet discovery only tests immediate cofaces. Give each submitted task
  // enough active candidates to amortize dispatch/wait, including shrinking
  // rounds on large plateaus. This does not gate local facet execution.
  static constexpr std::size_t kMinFacetDiscoveryTaskSize = 4096;
  // A cheap work proxy, not a time prediction: sparse closure entries visited,
  // live packed cells, or full-list candidates in the graph fallback. Avoid
  // dispatch/wait for small rounds and scale the task budget as rounds shrink.
  static constexpr std::size_t kMinFacetExecutionTaskWork = 1024;
  using PackedMask =
      std::array<std::uint64_t, (kPackedClosureBucketCapacity + 63) / 64>;

  template <typename T, std::size_t InlineCapacity>
  class InlineVector {
   public:
    void push_back(T value) {
      if (!using_overflow_ && size_ < InlineCapacity) {
        inline_entries_[size_++] = std::move(value);
        return;
      }
      if (!using_overflow_) {
        overflow_entries_.reserve(2 * InlineCapacity);
        for (std::size_t index = 0; index < size_; ++index) {
          overflow_entries_.push_back(std::move(inline_entries_[index]));
        }
        using_overflow_ = true;
      }
      overflow_entries_.push_back(std::move(value));
      size_ = overflow_entries_.size();
    }

    std::size_t size() const { return size_; }
    bool uses_overflow() const { return using_overflow_; }

    const T& operator[](std::size_t index) const {
      return using_overflow_ ? overflow_entries_[index]
                             : inline_entries_[index];
    }

    template <bool CountComparisons = false>
    std::size_t index_of(const T& value, std::size_t* comparisons = nullptr) const {
      for (std::size_t index = 0; index < size_; ++index) {
        if constexpr (CountComparisons) {
          ++*comparisons;
        }
        if ((*this)[index] == value) {
          return index;
        }
      }
      return size_;
    }

    template <typename Destination>
    void append_to(Destination& destination) const {
      destination.reserve(destination.size() + size_);
      for (std::size_t index = 0; index < size_; ++index) {
        destination.push_back((*this)[index]);
      }
    }

   private:
    std::array<T, InlineCapacity> inline_entries_{};
    std::vector<T> overflow_entries_;
    std::size_t size_ = 0;
    bool using_overflow_ = false;
  };

  struct EmptyFacetKernelDiagnostics {};

  struct FacetKernelDiagnostics {
    std::uint64_t core_nanoseconds = 0;
    std::uint64_t local_reduction_nanoseconds = 0;
    std::size_t facet_cell_visits = 0;
    std::size_t local_candidate_visits = 0;
    std::size_t local_coboundary_visits = 0;
    std::size_t local_coboundary_mask_tests = 0;
    std::size_t local_membership_tests = 0;
    std::size_t local_membership_comparisons = 0;
    std::size_t local_large_membership_tests = 0;
    std::size_t local_large_membership_comparisons = 0;
    std::size_t local_sparse_scan_passes = 0;
    std::size_t local_sparse_candidate_visits = 0;
    std::size_t local_removed_candidate_visits = 0;
    std::size_t local_protected_candidate_visits = 0;
    std::size_t inline_cell_overflows = 0;
    std::size_t inline_event_overflows = 0;
  };

  template <bool CollectMetrics>
  struct FacetKernelResult
      : std::conditional_t<CollectMetrics, FacetKernelDiagnostics,
                           EmptyFacetKernelDiagnostics> {
    InlineVector<ReductionKernelEvent, kInlineEventCapacity> events;
  };

  static_assert(sizeof(FacetKernelResult<false>) <
                sizeof(FacetKernelResult<true>));

  struct LevelCells {
    bool enabled = false;
    std::vector<SimplexId> entries;
    std::vector<std::pair<std::size_t, std::size_t>> ranges;
    const std::vector<std::uint64_t>* packed_masks = nullptr;
    const std::vector<std::uint64_t>* packed_cofaces = nullptr;
    const std::vector<SimplexId>* packed_bucket = nullptr;
    std::size_t packed_block_count = 0;
    PackedMask packed_active{};
    PackedMask packed_unique{};
    const std::vector<SimplexId>* cached_entries = nullptr;
    const std::vector<std::pair<std::size_t, std::size_t>>* cached_ranges =
        nullptr;
  };

  struct LevelScratch {
    void prepare(std::size_t bucket_size, bool collect_metrics) {
      facet_flags.resize(bucket_size);
      if (collect_metrics) {
        coboundary_visits.resize(bucket_size);
      } else {
        coboundary_visits.clear();
      }
      facets.clear();
      if (facets.capacity() < bucket_size) {
        facets.reserve(bucket_size);
      }
      active_simplices.clear();
      if (active_simplices.capacity() < bucket_size) {
        active_simplices.reserve(bucket_size);
      }
      facet_results.clear();
      diagnostic_facet_results.clear();
      round_events.clear();
      if (round_events.capacity() < bucket_size) {
        round_events.reserve(bucket_size);
      }
      included.resize(bucket_size);
      closure_masks.clear();
      coface_masks.clear();
      cell_indices.clear();
      if (cell_indices.capacity() < kInlineCellCapacity) {
        cell_indices.reserve(kInlineCellCapacity);
      }
      level_cells.enabled = false;
      level_cells.entries.clear();
      level_cells.ranges.clear();
      level_cells.packed_masks = nullptr;
      level_cells.packed_cofaces = nullptr;
      level_cells.packed_bucket = nullptr;
      level_cells.packed_block_count = 0;
      level_cells.packed_active.fill(0);
      level_cells.packed_unique.fill(0);
      level_cells.cached_entries = nullptr;
      level_cells.cached_ranges = nullptr;
    }

    std::vector<std::uint8_t> facet_flags;
    std::vector<std::size_t> coboundary_visits;
    std::vector<SimplexId> facets;
    std::vector<SimplexId> active_simplices;
    std::vector<FacetKernelResult<false>> facet_results;
    std::vector<FacetKernelResult<true>> diagnostic_facet_results;
    std::vector<ReductionKernelEvent> round_events;
    LevelCells level_cells;
    std::vector<std::uint8_t> included;
    std::vector<std::size_t> cell_indices;
    std::vector<std::uint64_t> closure_masks;
    std::vector<std::uint64_t> coface_masks;
  };

  class FixedEventBuffer {
   public:
    FixedEventBuffer(ReductionKernelEvent* entries, std::size_t capacity)
        : entries_(entries), capacity_(capacity) {}

    void reserve(std::size_t requested) const {
      if (requested > capacity_) {
        throw std::length_error(
            "Reduction-kernel level event capacity exceeded.");
      }
    }

    void push_back(const ReductionKernelEvent& event) {
      reserve(size_ + 1);
      ::new (static_cast<void*>(entries_ + size_)) ReductionKernelEvent(event);
      ++size_;
    }

    std::size_t size() const { return size_; }

   private:
    ReductionKernelEvent* entries_ = nullptr;
    std::size_t capacity_ = 0;
    std::size_t size_ = 0;
  };

  static std::uint64_t elapsed_nanoseconds(Clock::time_point start,
                                           Clock::time_point stop) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start)
            .count());
  }

  static std::size_t trailing_zero_count(std::uint64_t bits) {
#if defined(__clang__) || defined(__GNUC__)
    return static_cast<std::size_t>(
        __builtin_ctzll(static_cast<unsigned long long>(bits)));
#else
    std::size_t count = 0;
    while ((bits & std::uint64_t{1}) == 0) {
      bits >>= 1;
      ++count;
    }
    return count;
#endif
  }

  static std::size_t population_count(std::uint64_t bits) {
#if defined(__clang__) || defined(__GNUC__)
    return static_cast<std::size_t>(
        __builtin_popcountll(static_cast<unsigned long long>(bits)));
#else
    std::size_t count = 0;
    while (bits != 0) {
      bits &= bits - 1;
      ++count;
    }
    return count;
#endif
  }

  template <bool CollectMetrics>
  Clock::time_point profile_start() const {
    if constexpr (CollectMetrics) {
      return Clock::now();
    }
    return Clock::time_point{};
  }

  template <bool CollectMetrics>
  void profile_add(std::uint64_t& destination,
                   Clock::time_point start) const {
    if constexpr (CollectMetrics) {
      destination += elapsed_nanoseconds(start, Clock::now());
    }
  }

 public:
  explicit ReductionKernelWorkspace(
      const ComplexView& complex,
      ReductionKernelExecutionOptions options = {},
      std::shared_ptr<BoundedTaskExecutor> executor = {})
      : complex_(complex),
        options_(options),
        executor_(std::move(executor)),
        active_(complex.size(), 0),
        facet_incidence_(complex.size(), 0),
        round_removed_(complex.size(), 0),
        bucket_index_(complex.size(), std::numeric_limits<std::size_t>::max()) {
    if (options_.policy == ReductionKernelExecutionPolicy::Parallel &&
        executor_ == nullptr) {
      executor_ =
          std::make_shared<BoundedTaskExecutor>(options_.max_workers);
    }
    const std::size_t scratch_count =
        executor_ == nullptr ? 1 : executor_->worker_count();
    level_scratch_.resize(scratch_count);
  }

  const ReductionKernelMetrics& metrics() const { return metrics_; }

  std::vector<ReductionKernelEvent> compute_level(LevelId level) {
    auto result = compute_level_isolated(level);
    accumulate_metrics(metrics_, result.metrics);
    return std::move(result.events);
  }

  // Each level owns disjoint entries of active_ and round_removed_. This
  // isolated form therefore supports Algorithm 2 level tasks without sharing
  // event buffers or counters between workers.
  ReductionKernelLevelResult compute_level_isolated(
      LevelId level, bool allow_intra_level_parallelism = true) {
    LevelScratch scratch;
    ReductionKernelLevelResult result;
    result.metrics = dispatch_level_with_scratch(
        level, allow_intra_level_parallelism, scratch, result.events);
    return result;
  }

  ReductionKernelLevelResult compute_level_isolated_reusing_scratch(
      LevelId level, std::size_t scratch_index,
      bool allow_intra_level_parallelism = true) {
    // A caller may process levels concurrently only when each long-lived task
    // owns a distinct scratch index.
    if (scratch_index >= level_scratch_.size()) {
      throw std::out_of_range(
          "Reduction-kernel scratch index exceeds executor workers.");
    }
    ReductionKernelLevelResult result;
    result.metrics = dispatch_level_with_scratch(
        level, allow_intra_level_parallelism, level_scratch_[scratch_index],
        result.events);
    return result;
  }

  ReductionKernelMetrics compute_level_isolated_into(
      LevelId level, std::size_t scratch_index,
      ReductionKernelEvent* event_storage, std::size_t event_capacity,
      std::size_t& event_count,
      bool allow_intra_level_parallelism = true) {
    // The caller owns this level's disjoint event slice and the scratch index
    // assigned to its long-lived task.
    if (scratch_index >= level_scratch_.size()) {
      throw std::out_of_range(
          "Reduction-kernel scratch index exceeds executor workers.");
    }
    FixedEventBuffer events(event_storage, event_capacity);
    auto metrics = dispatch_level_with_scratch(
        level, allow_intra_level_parallelism, level_scratch_[scratch_index],
        events);
    event_count = events.size();
    return metrics;
  }

  void compute_level_isolated_into_unprofiled(
      LevelId level, std::size_t scratch_index,
      ReductionKernelEvent* event_storage, std::size_t event_capacity,
      std::size_t& event_count,
      bool allow_intra_level_parallelism = true) {
    // This entry point keeps the ordinary construction path free of a
    // returned metrics object and instantiates only the metrics-free kernel.
    if (scratch_index >= level_scratch_.size()) {
      throw std::out_of_range(
          "Reduction-kernel scratch index exceeds executor workers.");
    }
    FixedEventBuffer events(event_storage, event_capacity);
    ReductionKernelMetrics ignored;
    compute_level_isolated_with_scratch<false>(
        level, allow_intra_level_parallelism, level_scratch_[scratch_index],
        events, ignored);
    event_count = events.size();
  }

 private:
  template <typename EventBuffer>
  ReductionKernelMetrics dispatch_level_with_scratch(
      LevelId level, bool allow_intra_level_parallelism,
      LevelScratch& scratch, EventBuffer& events) {
    ReductionKernelMetrics metrics;
    if (options_.collect_metrics) {
      compute_level_isolated_with_scratch<true>(
          level, allow_intra_level_parallelism, scratch, events, metrics);
    } else {
      compute_level_isolated_with_scratch<false>(
          level, allow_intra_level_parallelism, scratch, events, metrics);
    }
    return metrics;
  }

  template <bool CollectMetrics, typename EventBuffer>
  void compute_level_isolated_with_scratch(
      LevelId level, bool allow_intra_level_parallelism,
      LevelScratch& scratch, EventBuffer& events,
      ReductionKernelMetrics& metrics) {
    const auto& bucket = complex_.simplices_of_level(level);
    if constexpr (CollectMetrics) {
      metrics.executor_workers =
          executor_ == nullptr ? 1 : executor_->worker_count();
    }
    events.reserve(bucket.size());
    std::size_t remaining = bucket.size();
    if constexpr (CollectMetrics) {
      ++metrics.levels;
    }

    for (std::size_t index = 0; index < bucket.size(); ++index) {
      const SimplexId simplex = bucket[index];
      active_[simplex] = 1;
      bucket_index_[simplex] = index;
    }
    const auto closure_start = profile_start<CollectMetrics>();
    const bool cache_level_cells =
        std::any_of(bucket.begin(), bucket.end(), [&](SimplexId simplex) {
          return complex_.dimension(simplex) >= 2;
        });
    scratch.prepare(bucket.size(), CollectMetrics);
    scratch.active_simplices.assign(bucket.begin(), bucket.end());
    auto& level_cells = scratch.level_cells;
    build_level_cells<CollectMetrics>(bucket, cache_level_cells, scratch, level_cells, metrics);
    if constexpr (CollectMetrics) {
      const auto elapsed = elapsed_nanoseconds(closure_start, Clock::now());
      metrics.closure_nanoseconds += elapsed;
      metrics.closure_initial_nanoseconds += elapsed;
    }

    while (remaining > 0) {
      bool kernel_round_changed = false;

      do {
        if constexpr (CollectMetrics) {
          ++metrics.kernel_rounds;
        }
        kernel_round_changed = false;
        const auto facet_start = profile_start<CollectMetrics>();
        const auto& facets = active_facets<CollectMetrics>(
            level, bucket, remaining, scratch, metrics,
            allow_intra_level_parallelism);
        profile_add<CollectMetrics>(metrics.facet_nanoseconds, facet_start);
        // Complete all cell writes before incidence and local facet tasks read
        // them. Concurrent levels use disjoint scratch and bucket indices.
        const auto facet_closure_start = profile_start<CollectMetrics>();
        prepare_facet_cells<CollectMetrics>(facets, bucket, scratch, level_cells, metrics);
        profile_add<CollectMetrics>(metrics.closure_nanoseconds,
                                    facet_closure_start);
        if constexpr (CollectMetrics) {
          metrics.facet_kernels += facets.size();
        }
        const auto essential_start = profile_start<CollectMetrics>();
        compute_facet_incidence<CollectMetrics>(
            facets, scratch.active_simplices, level_cells, metrics);
        profile_add<CollectMetrics>(metrics.essential_nanoseconds,
                                    essential_start);
        scratch.round_events.clear();
        auto record_facet_events = [&](const auto& result) {
          for (std::size_t index = 0; index < result.events.size(); ++index) {
            const auto& event = result.events[index];
            if (event.is_perforation()) {
              throw std::logic_error(
                  "A facet kernel returned a non-reduction event.");
            }
            if (round_removed_[event.sigma] || round_removed_[event.tau]) {
              throw std::logic_error(
                  "Facet reduction kernels removed the same simplex twice.");
            }
            round_removed_[event.sigma] = 1;
            round_removed_[event.tau] = 1;
            scratch.round_events.push_back(event);
          }
        };

        const auto execution_start = profile_start<CollectMetrics>();
        const std::size_t facet_tasks = facet_task_count(
            facets, scratch.active_simplices, level_cells,
            allow_intra_level_parallelism);
        if constexpr (!CollectMetrics) {
          if (facet_tasks <= 1) {
            for (SimplexId facet : facets) {
              const auto result = compute_facet_kernel<false>(
                  level, facet, scratch.active_simplices, level_cells);
              record_facet_events(result);
            }
          } else {
            const auto& facet_results = execute_facets<false>(
                level, facets, scratch.active_simplices, level_cells,
                scratch.facet_results, metrics, facet_tasks);
            for (const auto& result : facet_results) {
              record_facet_events(result);
            }
          }
        } else {
          const auto& facet_results = execute_facets<true>(
              level, facets, scratch.active_simplices, level_cells,
              scratch.diagnostic_facet_results, metrics, facet_tasks);
          profile_add<true>(metrics.facet_execution_nanoseconds, execution_start);
          const auto aggregation_start = profile_start<true>();
          if (facet_results.size() > 1) {
            ++metrics.aggregation_rounds;
          }
          for (const auto& facet_result : facet_results) {
            metrics.core_nanoseconds += facet_result.core_nanoseconds;
            metrics.local_reduction_nanoseconds +=
                facet_result.local_reduction_nanoseconds;
            metrics.facet_cell_visits += facet_result.facet_cell_visits;
            metrics.local_candidate_visits +=
                facet_result.local_candidate_visits;
            metrics.local_coboundary_visits +=
                facet_result.local_coboundary_visits;
            metrics.local_coboundary_mask_tests +=
                facet_result.local_coboundary_mask_tests;
            metrics.local_membership_tests +=
                facet_result.local_membership_tests;
            metrics.local_membership_comparisons += facet_result.local_membership_comparisons;
            metrics.local_large_membership_tests += facet_result.local_large_membership_tests;
            metrics.local_large_membership_comparisons += facet_result.local_large_membership_comparisons;
            metrics.local_sparse_scan_passes += facet_result.local_sparse_scan_passes;
            metrics.local_sparse_candidate_visits += facet_result.local_sparse_candidate_visits;
            metrics.local_removed_candidate_visits += facet_result.local_removed_candidate_visits;
            metrics.local_protected_candidate_visits += facet_result.local_protected_candidate_visits;
            metrics.inline_cell_overflows +=
                facet_result.inline_cell_overflows;
            metrics.inline_event_overflows +=
                facet_result.inline_event_overflows;
            record_facet_events(facet_result);
          }
          profile_add<true>(metrics.aggregation_nanoseconds,
                            aggregation_start);
        }

        const auto merge_start = profile_start<CollectMetrics>();
        for (const auto& event : scratch.round_events) {
          if (!active_[event.sigma] || !active_[event.tau]) {
            throw std::logic_error(
                "A reduction-kernel round removed an inactive simplex.");
          }
          active_[event.sigma] = 0;
          active_[event.tau] = 0;
          if (level_cells.packed_masks != nullptr) {
            for (SimplexId simplex : {event.sigma, event.tau}) {
              const std::size_t index = bucket_index_[simplex];
              level_cells.packed_active[index / 64] &=
                  ~(std::uint64_t{1} << (index % 64));
            }
          }
          round_removed_[event.sigma] = 0;
          round_removed_[event.tau] = 0;
          remaining -= 2;
          kernel_round_changed = true;
          events.push_back(event);
          if constexpr (CollectMetrics) {
            ++metrics.reductions;
          }
        }
        profile_add<CollectMetrics>(metrics.merge_nanoseconds, merge_start);
      } while (kernel_round_changed);

      if (remaining == 0) {
        break;
      }

      const auto facet_start = profile_start<CollectMetrics>();
      const auto& facets = active_facets<CollectMetrics>(
          level, bucket, remaining, scratch, metrics,
          allow_intra_level_parallelism);
      profile_add<CollectMetrics>(metrics.facet_nanoseconds, facet_start);
      if (facets.empty()) {
        throw std::logic_error(
            "A nonempty reduction-kernel section has no facet.");
      }
      const SimplexId critical = facets.front();
      events.push_back(ReductionKernelEvent{critical, kInvalidSimplex});
      active_[critical] = 0;
      if (level_cells.packed_masks != nullptr) {
        const std::size_t index = bucket_index_[critical];
        level_cells.packed_active[index / 64] &=
            ~(std::uint64_t{1} << (index % 64));
      }
      --remaining;
      if constexpr (CollectMetrics) {
        ++metrics.perforations;
      }
    }

  }

 public:
  static void accumulate_metrics(ReductionKernelMetrics& destination,
                                 const ReductionKernelMetrics& source) {
    destination.facet_nanoseconds += source.facet_nanoseconds;
    destination.essential_nanoseconds += source.essential_nanoseconds;
    destination.core_nanoseconds += source.core_nanoseconds;
    destination.local_reduction_nanoseconds +=
        source.local_reduction_nanoseconds;
    destination.aggregation_nanoseconds += source.aggregation_nanoseconds;
    destination.merge_nanoseconds += source.merge_nanoseconds;
    destination.closure_nanoseconds += source.closure_nanoseconds;
    destination.closure_initial_nanoseconds += source.closure_initial_nanoseconds;
    destination.closure_packed_nanoseconds += source.closure_packed_nanoseconds;
    destination.closure_traversal_nanoseconds += source.closure_traversal_nanoseconds;
    destination.closure_sort_nanoseconds += source.closure_sort_nanoseconds;
    destination.closure_materialize_nanoseconds += source.closure_materialize_nanoseconds;
    destination.closure_sparse_cells += source.closure_sparse_cells;
    destination.closure_sparse_entries += source.closure_sparse_entries;
    destination.closure_boundary_visits += source.closure_boundary_visits;
    destination.closure_duplicate_faces += source.closure_duplicate_faces;
    destination.closure_index_growths += source.closure_index_growths;
    destination.closure_entry_growths += source.closure_entry_growths;
    destination.facet_execution_nanoseconds += source.facet_execution_nanoseconds;
    destination.levels += source.levels;
    destination.kernel_rounds += source.kernel_rounds;
    destination.facet_kernels += source.facet_kernels;
    destination.reductions += source.reductions;
    destination.perforations += source.perforations;
    destination.parallel_batches += source.parallel_batches;
    destination.facet_parallel_tasks += source.facet_parallel_tasks;
    destination.max_parallel_facets =
        std::max(destination.max_parallel_facets,
                 source.max_parallel_facets);
    destination.parallel_level_batches += source.parallel_level_batches;
    destination.max_parallel_levels =
        std::max(destination.max_parallel_levels,
                 source.max_parallel_levels);
    destination.executor_workers =
        std::max(destination.executor_workers, source.executor_workers);
    destination.facet_discovery_parallel_tasks +=
        source.facet_discovery_parallel_tasks;
    destination.essential_parallel_tasks += source.essential_parallel_tasks;
    destination.aggregation_rounds += source.aggregation_rounds;
    destination.aggregation_parallel_tasks +=
        source.aggregation_parallel_tasks;
    destination.facet_discovery_coboundary_visits +=
        source.facet_discovery_coboundary_visits;
    destination.facet_discovery_mask_tests +=
        source.facet_discovery_mask_tests;
    destination.incidence_cell_visits += source.incidence_cell_visits;
    destination.facet_cell_visits += source.facet_cell_visits;
    destination.local_candidate_visits += source.local_candidate_visits;
    destination.local_coboundary_visits += source.local_coboundary_visits;
    destination.local_coboundary_mask_tests += source.local_coboundary_mask_tests;
    destination.local_membership_tests += source.local_membership_tests;
    destination.local_membership_comparisons += source.local_membership_comparisons;
    destination.local_large_membership_tests += source.local_large_membership_tests;
    destination.local_large_membership_comparisons += source.local_large_membership_comparisons;
    destination.local_sparse_scan_passes += source.local_sparse_scan_passes;
    destination.local_sparse_candidate_visits += source.local_sparse_candidate_visits;
    destination.local_removed_candidate_visits += source.local_removed_candidate_visits;
    destination.local_protected_candidate_visits += source.local_protected_candidate_visits;
    destination.inline_cell_overflows += source.inline_cell_overflows;
    destination.inline_event_overflows += source.inline_event_overflows;
  }

 private:
  bool is_face_of(SimplexId face, SimplexId simplex) const {
    const auto& face_vertices = complex_.vertices(face);
    const auto& simplex_vertices = complex_.vertices(simplex);
    return std::includes(simplex_vertices.begin(), simplex_vertices.end(),
                         face_vertices.begin(), face_vertices.end());
  }

  const std::vector<SimplexId>& cell_entries(
      const LevelCells& cells) const {
    return cells.cached_entries == nullptr ? cells.entries
                                           : *cells.cached_entries;
  }

  std::pair<std::size_t, std::size_t> cell_range(
      const LevelCells& cells, SimplexId simplex) const {
    if (cells.cached_ranges != nullptr) {
      return (*cells.cached_ranges)[simplex];
    }
    return cells.ranges[bucket_index_[simplex]];
  }

  bool use_precomputed_cache() const {
    if constexpr (HasSameLevelClosureCache<ComplexView>::value) {
      const bool parallel_levels =
          options_.policy == ReductionKernelExecutionPolicy::Parallel &&
          executor_ != nullptr && executor_->worker_count() > 1;
      return !parallel_levels && complex_.has_same_level_closure_cache();
    }
    return false;
  }

  template <typename Function>
  void visit_same_level_coboundary(
      SimplexId simplex, LevelId level, bool use_cache,
      Function&& function) const {
    if constexpr (HasSameLevelClosureCache<ComplexView>::value) {
      if (use_cache) {
        const auto& ranges = complex_.same_level_coboundary_ranges();
        const auto& entries = complex_.same_level_coboundary_entries();
        const auto [first, last] = ranges[simplex];
        for (std::size_t index = first; index < last; ++index) {
          if (!function(entries[index])) {
            return;
          }
        }
        return;
      }
    }
    for (SimplexId coface : complex_.coboundary(simplex)) {
      if (complex_.level(coface) == level && !function(coface)) {
        return;
      }
    }
  }

  template <bool CollectMetrics>
  void build_level_cells(
      const std::vector<SimplexId>& bucket, bool enabled,
      LevelScratch& scratch, LevelCells& cells, ReductionKernelMetrics& metrics) const {
    cells.enabled = enabled;
    if (!enabled) {
      return;
    }
    if constexpr (HasSameLevelClosureCache<ComplexView>::value) {
      if (use_precomputed_cache()) {
        cells.cached_entries = &complex_.same_level_closure_entries();
        cells.cached_ranges = &complex_.same_level_closure_ranges();
        return;
      }
    }
    if (bucket.size() <= kPackedClosureBucketCapacity) {
      const auto packed_start = profile_start<CollectMetrics>();
      const std::size_t block_count = (bucket.size() + 63) / 64;
      auto& masks = scratch.closure_masks;
      masks.assign(bucket.size() * block_count, 0);
      auto& cofaces = scratch.coface_masks;
      cofaces.assign(bucket.size() * block_count, 0);
      cells.packed_masks = &masks;
      cells.packed_cofaces = &cofaces;
      cells.packed_bucket = &bucket;
      cells.packed_block_count = block_count;
      for (std::size_t simplex_index = 0; simplex_index < bucket.size();
           ++simplex_index) {
        const SimplexId simplex = bucket[simplex_index];
        auto* simplex_mask = masks.data() + simplex_index * block_count;
        simplex_mask[simplex_index / 64] |=
            std::uint64_t{1} << (simplex_index % 64);
        cells.packed_active[simplex_index / 64] |=
            std::uint64_t{1} << (simplex_index % 64);
        for (SimplexId face : complex_.boundary(simplex)) {
          if (complex_.level(face) != complex_.level(simplex)) {
            continue;
          }
          const std::size_t face_index = bucket_index_[face];
          if (face_index >= simplex_index) {
            throw std::logic_error(
                "Reduction-kernel level bucket is not face-first.");
          }
          // Only immediate cofaces belong in this mask, not all containing
          // simplices. The same boundary visit also builds transitive closure.
          cofaces[face_index * block_count + simplex_index / 64] |=
              std::uint64_t{1} << (simplex_index % 64);
          const auto* face_mask = masks.data() + face_index * block_count;
          for (std::size_t block = 0; block < block_count; ++block) {
            simplex_mask[block] |= face_mask[block];
          }
        }
      }
      profile_add<CollectMetrics>(metrics.closure_packed_nanoseconds, packed_start);
      return;
    }
    if (cells.entries.capacity() < 4 * bucket.size()) {
      cells.entries.reserve(4 * bucket.size());
    }
    // Sparse closures are needed only for simplices exposed as facets. An
    // empty range marks an unprepared cell (every actual cell contains itself).
    cells.ranges.assign(bucket.size(), {0, 0});
    std::fill(scratch.included.begin(), scratch.included.end(), 0);
  }

  template <bool CollectMetrics>
  void prepare_facet_cells(
      const std::vector<SimplexId>& facets,
      const std::vector<SimplexId>& bucket,
      LevelScratch& scratch, LevelCells& cells, ReductionKernelMetrics& metrics) const {
    if (!cells.enabled || cells.packed_masks != nullptr ||
        cells.cached_entries != nullptr) {
      return;
    }
    auto& included = scratch.included;
    auto& cell_indices = scratch.cell_indices;
    for (SimplexId facet : facets) {
      const std::size_t facet_index = bucket_index_[facet];
      if (cells.ranges[facet_index].first != cells.ranges[facet_index].second) {
        continue;
      }
      const std::size_t first = cells.entries.size();
      const auto traversal_start = profile_start<CollectMetrics>();
      if constexpr (CollectMetrics) {
        ++metrics.closure_sparse_cells;
      }
      cell_indices.clear();
      cell_indices.push_back(facet_index);
      included[facet_index] = 1;
      // Visit each same-level face once, including inactive faces so the cell
      // remains the same immutable closure as in the eager implementation.
      for (std::size_t next = 0; next < cell_indices.size(); ++next) {
        const auto simplex_index = cell_indices[next];
        for (SimplexId face : complex_.boundary(bucket[simplex_index])) {
          if constexpr (CollectMetrics) {
            ++metrics.closure_boundary_visits;
          }
          if (complex_.level(face) != complex_.level(facet)) {
            continue;
          }
          const std::size_t face_index = bucket_index_[face];
          if (face_index >= simplex_index) {
            throw std::logic_error(
                "Reduction-kernel level bucket is not face-first.");
          }
          if (!included[face_index]) {
            included[face_index] = 1;
            if constexpr (CollectMetrics) {
              metrics.closure_index_growths += cell_indices.size() == cell_indices.capacity();
            }
            cell_indices.push_back(face_index);
          } else if constexpr (CollectMetrics) {
            ++metrics.closure_duplicate_faces;
          }
        }
      }
      profile_add<CollectMetrics>(metrics.closure_traversal_nanoseconds, traversal_start);
      if constexpr (CollectMetrics) {
        metrics.closure_sparse_entries += cell_indices.size();
      }
      const auto sort_start = profile_start<CollectMetrics>();
      std::sort(cell_indices.begin(), cell_indices.end());
      profile_add<CollectMetrics>(metrics.closure_sort_nanoseconds, sort_start);
      const auto materialize_start = profile_start<CollectMetrics>();
      for (std::size_t local_index : cell_indices) {
        if constexpr (CollectMetrics) {
          metrics.closure_entry_growths += cells.entries.size() == cells.entries.capacity();
        }
        cells.entries.push_back(bucket[local_index]);
        included[local_index] = 0;
      }
      cells.ranges[facet_index] = {first, cells.entries.size()};
      profile_add<CollectMetrics>(metrics.closure_materialize_nanoseconds, materialize_start);
    }
  }

  template <typename Function>
  std::size_t parallel_for_indices(std::size_t count,
                                   std::size_t task_count,
                                   Function&& function) const {
    const std::size_t chunk_size = count / task_count;
    const std::size_t extra = count % task_count;
    std::vector<std::future<void>> futures;
    futures.reserve(task_count);
    try {
      for (std::size_t task = 0; task < task_count; ++task) {
        const std::size_t first = task * chunk_size + std::min(task, extra);
        const std::size_t last = first + chunk_size + (task < extra ? 1 : 0);
        futures.push_back(executor_->submit([first, last, &function]() {
          for (std::size_t index = first; index < last; ++index) {
            function(index);
          }
        }));
      }
      for (auto& future : futures) {
        executor_->get(future);
      }
    } catch (...) {
      // A view accessor or submission may throw. Every task captures function
      // and its round buffers, so drain before those captures leave scope.
      for (auto& future : futures) {
        if (future.valid()) {
          try {
            executor_->get(future);
          } catch (...) {
          }
        }
      }
      throw;
    }
    return futures.size();
  }

  template <bool CollectMetrics>
  const std::vector<SimplexId>& active_facets(
      LevelId level, const std::vector<SimplexId>& bucket,
      std::size_t remaining,
      LevelScratch& scratch,
      ReductionKernelMetrics& metrics,
      bool allow_parallelism) const {
    auto& facets = scratch.facets;
    facets.clear();
    const auto& cells = scratch.level_cells;
    if (cells.packed_cofaces != nullptr) {
      // The active mask is updated only after all facet tasks finish. Its
      // set bits retain the original bucket order across every kernel round.
      for (std::size_t block = 0; block < cells.packed_block_count; ++block) {
        std::uint64_t candidates = cells.packed_active[block];
        while (candidates != 0) {
          const std::size_t index = 64 * block + trailing_zero_count(candidates);
          candidates &= candidates - 1;
          const auto* cofaces = cells.packed_cofaces->data() +
                                index * cells.packed_block_count;
          bool has_active_coface = false;
          for (std::size_t b = 0; b < cells.packed_block_count; ++b) {
            if constexpr (CollectMetrics) {
              ++metrics.facet_discovery_mask_tests;
            }
            if ((cofaces[b] & cells.packed_active[b]) != 0) {
              has_active_coface = true;
              break;
            }
          }
          if (!has_active_coface) {
            facets.push_back(bucket[index]);
          }
        }
      }
      return facets;
    }
    const std::size_t workers =
        executor_ == nullptr ? 1 : executor_->worker_count();
    const std::size_t task_count =
        std::min(workers, remaining / kMinFacetDiscoveryTaskSize);
    const bool sequential_scan =
        !allow_parallelism ||
        options_.policy == ReductionKernelExecutionPolicy::Sequential ||
        task_count <= 1;
    const bool use_cache = use_precomputed_cache();
    if (sequential_scan) {
      std::size_t active_count = 0;
      for (SimplexId simplex : scratch.active_simplices) {
        if (!active_[simplex]) {
          continue;
        }
        scratch.active_simplices[active_count++] = simplex;
        bool has_active_coface = false;
        visit_same_level_coboundary(
            simplex, level, use_cache, [&](SimplexId coface) {
              if constexpr (CollectMetrics) {
                ++metrics.facet_discovery_coboundary_visits;
              }
              if (active_[coface]) {
                has_active_coface = true;
                return false;
              }
              return true;
            });
        if (!has_active_coface) {
          facets.push_back(simplex);
        }
      }
      scratch.active_simplices.resize(active_count);
      return facets;
    }

    // Preserve canonical order while discarding removals before dispatch.
    // No task sees or mutates this list until all discovery tasks have joined.
    auto& candidates = scratch.active_simplices;
    candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
                                   [this](SimplexId simplex) {
                                     return !active_[simplex];
                                   }), candidates.end());
    auto& facet_flags = scratch.facet_flags;
    auto& coboundary_visits = scratch.coboundary_visits;
    const std::size_t parallel_tasks = parallel_for_indices(
        candidates.size(), task_count,
        [this, level, &candidates, &facet_flags,
         &coboundary_visits, use_cache](std::size_t index) {
          const SimplexId simplex = candidates[index];
          if constexpr (CollectMetrics) {
            coboundary_visits[index] = 0;
          }
          bool has_active_coface = false;
          visit_same_level_coboundary(
              simplex, level, use_cache, [&](SimplexId coface) {
                if constexpr (CollectMetrics) {
                  ++coboundary_visits[index];
                }
                if (active_[coface]) {
                  has_active_coface = true;
                  return false;
                }
                return true;
              });
          facet_flags[index] = !has_active_coface;
        });
    if constexpr (CollectMetrics) {
      metrics.facet_discovery_parallel_tasks += parallel_tasks;
    }

    for (std::size_t index = 0; index < candidates.size(); ++index) {
      if constexpr (CollectMetrics) {
        metrics.facet_discovery_coboundary_visits += coboundary_visits[index];
      }
      if (facet_flags[index]) {
        facets.push_back(candidates[index]);
      }
    }
    return facets;
  }

  template <bool CollectMetrics>
  void compute_facet_incidence(
      const std::vector<SimplexId>& facets,
      const std::vector<SimplexId>& bucket,
      LevelCells& level_cells,
      ReductionKernelMetrics& metrics) {
    if (level_cells.packed_masks != nullptr) {
      PackedMask seen{};
      PackedMask shared{};
      const std::size_t block_count = level_cells.packed_block_count;
      // A face is protected exactly when at least two current facets contain
      // it. Accumulate that predicate one word at a time instead of counting
      // incidence simplex by simplex. These masks are immutable during the
      // subsequent facet tasks, including intra-level parallel execution.
      for (SimplexId facet : facets) {
        const auto* cell_mask = level_cells.packed_masks->data() +
                                bucket_index_[facet] * block_count;
        for (std::size_t block = 0; block < block_count; ++block) {
          shared[block] |= seen[block] & cell_mask[block];
          seen[block] |= cell_mask[block];
        }
      }
      for (std::size_t block = 0; block < block_count; ++block) {
        level_cells.packed_unique[block] =
            seen[block] & ~shared[block] & level_cells.packed_active[block];
      }
      return;
    }
    // The coordinator completes incidence before facet tasks read it. Walking
    // cached closures is linear in their entries for every execution policy;
    // parallel simplex-versus-facet scans performed much more work on plateaus.
    if (level_cells.enabled) {
      for (SimplexId simplex : bucket) {
        facet_incidence_[simplex] = 0;
      }
      for (SimplexId facet : facets) {
        const auto [first, last] = cell_range(level_cells, facet);
        const auto& entries = cell_entries(level_cells);
        for (std::size_t index = first; index < last; ++index) {
          if constexpr (CollectMetrics) {
            ++metrics.incidence_cell_visits;
          }
          const SimplexId simplex = entries[index];
          if (active_[simplex] && facet_incidence_[simplex] < 2) {
            ++facet_incidence_[simplex];
          }
        }
      }
      return;
    }
    // Only dimension-0/1 levels omit closure storage. Their same-level closure
    // consists of the facet itself and its immediate boundary. Filter levels
    // before reading active_ so concurrently processed levels remain disjoint.
    for (SimplexId simplex : bucket) {
      facet_incidence_[simplex] = 0;
    }
    const auto visit = [&](SimplexId simplex) {
      if constexpr (CollectMetrics) {
        ++metrics.incidence_cell_visits;
      }
      if (active_[simplex] && facet_incidence_[simplex] < 2) {
        ++facet_incidence_[simplex];
      }
    };
    for (SimplexId facet : facets) {
      visit(facet);
      for (SimplexId face : complex_.boundary(facet)) {
        if (complex_.level(face) == complex_.level(facet)) {
          visit(face);
        }
      }
    }
  }

  template <bool CollectMetrics>
  FacetKernelResult<CollectMetrics> compute_packed_facet_kernel(
      SimplexId facet,
      const LevelCells& level_cells) const {
    FacetKernelResult<CollectMetrics> result;
    const auto& bucket = *level_cells.packed_bucket;
    const std::size_t block_count = level_cells.packed_block_count;
    const std::size_t facet_index = bucket_index_[facet];
    const auto* cell_mask =
        level_cells.packed_masks->data() + facet_index * block_count;

    if constexpr (CollectMetrics) {
      const auto core_start = profile_start<CollectMetrics>();
      for (std::size_t block = 0; block < block_count; ++block) {
        std::uint64_t entries = cell_mask[block];
        while (entries != 0) {
          ++result.facet_cell_visits;
          entries &= entries - 1;
        }
      }
      result.core_nanoseconds =
          elapsed_nanoseconds(core_start, Clock::now());
    }

    const auto reduction_start = profile_start<CollectMetrics>();
    PackedMask live_cell{};
    for (std::size_t block = 0; block < block_count; ++block) {
      live_cell[block] = cell_mask[block] & level_cells.packed_active[block];
    }
    while (true) {
      SimplexId reduction_sigma = kInvalidSimplex;
      SimplexId reduction_tau = kInvalidSimplex;
      std::size_t reduction_sigma_index = bucket.size();
      std::size_t reduction_tau_index = bucket.size();

      // Iterating set bits from low to high preserves the canonical bucket
      // order used by the sparse kernel without materializing a local cell.
      bool found_reduction = false;
      for (std::size_t block = 0;
           block < block_count && !found_reduction; ++block) {
        std::uint64_t candidates =
            live_cell[block] & level_cells.packed_unique[block];
        while (candidates != 0) {
          const std::size_t offset = trailing_zero_count(candidates);
          const std::size_t sigma_index = 64 * block + offset;
          candidates &= candidates - 1;
          if constexpr (CollectMetrics) {
            ++result.local_candidate_visits;
          }
          const SimplexId sigma = bucket[sigma_index];

          std::size_t unique_coface_index = bucket.size();
          const auto* cofaces = level_cells.packed_cofaces->data() +
                                sigma_index * block_count;
          for (std::size_t b = 0; b < block_count; ++b) {
            if constexpr (CollectMetrics) {
              ++result.local_coboundary_mask_tests;
            }
            const std::uint64_t live_cofaces = cofaces[b] & live_cell[b];
            if (live_cofaces == 0) {
              continue;
            }
            // Reject two bits in one word, or a second nonempty word. The
            // intersection includes protected cofaces in the uniqueness test.
            if ((live_cofaces & (live_cofaces - 1)) != 0 ||
                unique_coface_index != bucket.size()) {
              unique_coface_index = bucket.size();
              break;
            }
            unique_coface_index = 64 * b + trailing_zero_count(live_cofaces);
          }
          if (unique_coface_index != bucket.size() &&
              (level_cells.packed_unique[unique_coface_index / 64] &
               (std::uint64_t{1} << (unique_coface_index % 64))) != 0) {
            reduction_sigma = sigma;
            reduction_tau = bucket[unique_coface_index];
            reduction_sigma_index = sigma_index;
            reduction_tau_index = unique_coface_index;
            found_reduction = true;
            break;
          }
        }
      }

      if (reduction_sigma == kInvalidSimplex) {
        break;
      }
      live_cell[reduction_sigma_index / 64] &=
          ~(std::uint64_t{1} << (reduction_sigma_index % 64));
      live_cell[reduction_tau_index / 64] &=
          ~(std::uint64_t{1} << (reduction_tau_index % 64));
      result.events.push_back(ReductionKernelEvent{
          reduction_sigma, reduction_tau});
    }
    if constexpr (CollectMetrics) {
      result.local_reduction_nanoseconds =
          elapsed_nanoseconds(reduction_start, Clock::now());
      result.inline_event_overflows = result.events.uses_overflow() ? 1 : 0;
    }
    return result;
  }

  template <bool CollectMetrics>
  FacetKernelResult<CollectMetrics> compute_facet_kernel(
      LevelId level, SimplexId facet,
      const std::vector<SimplexId>& bucket,
      const LevelCells& level_cells) const {
    if (level_cells.packed_masks != nullptr) {
      return compute_packed_facet_kernel<CollectMetrics>(
          facet, level_cells);
    }
    FacetKernelResult<CollectMetrics> result;
    const auto core_start = profile_start<CollectMetrics>();
    InlineVector<SimplexId, kInlineCellCapacity> cell;
    if (level_cells.enabled) {
      const auto [first, last] = cell_range(level_cells, facet);
      const auto& entries = cell_entries(level_cells);
      for (std::size_t index = first; index < last; ++index) {
        if constexpr (CollectMetrics) {
          ++result.facet_cell_visits;
        }
        const SimplexId simplex = entries[index];
        if (active_[simplex]) {
          cell.push_back(simplex);
        }
      }
    } else {
      for (std::size_t index = 0; index < bucket.size(); ++index) {
        const SimplexId simplex = bucket[index];
        if constexpr (CollectMetrics) {
          ++result.facet_cell_visits;
        }
        if (active_[simplex] && is_face_of(simplex, facet)) {
          cell.push_back(simplex);
        }
      }
    }
    if constexpr (CollectMetrics) {
      result.core_nanoseconds =
          elapsed_nanoseconds(core_start, Clock::now());
      result.inline_cell_overflows = cell.uses_overflow() ? 1 : 0;
    }

    const auto reduction_start = profile_start<CollectMetrics>();
    const bool use_cache = use_precomputed_cache();
    std::array<std::uint8_t, kInlineCellCapacity> inline_removed{};
    std::vector<std::uint8_t> overflow_removed(
        cell.size() > inline_removed.size() ? cell.size() : 0, 0);
    auto is_locally_removed = [&](std::size_t index) {
      return overflow_removed.empty() ? inline_removed[index]
                                      : overflow_removed[index];
    };
    auto mark_locally_removed = [&](std::size_t index) {
      if (overflow_removed.empty()) {
        inline_removed[index] = 1;
      } else {
        overflow_removed[index] = 1;
      }
    };
    // Workspace cells inherit canonical bucket order, not numerical simplex
    // ID order. Preserve arbitrary external cache ordering via the old lookup.
    const bool ordered_lookup = cell.size() > kInlineCellCapacity &&
                                level_cells.cached_entries == nullptr;
    // Incidence is an immutable round snapshot: protected simplices cannot
    // become local reduction candidates. Filter only the scan, never the
    // full cell used for coface membership, and retain original cell indices.
    // Small and externally cached cells keep the allocation-free legacy scan.
    // Preparation (including any overflow allocation) is timed as local work.
    const bool filter_candidates = ordered_lookup;
    // Do not initialize the inline array on the unchanged small-cell path.
    std::optional<InlineVector<std::size_t, kInlineCellCapacity>> eligible_indices;
    if (filter_candidates) {
      eligible_indices.emplace();
      for (std::size_t index = 0; index < cell.size(); ++index) {
        if (facet_incidence_[cell[index]] <= 1) {
          eligible_indices->push_back(index);
        }
      }
    }
    const std::size_t candidate_count =
        filter_candidates ? eligible_indices->size() : cell.size();
    const auto find_coface_index = [&](SimplexId coface,
                                       std::size_t* comparisons = nullptr) {
      if (!ordered_lookup) {
        return cell.template index_of<CollectMetrics>(coface, comparisons);
      }
      const std::size_t rank = bucket_index_[coface];
      std::size_t first = 0;
      std::size_t last = cell.size();
      while (first < last) {
        const std::size_t middle = first + (last - first) / 2;
        if constexpr (CollectMetrics) {
          ++*comparisons;
        }
        if (bucket_index_[cell[middle]] < rank) {
          first = middle + 1;
        } else {
          last = middle;
        }
      }
      if (first < cell.size()) {
        if constexpr (CollectMetrics) {
          ++*comparisons;
        }
        if (cell[first] == coface) {
          return first;
        }
      }
      return cell.size();
    };
    while (true) {
      if constexpr (CollectMetrics) {
        ++result.local_sparse_scan_passes;
      }
      SimplexId reduction_sigma = kInvalidSimplex;
      SimplexId reduction_tau = kInvalidSimplex;
      std::size_t reduction_sigma_index = cell.size();
      std::size_t reduction_tau_index = cell.size();

      // Cell order is inherited from the level bucket, so local choices and
      // the merged event order are identical under both execution policies.
      for (std::size_t candidate = 0; candidate < candidate_count; ++candidate) {
        const std::size_t sigma_index =
            filter_candidates ? (*eligible_indices)[candidate] : candidate;
        if constexpr (CollectMetrics) {
          ++result.local_candidate_visits;
          ++result.local_sparse_candidate_visits;
        }
        const SimplexId sigma = cell[sigma_index];
        if (is_locally_removed(sigma_index) ||
            (!filter_candidates && facet_incidence_[sigma] > 1)) {
          if constexpr (CollectMetrics) {
            if (is_locally_removed(sigma_index)) {
              ++result.local_removed_candidate_visits;
            } else {
              ++result.local_protected_candidate_visits;
            }
          }
          continue;
        }
        SimplexId unique_coface = kInvalidSimplex;
        std::size_t unique_coface_index = cell.size();
        std::size_t coface_count = 0;
        visit_same_level_coboundary(
            sigma, level, use_cache, [&](SimplexId coface) {
              if constexpr (CollectMetrics) {
                ++result.local_coboundary_visits;
              }
              if (!active_[coface]) {
                return true;
              }
              if constexpr (CollectMetrics) {
                ++result.local_membership_tests;
              }
              std::size_t coface_index;
              if constexpr (CollectMetrics) {
                std::size_t comparisons = 0;
                coface_index = find_coface_index(coface, &comparisons);
                result.local_membership_comparisons += comparisons;
                if (cell.size() > kInlineCellCapacity) {
                  ++result.local_large_membership_tests;
                  result.local_large_membership_comparisons += comparisons;
                }
              } else {
                coface_index = find_coface_index(coface);
              }
              if (coface_index == cell.size()) {
                return true;
              }
              if (is_locally_removed(coface_index)) {
                return true;
              }
              unique_coface = coface;
              unique_coface_index = coface_index;
              ++coface_count;
              return coface_count <= 1;
            });
        if (coface_count == 1 && facet_incidence_[unique_coface] == 1) {
          reduction_sigma = sigma;
          reduction_tau = unique_coface;
          reduction_sigma_index = sigma_index;
          reduction_tau_index = unique_coface_index;
          break;
        }
      }

      if (reduction_sigma == kInvalidSimplex) {
        break;
      }
      mark_locally_removed(reduction_sigma_index);
      mark_locally_removed(reduction_tau_index);
      result.events.push_back(ReductionKernelEvent{
          reduction_sigma, reduction_tau});
    }
    if constexpr (CollectMetrics) {
      result.local_reduction_nanoseconds =
          elapsed_nanoseconds(reduction_start, Clock::now());
      result.inline_event_overflows = result.events.uses_overflow() ? 1 : 0;
    }
    return result;
  }

  std::size_t facet_task_count(
      const std::vector<SimplexId>& facets,
      const std::vector<SimplexId>& bucket,
      const LevelCells& cells, bool allow_parallelism) const {
    const std::size_t workers =
        executor_ == nullptr ? 1 : executor_->worker_count();
    if (!allow_parallelism ||
        options_.policy == ReductionKernelExecutionPolicy::Sequential ||
        workers <= 1 || facets.size() <= 1) {
      return 1;
    }
    const std::size_t max_tasks = std::min(workers, facets.size());
    std::size_t work = 0;
    for (SimplexId facet : facets) {
      std::size_t cost = 0;
      if (cells.packed_masks != nullptr) {
        const auto* mask = cells.packed_masks->data() +
                           bucket_index_[facet] * cells.packed_block_count;
        for (std::size_t block = 0; block < cells.packed_block_count; ++block) {
          cost += population_count(mask[block] & cells.packed_active[block]);
        }
      } else if (cells.enabled) {
        const auto [first, last] = cell_range(cells, facet);
        // The local kernel scans this full range even if some cells were
        // removed in earlier rounds. Do not walk it again for an exact count.
        cost = last - first;
      } else {
        cost = bucket.size();
      }
      // Saturate rather than overflowing for unusually large generic views.
      work += std::min(cost, std::numeric_limits<std::size_t>::max() - work);
      if (work / kMinFacetExecutionTaskWork >= max_tasks) {
        return max_tasks;
      }
    }
    return std::max<std::size_t>(1, work / kMinFacetExecutionTaskWork);
  }

  template <bool CollectMetrics>
  const std::vector<FacetKernelResult<CollectMetrics>>& execute_facets(
      LevelId level, const std::vector<SimplexId>& facets,
      const std::vector<SimplexId>& bucket,
      const LevelCells& level_cells,
      std::vector<FacetKernelResult<CollectMetrics>>& results,
      ReductionKernelMetrics& metrics,
      std::size_t task_count) const {
    results.clear();
    if (results.capacity() < facets.size()) {
      results.reserve(facets.size());
    }

    if (task_count <= 1) {
      for (SimplexId facet : facets) {
        results.push_back(
            compute_facet_kernel<CollectMetrics>(
                level, facet, bucket, level_cells));
      }
      return results;
    }

    // Allocate before dispatch: tasks write disjoint, stable slots, and the
    // coordinator consumes them in canonical facet order after all tasks join.
    results.resize(facets.size());
    // Roughly four chunks per worker balance unequal facet costs without a
    // future, queue lock, notification, and barrier for every individual facet.
    const std::size_t chunk_size =
        std::max<std::size_t>(1, facets.size() / task_count / 4);
    std::atomic<std::size_t> next_facet{0};
    std::vector<std::future<void>> futures;
    futures.reserve(task_count);
    if constexpr (CollectMetrics) {
      ++metrics.parallel_batches;
      metrics.facet_parallel_tasks += task_count;
      metrics.max_parallel_facets =
          std::max(metrics.max_parallel_facets, task_count);
    }
    try {
      for (std::size_t task = 0; task < task_count; ++task) {
        futures.push_back(executor_->submit([&, this, level]() {
          while (true) {
            const std::size_t first =
                next_facet.fetch_add(chunk_size, std::memory_order_relaxed);
            if (first >= facets.size()) {
              break;
            }
            const std::size_t last =
                first + std::min(chunk_size, facets.size() - first);
            for (std::size_t index = first; index < last; ++index) {
              results[index] = compute_facet_kernel<CollectMetrics>(
                  level, facets[index], bucket, level_cells);
            }
          }
        }));
      }
      for (auto& future : futures) {
        executor_->get(future);
      }
    } catch (...) {
      // Submission/allocation or a local kernel may throw. Drain every task
      // before captured state (especially next_facet) can leave scope, and
      // preserve the original exception if other tasks also fail.
      for (auto& future : futures) {
        if (future.valid()) {
          try {
            executor_->get(future);
          } catch (...) {
          }
        }
      }
      throw;
    }
    return results;
  }

  const ComplexView& complex_;
  ReductionKernelExecutionOptions options_;
  std::shared_ptr<BoundedTaskExecutor> executor_;
  std::vector<std::uint8_t> active_;
  // Only the categories zero, one, and more than one are required.
  std::vector<std::uint8_t> facet_incidence_;
  std::vector<std::uint8_t> round_removed_;
  std::vector<std::size_t> bucket_index_;
  std::vector<LevelScratch> level_scratch_;
  ReductionKernelMetrics metrics_;
};

}  // namespace morseframes
