#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <future>
#include <limits>
#include <memory>
#include <new>
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
  std::size_t levels = 0;
  std::size_t kernel_rounds = 0;
  std::size_t facet_kernels = 0;
  std::size_t reductions = 0;
  std::size_t perforations = 0;
  std::size_t parallel_batches = 0;
  std::size_t max_parallel_facets = 0;
  std::size_t parallel_level_batches = 0;
  std::size_t max_parallel_levels = 0;
  std::size_t executor_workers = 1;
  std::size_t facet_discovery_parallel_tasks = 0;
  std::size_t essential_parallel_tasks = 0;
  std::size_t aggregation_rounds = 0;
  std::size_t aggregation_parallel_tasks = 0;
  std::size_t facet_discovery_coboundary_visits = 0;
  std::size_t incidence_cell_visits = 0;
  std::size_t facet_cell_visits = 0;
  std::size_t local_candidate_visits = 0;
  std::size_t local_coboundary_visits = 0;
  std::size_t local_membership_tests = 0;
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

    std::size_t index_of(const T& value) const {
      for (std::size_t index = 0; index < size_; ++index) {
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

  struct FacetKernelResult {
    InlineVector<ReductionKernelEvent, kInlineEventCapacity> events;
    std::uint64_t core_nanoseconds = 0;
    std::uint64_t local_reduction_nanoseconds = 0;
    std::size_t facet_cell_visits = 0;
    std::size_t local_candidate_visits = 0;
    std::size_t local_coboundary_visits = 0;
    std::size_t local_membership_tests = 0;
    std::size_t inline_cell_overflows = 0;
    std::size_t inline_event_overflows = 0;
  };

  struct LevelCells {
    bool enabled = false;
    std::vector<SimplexId> entries;
    std::vector<std::pair<std::size_t, std::size_t>> ranges;
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
      round_events.clear();
      if (round_events.capacity() < bucket_size) {
        round_events.reserve(bucket_size);
      }
      included.resize(bucket_size);
      cell_indices.clear();
      if (cell_indices.capacity() < kInlineCellCapacity) {
        cell_indices.reserve(kInlineCellCapacity);
      }
      level_cells.enabled = false;
      level_cells.entries.clear();
      level_cells.ranges.clear();
      level_cells.cached_entries = nullptr;
      level_cells.cached_ranges = nullptr;
    }

    std::vector<std::uint8_t> facet_flags;
    std::vector<std::size_t> coboundary_visits;
    std::vector<SimplexId> facets;
    std::vector<SimplexId> active_simplices;
    std::vector<FacetKernelResult> facet_results;
    std::vector<ReductionKernelEvent> round_events;
    LevelCells level_cells;
    std::vector<std::uint8_t> included;
    std::vector<std::size_t> cell_indices;
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
    build_level_cells(bucket, cache_level_cells, scratch, level_cells);
    profile_add<CollectMetrics>(metrics.closure_nanoseconds, closure_start);

    while (remaining > 0) {
      bool kernel_round_changed = false;

      do {
        if constexpr (CollectMetrics) {
          ++metrics.kernel_rounds;
        }
        kernel_round_changed = false;
        const auto facet_start = profile_start<CollectMetrics>();
        const auto& facets = active_facets<CollectMetrics>(
            level, bucket, scratch, metrics,
            allow_intra_level_parallelism);
        profile_add<CollectMetrics>(metrics.facet_nanoseconds, facet_start);
        if constexpr (CollectMetrics) {
          metrics.facet_kernels += facets.size();
        }
        const auto essential_start = profile_start<CollectMetrics>();
        compute_facet_incidence<CollectMetrics>(
            facets, scratch.active_simplices, level_cells, metrics,
            allow_intra_level_parallelism);
        profile_add<CollectMetrics>(metrics.essential_nanoseconds,
                                    essential_start);
        scratch.round_events.clear();
        auto record_facet_events = [&](const FacetKernelResult& result) {
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

        if constexpr (!CollectMetrics) {
          const std::size_t workers =
              executor_ == nullptr ? 1 : executor_->worker_count();
          const bool execute_sequentially =
              !allow_intra_level_parallelism ||
              options_.policy == ReductionKernelExecutionPolicy::Sequential ||
              workers <= 1 || facets.size() <= 1;
          if (execute_sequentially) {
            for (SimplexId facet : facets) {
              const auto result = compute_facet_kernel<false>(
                  level, facet, scratch.active_simplices, level_cells);
              record_facet_events(result);
            }
          } else {
            const auto& facet_results = execute_facets<false>(
                level, facets, scratch.active_simplices, level_cells,
                scratch.facet_results, metrics,
                allow_intra_level_parallelism);
            for (const auto& result : facet_results) {
              record_facet_events(result);
            }
          }
        } else {
          const auto& facet_results = execute_facets<true>(
              level, facets, scratch.active_simplices, level_cells,
              scratch.facet_results, metrics,
              allow_intra_level_parallelism);
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
            metrics.local_membership_tests +=
                facet_result.local_membership_tests;
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
          level, bucket, scratch, metrics,
          allow_intra_level_parallelism);
      profile_add<CollectMetrics>(metrics.facet_nanoseconds, facet_start);
      if (facets.empty()) {
        throw std::logic_error(
            "A nonempty reduction-kernel section has no facet.");
      }
      const SimplexId critical = facets.front();
      events.push_back(ReductionKernelEvent{critical, kInvalidSimplex});
      active_[critical] = 0;
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
    destination.levels += source.levels;
    destination.kernel_rounds += source.kernel_rounds;
    destination.facet_kernels += source.facet_kernels;
    destination.reductions += source.reductions;
    destination.perforations += source.perforations;
    destination.parallel_batches += source.parallel_batches;
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
    destination.incidence_cell_visits += source.incidence_cell_visits;
    destination.facet_cell_visits += source.facet_cell_visits;
    destination.local_candidate_visits += source.local_candidate_visits;
    destination.local_coboundary_visits += source.local_coboundary_visits;
    destination.local_membership_tests += source.local_membership_tests;
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

  void build_level_cells(
      const std::vector<SimplexId>& bucket, bool enabled,
      LevelScratch& scratch, LevelCells& cells) const {
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
    if (cells.entries.capacity() < 4 * bucket.size()) {
      cells.entries.reserve(4 * bucket.size());
    }
    if (cells.ranges.capacity() < bucket.size()) {
      cells.ranges.reserve(bucket.size());
    }
    auto& included = scratch.included;
    std::fill(included.begin(), included.end(), 0);
    auto& cell_indices = scratch.cell_indices;
    for (std::size_t facet_index = 0; facet_index < bucket.size();
         ++facet_index) {
      const SimplexId facet = bucket[facet_index];
      const std::size_t first = cells.entries.size();
      cell_indices.clear();
      included[facet_index] = 1;
      cell_indices.push_back(facet_index);
      for (SimplexId face : complex_.boundary(facet)) {
        if (complex_.level(face) != complex_.level(facet)) {
          continue;
        }
        const std::size_t face_index = bucket_index_[face];
        if (face_index >= facet_index) {
          throw std::logic_error(
              "Reduction-kernel level bucket is not face-first.");
        }
        const auto [face_first, face_last] =
            cells.ranges[face_index];
        for (std::size_t index = face_first; index < face_last; ++index) {
          const std::size_t local_index =
              bucket_index_[cells.entries[index]];
          if (!included[local_index]) {
            included[local_index] = 1;
            cell_indices.push_back(local_index);
          }
        }
      }
      std::sort(cell_indices.begin(), cell_indices.end());
      for (std::size_t local_index : cell_indices) {
        cells.entries.push_back(bucket[local_index]);
        included[local_index] = 0;
      }
      cells.ranges.emplace_back(first, cells.entries.size());
    }
  }

  template <typename Function>
  std::size_t parallel_for_indices(std::size_t count,
                                   Function&& function,
                                   bool allow_parallelism) const {
    const std::size_t workers =
        executor_ == nullptr ? 1 : executor_->worker_count();
    if (!allow_parallelism ||
        options_.policy == ReductionKernelExecutionPolicy::Sequential ||
        workers <= 1 || count <= 1) {
      for (std::size_t index = 0; index < count; ++index) {
        function(index);
      }
      return 0;
    }

    const std::size_t task_count = std::min(workers, count);
    const std::size_t chunk_size = (count + task_count - 1) / task_count;
    std::vector<std::future<void>> futures;
    futures.reserve(task_count);
    for (std::size_t first = 0; first < count; first += chunk_size) {
      const std::size_t last = std::min(count, first + chunk_size);
      futures.push_back(executor_->submit([first, last, &function]() {
        for (std::size_t index = first; index < last; ++index) {
          function(index);
        }
      }));
    }
    for (auto& future : futures) {
      executor_->get(future);
    }
    return futures.size();
  }

  template <bool CollectMetrics>
  const std::vector<SimplexId>& active_facets(
      LevelId level, const std::vector<SimplexId>& bucket,
      LevelScratch& scratch,
      ReductionKernelMetrics& metrics,
      bool allow_parallelism) const {
    auto& facets = scratch.facets;
    facets.clear();
    const std::size_t workers =
        executor_ == nullptr ? 1 : executor_->worker_count();
    const bool sequential_scan =
        !allow_parallelism ||
        options_.policy == ReductionKernelExecutionPolicy::Sequential ||
        workers <= 1;
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

    auto& facet_flags = scratch.facet_flags;
    auto& coboundary_visits = scratch.coboundary_visits;
    std::fill(facet_flags.begin(), facet_flags.end(), 0);
    std::fill(coboundary_visits.begin(), coboundary_visits.end(), 0);
    const std::size_t parallel_tasks = parallel_for_indices(
        bucket.size(), [this, level, &bucket, &facet_flags,
                        &coboundary_visits, use_cache](std::size_t index) {
          const SimplexId simplex = bucket[index];
          if (!active_[simplex]) {
            return;
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
          if (has_active_coface) {
            return;
          }
          facet_flags[index] = 1;
        }, allow_parallelism);
    if constexpr (CollectMetrics) {
      metrics.facet_discovery_parallel_tasks += parallel_tasks;
      for (std::size_t visits : coboundary_visits) {
        metrics.facet_discovery_coboundary_visits += visits;
      }
    }

    for (std::size_t index = 0; index < bucket.size(); ++index) {
      if (facet_flags[index]) {
        facets.push_back(bucket[index]);
      }
    }
    return facets;
  }

  template <bool CollectMetrics>
  void compute_facet_incidence(
      const std::vector<SimplexId>& facets,
      const std::vector<SimplexId>& bucket,
      const LevelCells& level_cells,
      ReductionKernelMetrics& metrics,
      bool allow_parallelism) {
    const bool use_cached_cells =
        level_cells.enabled &&
        (options_.policy == ReductionKernelExecutionPolicy::Sequential ||
         !allow_parallelism || executor_ == nullptr ||
         executor_->worker_count() <= 1);
    if (use_cached_cells) {
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
    std::vector<std::size_t> incidence_visits(
        CollectMetrics ? bucket.size() : 0, 0);
    const std::size_t parallel_tasks = parallel_for_indices(
        bucket.size(), [this, &facets, &bucket,
                        &incidence_visits](std::size_t index) {
          const SimplexId simplex = bucket[index];
          std::uint8_t incidence = 0;
          if (active_[simplex]) {
            for (SimplexId facet : facets) {
              if constexpr (CollectMetrics) {
                ++incidence_visits[index];
              }
              if (is_face_of(simplex, facet)) {
                ++incidence;
                if (incidence == 2) {
                  break;
                }
              }
            }
          }
          facet_incidence_[simplex] = incidence;
        }, allow_parallelism);
    if constexpr (CollectMetrics) {
      metrics.essential_parallel_tasks += parallel_tasks;
      for (std::size_t visits : incidence_visits) {
        metrics.incidence_cell_visits += visits;
      }
    }
  }

  template <bool CollectMetrics>
  FacetKernelResult compute_facet_kernel(
      LevelId level, SimplexId facet,
      const std::vector<SimplexId>& bucket,
      const LevelCells& level_cells) const {
    FacetKernelResult result;
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
    while (true) {
      SimplexId reduction_sigma = kInvalidSimplex;
      SimplexId reduction_tau = kInvalidSimplex;
      std::size_t reduction_sigma_index = cell.size();
      std::size_t reduction_tau_index = cell.size();

      // Cell order is inherited from the level bucket, so local choices and
      // the merged event order are identical under both execution policies.
      for (std::size_t sigma_index = 0; sigma_index < cell.size();
           ++sigma_index) {
        if constexpr (CollectMetrics) {
          ++result.local_candidate_visits;
        }
        const SimplexId sigma = cell[sigma_index];
        if (is_locally_removed(sigma_index) ||
            facet_incidence_[sigma] > 1) {
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
              const std::size_t coface_index = cell.index_of(coface);
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

  template <bool CollectMetrics>
  const std::vector<FacetKernelResult>& execute_facets(
      LevelId level, const std::vector<SimplexId>& facets,
      const std::vector<SimplexId>& bucket,
      const LevelCells& level_cells,
      std::vector<FacetKernelResult>& results,
      ReductionKernelMetrics& metrics,
      bool allow_parallelism) const {
    results.clear();
    if (results.capacity() < facets.size()) {
      results.reserve(facets.size());
    }

    const std::size_t workers =
        executor_ == nullptr ? 1 : executor_->worker_count();
    if (!allow_parallelism ||
        options_.policy == ReductionKernelExecutionPolicy::Sequential ||
        workers <= 1 || facets.size() <= 1) {
      for (SimplexId facet : facets) {
        results.push_back(
            compute_facet_kernel<CollectMetrics>(
                level, facet, bucket, level_cells));
      }
      return results;
    }

    for (std::size_t first = 0; first < facets.size(); first += workers) {
      const std::size_t count = std::min(workers, facets.size() - first);
      if constexpr (CollectMetrics) {
        ++metrics.parallel_batches;
        metrics.max_parallel_facets =
            std::max(metrics.max_parallel_facets, count);
      }
      std::vector<std::future<FacetKernelResult>> futures;
      futures.reserve(count);
      for (std::size_t offset = 0; offset < count; ++offset) {
        const SimplexId facet = facets[first + offset];
        futures.push_back(executor_->submit(
            [this, level, facet, &bucket, &level_cells]() {
              return compute_facet_kernel<CollectMetrics>(
                  level, facet, bucket, level_cells);
            }));
      }
      for (auto& future : futures) {
        results.push_back(executor_->get(future));
      }
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
