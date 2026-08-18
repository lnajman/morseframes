#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <future>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "morseframes/complex_view.hpp"
#include "morseframes/task_executor.hpp"

namespace morseframes {

enum class ReductionKernelEventType {
  Reduction,
  Perforation,
};

struct ReductionKernelEvent {
  ReductionKernelEventType type = ReductionKernelEventType::Perforation;
  SimplexId sigma = kInvalidSimplex;
  SimplexId tau = kInvalidSimplex;
};

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

    void append_to(std::vector<T>& destination) const {
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

  struct LevelScratch {
    explicit LevelScratch(std::size_t bucket_size, bool collect_metrics)
        : facet_flags(bucket_size, 0),
          coboundary_visits(collect_metrics ? bucket_size : 0, 0) {
      facets.reserve(bucket_size);
    }

    std::vector<std::uint8_t> facet_flags;
    std::vector<std::size_t> coboundary_visits;
    std::vector<SimplexId> facets;
    std::vector<FacetKernelResult> facet_results;
  };

  struct LevelCells {
    bool enabled = false;
    std::vector<SimplexId> entries;
    std::vector<std::pair<std::size_t, std::size_t>> ranges;
  };

  static std::uint64_t elapsed_nanoseconds(Clock::time_point start,
                                           Clock::time_point stop) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start)
            .count());
  }

  Clock::time_point profile_start() const {
    return options_.collect_metrics ? Clock::now() : Clock::time_point{};
  }

  void profile_add(std::uint64_t& destination,
                   Clock::time_point start) const {
    if (options_.collect_metrics) {
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
    const auto& bucket = complex_.simplices_of_level(level);
    ReductionKernelLevelResult result;
    auto& events = result.events;
    auto& metrics = result.metrics;
    metrics.executor_workers =
        executor_ == nullptr ? 1 : executor_->worker_count();
    events.reserve(bucket.size());
    std::size_t remaining = bucket.size();
    ++metrics.levels;

    for (std::size_t index = 0; index < bucket.size(); ++index) {
      const SimplexId simplex = bucket[index];
      active_[simplex] = 1;
      bucket_index_[simplex] = index;
    }
    const auto closure_start = profile_start();
    const bool cache_level_cells =
        std::any_of(bucket.begin(), bucket.end(), [&](SimplexId simplex) {
          return complex_.dimension(simplex) >= 2;
        });
    const auto level_cells = build_level_cells(bucket, cache_level_cells);
    profile_add(metrics.closure_nanoseconds, closure_start);
    LevelScratch scratch(bucket.size(), options_.collect_metrics);

    while (remaining > 0) {
      bool kernel_round_changed = false;

      do {
        ++metrics.kernel_rounds;
        kernel_round_changed = false;
        const auto facet_start = profile_start();
        const auto& facets = active_facets(
            level, bucket, scratch, metrics,
            allow_intra_level_parallelism);
        profile_add(metrics.facet_nanoseconds, facet_start);
        metrics.facet_kernels += facets.size();
        const auto essential_start = profile_start();
        compute_facet_incidence(facets, bucket, level_cells, metrics,
                                allow_intra_level_parallelism);
        profile_add(metrics.essential_nanoseconds, essential_start);
        for (SimplexId simplex : bucket) {
          round_removed_[simplex] = 0;
        }

        const auto& facet_results =
            execute_facets(level, facets, bucket, level_cells,
                           scratch.facet_results, metrics,
                           allow_intra_level_parallelism);
        const auto aggregation_start = profile_start();
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
          metrics.inline_cell_overflows += facet_result.inline_cell_overflows;
          metrics.inline_event_overflows +=
              facet_result.inline_event_overflows;
          for (std::size_t index = 0; index < facet_result.events.size();
               ++index) {
            const auto& event = facet_result.events[index];
            if (event.type != ReductionKernelEventType::Reduction) {
              throw std::logic_error(
                  "A facet kernel returned a non-reduction event.");
            }
            if (round_removed_[event.sigma] || round_removed_[event.tau]) {
              throw std::logic_error(
                  "Facet reduction kernels removed the same simplex twice.");
            }
            round_removed_[event.sigma] = 1;
            round_removed_[event.tau] = 1;
          }
        }
        profile_add(metrics.aggregation_nanoseconds, aggregation_start);

        const auto merge_start = profile_start();
        for (SimplexId simplex : bucket) {
          if (!round_removed_[simplex]) {
            continue;
          }
          if (!active_[simplex]) {
            throw std::logic_error(
                "A reduction-kernel round removed an inactive simplex.");
          }
          active_[simplex] = 0;
          --remaining;
          kernel_round_changed = true;
        }
        for (const auto& facet_result : facet_results) {
          facet_result.events.append_to(events);
          metrics.reductions += facet_result.events.size();
        }
        profile_add(metrics.merge_nanoseconds, merge_start);
      } while (kernel_round_changed);

      if (remaining == 0) {
        break;
      }

      const auto facet_start = profile_start();
      const auto& facets = active_facets(
          level, bucket, scratch, metrics,
          allow_intra_level_parallelism);
      profile_add(metrics.facet_nanoseconds, facet_start);
      if (facets.empty()) {
        throw std::logic_error(
            "A nonempty reduction-kernel section has no facet.");
      }
      const SimplexId critical = facets.front();
      events.push_back(ReductionKernelEvent{
          ReductionKernelEventType::Perforation, critical, kInvalidSimplex});
      active_[critical] = 0;
      --remaining;
      ++metrics.perforations;
    }

    return result;
  }

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

  LevelCells build_level_cells(
      const std::vector<SimplexId>& bucket, bool enabled) const {
    LevelCells cells;
    cells.enabled = enabled;
    if (!enabled) {
      return cells;
    }
    cells.entries.reserve(4 * bucket.size());
    cells.ranges.reserve(bucket.size());
    std::vector<std::uint8_t> included(bucket.size(), 0);
    std::vector<std::size_t> cell_indices;
    cell_indices.reserve(16);
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
    return cells;
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

  const std::vector<SimplexId>& active_facets(
      LevelId level, const std::vector<SimplexId>& bucket,
      LevelScratch& scratch,
      ReductionKernelMetrics& metrics,
      bool allow_parallelism) const {
    auto& facet_flags = scratch.facet_flags;
    auto& coboundary_visits = scratch.coboundary_visits;
    std::fill(facet_flags.begin(), facet_flags.end(), 0);
    std::fill(coboundary_visits.begin(), coboundary_visits.end(), 0);
    metrics.facet_discovery_parallel_tasks += parallel_for_indices(
        bucket.size(), [this, level, &bucket, &facet_flags,
                        &coboundary_visits](std::size_t index) {
          const SimplexId simplex = bucket[index];
          if (!active_[simplex]) {
            return;
          }
          for (SimplexId coface : complex_.coboundary(simplex)) {
            if (options_.collect_metrics) {
              ++coboundary_visits[index];
            }
            if (complex_.level(coface) == level && active_[coface]) {
              return;
            }
          }
          facet_flags[index] = 1;
        }, allow_parallelism);
    for (std::size_t visits : coboundary_visits) {
      metrics.facet_discovery_coboundary_visits += visits;
    }

    auto& facets = scratch.facets;
    facets.clear();
    for (std::size_t index = 0; index < bucket.size(); ++index) {
      if (facet_flags[index]) {
        facets.push_back(bucket[index]);
      }
    }
    return facets;
  }

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
        const auto [first, last] =
            level_cells.ranges[bucket_index_[facet]];
        for (std::size_t index = first; index < last; ++index) {
          if (options_.collect_metrics) {
            ++metrics.incidence_cell_visits;
          }
          const SimplexId simplex = level_cells.entries[index];
          if (active_[simplex] && facet_incidence_[simplex] < 2) {
            ++facet_incidence_[simplex];
          }
        }
      }
      return;
    }
    std::vector<std::size_t> incidence_visits(
        options_.collect_metrics ? bucket.size() : 0, 0);
    metrics.essential_parallel_tasks += parallel_for_indices(
        bucket.size(), [this, &facets, &bucket,
                        &incidence_visits](std::size_t index) {
          const SimplexId simplex = bucket[index];
          std::uint8_t incidence = 0;
          if (active_[simplex]) {
            for (SimplexId facet : facets) {
              if (options_.collect_metrics) {
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
    for (std::size_t visits : incidence_visits) {
      metrics.incidence_cell_visits += visits;
    }
  }

  FacetKernelResult compute_facet_kernel(
      LevelId level, SimplexId facet,
      const std::vector<SimplexId>& bucket,
      const LevelCells& level_cells) const {
    FacetKernelResult result;
    const auto core_start = profile_start();
    InlineVector<SimplexId, kInlineCellCapacity> cell;
    if (level_cells.enabled) {
      const auto [first, last] = level_cells.ranges[bucket_index_[facet]];
      for (std::size_t index = first; index < last; ++index) {
        if (options_.collect_metrics) {
          ++result.facet_cell_visits;
        }
        const SimplexId simplex = level_cells.entries[index];
        if (active_[simplex]) {
          cell.push_back(simplex);
        }
      }
    } else {
      for (std::size_t index = 0; index < bucket.size(); ++index) {
        const SimplexId simplex = bucket[index];
        if (options_.collect_metrics) {
          ++result.facet_cell_visits;
        }
        if (active_[simplex] && is_face_of(simplex, facet)) {
          cell.push_back(simplex);
        }
      }
    }
    if (options_.collect_metrics) {
      result.core_nanoseconds =
          elapsed_nanoseconds(core_start, Clock::now());
      result.inline_cell_overflows = cell.uses_overflow() ? 1 : 0;
    }

    const auto reduction_start = profile_start();
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
        if (options_.collect_metrics) {
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
        for (SimplexId coface : complex_.coboundary(sigma)) {
          if (options_.collect_metrics) {
            ++result.local_coboundary_visits;
          }
          if (complex_.level(coface) != level || !active_[coface]) {
            continue;
          }
          if (options_.collect_metrics) {
            ++result.local_membership_tests;
          }
          const std::size_t coface_index = cell.index_of(coface);
          if (coface_index == cell.size()) {
            continue;
          }
          if (is_locally_removed(coface_index)) {
            continue;
          }
          unique_coface = coface;
          unique_coface_index = coface_index;
          ++coface_count;
          if (coface_count > 1) {
            break;
          }
        }
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
          ReductionKernelEventType::Reduction, reduction_sigma,
          reduction_tau});
    }
    if (options_.collect_metrics) {
      result.local_reduction_nanoseconds =
          elapsed_nanoseconds(reduction_start, Clock::now());
      result.inline_event_overflows = result.events.uses_overflow() ? 1 : 0;
    }
    return result;
  }

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
            compute_facet_kernel(level, facet, bucket, level_cells));
      }
      return results;
    }

    for (std::size_t first = 0; first < facets.size(); first += workers) {
      const std::size_t count = std::min(workers, facets.size() - first);
      ++metrics.parallel_batches;
      metrics.max_parallel_facets =
          std::max(metrics.max_parallel_facets, count);
      std::vector<std::future<FacetKernelResult>> futures;
      futures.reserve(count);
      for (std::size_t offset = 0; offset < count; ++offset) {
        const SimplexId facet = facets[first + offset];
        futures.push_back(executor_->submit(
            [this, level, facet, &bucket, &level_cells]() {
              return compute_facet_kernel(level, facet, bucket, level_cells);
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
  ReductionKernelMetrics metrics_;
};

}  // namespace morseframes
