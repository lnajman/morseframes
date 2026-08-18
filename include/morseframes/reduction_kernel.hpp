#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <future>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <unordered_set>
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
};

struct ReductionKernelMetrics {
  std::uint64_t facet_nanoseconds = 0;
  std::uint64_t essential_nanoseconds = 0;
  std::uint64_t core_nanoseconds = 0;
  std::uint64_t local_reduction_nanoseconds = 0;
  std::uint64_t aggregation_nanoseconds = 0;
  std::uint64_t merge_nanoseconds = 0;
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

  struct FacetKernelResult {
    std::vector<ReductionKernelEvent> events;
    std::uint64_t core_nanoseconds = 0;
    std::uint64_t local_reduction_nanoseconds = 0;
  };

  static std::uint64_t elapsed_nanoseconds(Clock::time_point start,
                                           Clock::time_point stop) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start)
            .count());
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
        round_removed_(complex.size(), 0) {
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

    for (SimplexId simplex : bucket) {
      active_[simplex] = 1;
    }

    while (remaining > 0) {
      bool kernel_round_changed = false;

      do {
        ++metrics.kernel_rounds;
        kernel_round_changed = false;
        const auto facet_start = Clock::now();
        const auto facets = active_facets(
            level, bucket, metrics, allow_intra_level_parallelism);
        metrics.facet_nanoseconds +=
            elapsed_nanoseconds(facet_start, Clock::now());
        metrics.facet_kernels += facets.size();
        const auto essential_start = Clock::now();
        compute_facet_incidence(
            facets, bucket, metrics, allow_intra_level_parallelism);
        metrics.essential_nanoseconds +=
            elapsed_nanoseconds(essential_start, Clock::now());
        for (SimplexId simplex : bucket) {
          round_removed_[simplex] = 0;
        }

        auto facet_results = execute_facets(
            level, facets, bucket, metrics, allow_intra_level_parallelism);
        const auto aggregation_start = Clock::now();
        auto round_result =
            aggregate_facet_results(
                std::move(facet_results), metrics,
                allow_intra_level_parallelism);
        metrics.aggregation_nanoseconds +=
            elapsed_nanoseconds(aggregation_start, Clock::now());
        metrics.core_nanoseconds += round_result.core_nanoseconds;
        metrics.local_reduction_nanoseconds +=
            round_result.local_reduction_nanoseconds;
        auto& round_events = round_result.events;

        for (const auto& event : round_events) {
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

        const auto merge_start = Clock::now();
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
        events.insert(events.end(), round_events.begin(), round_events.end());
        metrics.reductions += round_events.size();
        metrics.merge_nanoseconds +=
            elapsed_nanoseconds(merge_start, Clock::now());
      } while (kernel_round_changed);

      if (remaining == 0) {
        break;
      }

      const auto facet_start = Clock::now();
      const auto facets = active_facets(
          level, bucket, metrics, allow_intra_level_parallelism);
      metrics.facet_nanoseconds +=
          elapsed_nanoseconds(facet_start, Clock::now());
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
  }

 private:
  bool is_face_of(SimplexId face, SimplexId simplex) const {
    const auto& face_vertices = complex_.vertices(face);
    const auto& simplex_vertices = complex_.vertices(simplex);
    return std::includes(simplex_vertices.begin(), simplex_vertices.end(),
                         face_vertices.begin(), face_vertices.end());
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

  std::vector<SimplexId> active_facets(
      LevelId level, const std::vector<SimplexId>& bucket,
      ReductionKernelMetrics& metrics,
      bool allow_parallelism) const {
    std::vector<std::uint8_t> facet_flags(bucket.size(), 0);
    metrics.facet_discovery_parallel_tasks += parallel_for_indices(
        bucket.size(), [this, level, &bucket, &facet_flags](std::size_t index) {
          const SimplexId simplex = bucket[index];
          if (!active_[simplex]) {
            return;
          }
          for (SimplexId coface : complex_.coboundary(simplex)) {
            if (complex_.level(coface) == level && active_[coface]) {
              return;
            }
          }
          facet_flags[index] = 1;
        }, allow_parallelism);

    std::vector<SimplexId> facets;
    for (std::size_t index = 0; index < bucket.size(); ++index) {
      if (facet_flags[index]) {
        facets.push_back(bucket[index]);
      }
    }
    return facets;
  }

  void compute_facet_incidence(const std::vector<SimplexId>& facets,
                               const std::vector<SimplexId>& bucket,
                               ReductionKernelMetrics& metrics,
                               bool allow_parallelism) {
    metrics.essential_parallel_tasks += parallel_for_indices(
        bucket.size(), [this, &facets, &bucket](std::size_t index) {
          const SimplexId simplex = bucket[index];
          std::uint8_t incidence = 0;
          if (active_[simplex]) {
            for (SimplexId facet : facets) {
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
  }

  FacetKernelResult compute_facet_kernel(
      LevelId level, SimplexId facet,
      const std::vector<SimplexId>& bucket) const {
    FacetKernelResult result;
    const auto core_start = Clock::now();
    std::vector<SimplexId> cell;
    for (SimplexId simplex : bucket) {
      if (active_[simplex] && is_face_of(simplex, facet)) {
        cell.push_back(simplex);
      }
    }
    result.core_nanoseconds =
        elapsed_nanoseconds(core_start, Clock::now());

    const auto reduction_start = Clock::now();
    std::unordered_set<SimplexId> locally_removed;
    while (true) {
      SimplexId reduction_sigma = kInvalidSimplex;
      SimplexId reduction_tau = kInvalidSimplex;

      // Cell order is inherited from the level bucket, so local choices and
      // the merged event order are identical under both execution policies.
      for (SimplexId sigma : cell) {
        if (locally_removed.count(sigma) != 0 ||
            facet_incidence_[sigma] > 1) {
          continue;
        }
        SimplexId unique_coface = kInvalidSimplex;
        std::size_t coface_count = 0;
        for (SimplexId coface : complex_.coboundary(sigma)) {
          if (complex_.level(coface) != level || !active_[coface] ||
              locally_removed.count(coface) != 0 ||
              !is_face_of(coface, facet)) {
            continue;
          }
          unique_coface = coface;
          ++coface_count;
          if (coface_count > 1) {
            break;
          }
        }
        if (coface_count == 1 && facet_incidence_[unique_coface] == 1) {
          reduction_sigma = sigma;
          reduction_tau = unique_coface;
          break;
        }
      }

      if (reduction_sigma == kInvalidSimplex) {
        break;
      }
      locally_removed.insert(reduction_sigma);
      locally_removed.insert(reduction_tau);
      result.events.push_back(ReductionKernelEvent{
          ReductionKernelEventType::Reduction, reduction_sigma,
          reduction_tau});
    }
    result.local_reduction_nanoseconds =
        elapsed_nanoseconds(reduction_start, Clock::now());
    return result;
  }

  std::vector<FacetKernelResult> execute_facets(
      LevelId level, const std::vector<SimplexId>& facets,
      const std::vector<SimplexId>& bucket,
      ReductionKernelMetrics& metrics,
      bool allow_parallelism) const {
    std::vector<FacetKernelResult> results;
    results.reserve(facets.size());

    const std::size_t workers =
        executor_ == nullptr ? 1 : executor_->worker_count();
    if (!allow_parallelism ||
        options_.policy == ReductionKernelExecutionPolicy::Sequential ||
        workers <= 1 || facets.size() <= 1) {
      for (SimplexId facet : facets) {
        results.push_back(compute_facet_kernel(level, facet, bucket));
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
            [this, level, facet, &bucket]() {
              return compute_facet_kernel(level, facet, bucket);
            }));
      }
      for (auto& future : futures) {
        results.push_back(executor_->get(future));
      }
    }
    return results;
  }

  static FacetKernelResult combine_facet_results(FacetKernelResult left,
                                                  FacetKernelResult right) {
    left.events.reserve(left.events.size() + right.events.size());
    left.events.insert(left.events.end(),
                       std::make_move_iterator(right.events.begin()),
                       std::make_move_iterator(right.events.end()));
    left.core_nanoseconds += right.core_nanoseconds;
    left.local_reduction_nanoseconds +=
        right.local_reduction_nanoseconds;
    return left;
  }

  FacetKernelResult aggregate_facet_results(
      std::vector<FacetKernelResult> results,
      ReductionKernelMetrics& metrics,
      bool allow_parallelism) const {
    while (results.size() > 1) {
      ++metrics.aggregation_rounds;
      const std::size_t pair_count = results.size() / 2;
      std::vector<FacetKernelResult> next;
      next.reserve((results.size() + 1) / 2);

      const bool parallel =
          allow_parallelism &&
          options_.policy == ReductionKernelExecutionPolicy::Parallel &&
          executor_ != nullptr && executor_->worker_count() > 1 &&
          pair_count > 1;
      if (parallel) {
        std::vector<std::future<FacetKernelResult>> futures;
        futures.reserve(pair_count);
        for (std::size_t pair = 0; pair < pair_count; ++pair) {
          auto left = std::move(results[2 * pair]);
          auto right = std::move(results[2 * pair + 1]);
          futures.push_back(executor_->submit(
              [left = std::move(left), right = std::move(right)]() mutable {
                return combine_facet_results(std::move(left),
                                             std::move(right));
              }));
        }
        metrics.aggregation_parallel_tasks += futures.size();
        for (auto& future : futures) {
          next.push_back(executor_->get(future));
        }
      } else {
        for (std::size_t pair = 0; pair < pair_count; ++pair) {
          next.push_back(combine_facet_results(
              std::move(results[2 * pair]),
              std::move(results[2 * pair + 1])));
        }
      }

      if (results.size() % 2 != 0) {
        next.push_back(std::move(results.back()));
      }
      results = std::move(next);
    }
    return results.empty() ? FacetKernelResult{} : std::move(results.front());
  }

  const ComplexView& complex_;
  ReductionKernelExecutionOptions options_;
  std::shared_ptr<BoundedTaskExecutor> executor_;
  std::vector<std::uint8_t> active_;
  // Only the categories zero, one, and more than one are required.
  std::vector<std::uint8_t> facet_incidence_;
  std::vector<std::uint8_t> round_removed_;
  ReductionKernelMetrics metrics_;
};

}  // namespace morseframes
