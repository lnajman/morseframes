#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <future>
#include <stdexcept>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

#include "morseframes/complex_view.hpp"

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
  std::uint64_t core_nanoseconds = 0;
  std::uint64_t local_reduction_nanoseconds = 0;
  std::uint64_t merge_nanoseconds = 0;
  std::size_t levels = 0;
  std::size_t kernel_rounds = 0;
  std::size_t facet_kernels = 0;
  std::size_t reductions = 0;
  std::size_t perforations = 0;
  std::size_t parallel_batches = 0;
  std::size_t max_parallel_facets = 0;
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
      ReductionKernelExecutionOptions options = {})
      : complex_(complex),
        options_(options),
        active_(complex.size(), 0),
        round_removed_(complex.size(), 0) {}

  const ReductionKernelMetrics& metrics() const { return metrics_; }

  std::vector<ReductionKernelEvent> compute_level(LevelId level) {
    const auto& bucket = complex_.simplices_of_level(level);
    std::vector<ReductionKernelEvent> events;
    events.reserve(bucket.size());
    std::size_t remaining = bucket.size();
    ++metrics_.levels;

    for (SimplexId simplex : bucket) {
      active_[simplex] = 1;
    }

    while (remaining > 0) {
      bool kernel_round_changed = false;

      do {
        ++metrics_.kernel_rounds;
        kernel_round_changed = false;
        const auto facet_start = Clock::now();
        const auto facets = active_facets(level, bucket);
        metrics_.facet_nanoseconds +=
            elapsed_nanoseconds(facet_start, Clock::now());
        metrics_.facet_kernels += facets.size();
        std::fill(round_removed_.begin(), round_removed_.end(), 0);

        const auto facet_results = execute_facets(facets, bucket);
        std::vector<ReductionKernelEvent> round_events;
        for (const auto& result : facet_results) {
          metrics_.core_nanoseconds += result.core_nanoseconds;
          metrics_.local_reduction_nanoseconds +=
              result.local_reduction_nanoseconds;
          round_events.insert(round_events.end(), result.events.begin(),
                              result.events.end());
        }

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
        metrics_.reductions += round_events.size();
        metrics_.merge_nanoseconds +=
            elapsed_nanoseconds(merge_start, Clock::now());
      } while (kernel_round_changed);

      if (remaining == 0) {
        break;
      }

      const auto facet_start = Clock::now();
      const auto facets = active_facets(level, bucket);
      metrics_.facet_nanoseconds +=
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
      ++metrics_.perforations;
    }

    return events;
  }

 private:
  bool is_face_of(SimplexId face, SimplexId simplex) const {
    const auto& face_vertices = complex_.vertices(face);
    const auto& simplex_vertices = complex_.vertices(simplex);
    return std::includes(simplex_vertices.begin(), simplex_vertices.end(),
                         face_vertices.begin(), face_vertices.end());
  }

  std::vector<SimplexId> active_facets(
      LevelId level, const std::vector<SimplexId>& bucket) const {
    std::vector<SimplexId> facets;
    for (SimplexId simplex : bucket) {
      if (!active_[simplex]) {
        continue;
      }
      bool has_active_coface = false;
      for (SimplexId coface : complex_.coboundary(simplex)) {
        if (complex_.level(coface) == level && active_[coface]) {
          has_active_coface = true;
          break;
        }
      }
      if (!has_active_coface) {
        facets.push_back(simplex);
      }
    }
    return facets;
  }

  FacetKernelResult compute_facet_kernel(
      SimplexId facet, const std::vector<SimplexId>& facets,
      const std::vector<SimplexId>& bucket) const {
    FacetKernelResult result;
    std::vector<SimplexId> cell;
    for (SimplexId simplex : bucket) {
      if (active_[simplex] && is_face_of(simplex, facet)) {
        cell.push_back(simplex);
      }
    }

    const auto core_start = Clock::now();
    std::unordered_set<SimplexId> protected_core;
    for (SimplexId simplex : cell) {
      for (SimplexId other_facet : facets) {
        if (other_facet != facet && is_face_of(simplex, other_facet)) {
          protected_core.insert(simplex);
          break;
        }
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
            protected_core.count(sigma) != 0) {
          continue;
        }
        SimplexId unique_coface = kInvalidSimplex;
        std::size_t coface_count = 0;
        for (SimplexId coface : complex_.coboundary(sigma)) {
          if (!active_[coface] || locally_removed.count(coface) != 0 ||
              !is_face_of(coface, facet)) {
            continue;
          }
          unique_coface = coface;
          ++coface_count;
          if (coface_count > 1) {
            break;
          }
        }
        if (coface_count == 1 &&
            protected_core.count(unique_coface) == 0) {
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
      const std::vector<SimplexId>& facets,
      const std::vector<SimplexId>& bucket) {
    std::vector<FacetKernelResult> results;
    results.reserve(facets.size());

    std::size_t workers = options_.max_workers;
    if (workers == 0) {
      workers = std::max<std::size_t>(1, std::thread::hardware_concurrency());
    }
    if (options_.policy == ReductionKernelExecutionPolicy::Sequential ||
        workers <= 1 || facets.size() <= 1) {
      for (SimplexId facet : facets) {
        results.push_back(compute_facet_kernel(facet, facets, bucket));
      }
      return results;
    }

    workers = std::min(workers, facets.size());
    for (std::size_t first = 0; first < facets.size(); first += workers) {
      const std::size_t count = std::min(workers, facets.size() - first);
      ++metrics_.parallel_batches;
      metrics_.max_parallel_facets =
          std::max(metrics_.max_parallel_facets, count);
      std::vector<std::future<FacetKernelResult>> futures;
      futures.reserve(count);
      for (std::size_t offset = 0; offset < count; ++offset) {
        const SimplexId facet = facets[first + offset];
        futures.push_back(std::async(
            std::launch::async,
            [this, facet, &facets, &bucket]() {
              return compute_facet_kernel(facet, facets, bucket);
            }));
      }
      for (auto& future : futures) {
        results.push_back(future.get());
      }
    }
    return results;
  }

  const ComplexView& complex_;
  ReductionKernelExecutionOptions options_;
  std::vector<std::uint8_t> active_;
  std::vector<std::uint8_t> round_removed_;
  ReductionKernelMetrics metrics_;
};

}  // namespace morseframes
