#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
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
};

// Mutable scratch space for Algorithm 1 on immutable ComplexView data.  The
// sequential implementation deliberately exposes facet-local work as a
// separate phase so that a later implementation can schedule those cells in
// parallel without changing FSequenceBuilder or MorseSequence.
template <class ComplexView>
class ReductionKernelWorkspace {
  static_assert(is_complex_view_v<ComplexView>,
                "ReductionKernelWorkspace requires a Morse complex-view type.");

 private:
  using Clock = std::chrono::steady_clock;

  static std::uint64_t elapsed_nanoseconds(Clock::time_point start,
                                           Clock::time_point stop) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start)
            .count());
  }

 public:
  explicit ReductionKernelWorkspace(const ComplexView& complex)
      : complex_(complex),
        active_(complex.size(), 0),
        locally_removed_(complex.size(), 0),
        protected_core_(complex.size(), 0),
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
        std::vector<ReductionKernelEvent> round_events;

        for (SimplexId facet : facets) {
          std::vector<SimplexId> cell;
          for (SimplexId simplex : bucket) {
            if (active_[simplex] && is_face_of(simplex, facet)) {
              cell.push_back(simplex);
            }
          }

          std::fill(locally_removed_.begin(), locally_removed_.end(), 0);
          std::fill(protected_core_.begin(), protected_core_.end(), 0);

          const auto core_start = Clock::now();
          // Proposition 1 identifies the core with the attachment: a face of
          // this facet cell is protected exactly when another current facet
          // contains it.
          for (SimplexId simplex : cell) {
            for (SimplexId other_facet : facets) {
              if (other_facet != facet && is_face_of(simplex, other_facet)) {
                protected_core_[simplex] = 1;
                break;
              }
            }
          }
          metrics_.core_nanoseconds +=
              elapsed_nanoseconds(core_start, Clock::now());

          const auto reduction_start = Clock::now();
          while (true) {
            SimplexId reduction_sigma = kInvalidSimplex;
            SimplexId reduction_tau = kInvalidSimplex;

            // The bucket order is dimension/lexicographic, providing a
            // deterministic choice when several local kernels exist.
            for (SimplexId sigma : cell) {
              if (locally_removed_[sigma] || protected_core_[sigma]) {
                continue;
              }
              SimplexId unique_coface = kInvalidSimplex;
              std::size_t coface_count = 0;
              for (SimplexId coface : complex_.coboundary(sigma)) {
                if (!active_[coface] || locally_removed_[coface] ||
                    !is_face_of(coface, facet)) {
                  continue;
                }
                unique_coface = coface;
                ++coface_count;
                if (coface_count > 1) {
                  break;
                }
              }
              if (coface_count == 1 && !protected_core_[unique_coface]) {
                reduction_sigma = sigma;
                reduction_tau = unique_coface;
                break;
              }
            }

            if (reduction_sigma == kInvalidSimplex) {
              break;
            }
            locally_removed_[reduction_sigma] = 1;
            locally_removed_[reduction_tau] = 1;
            round_events.push_back(ReductionKernelEvent{
                ReductionKernelEventType::Reduction, reduction_sigma,
                reduction_tau});
            ++metrics_.reductions;
          }
          metrics_.local_reduction_nanoseconds +=
              elapsed_nanoseconds(reduction_start, Clock::now());

          for (SimplexId simplex : cell) {
            if (!locally_removed_[simplex]) {
              continue;
            }
            if (round_removed_[simplex]) {
              throw std::logic_error(
                  "Facet reduction kernels removed the same simplex twice.");
            }
            round_removed_[simplex] = 1;
          }
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

  const ComplexView& complex_;
  std::vector<std::uint8_t> active_;
  std::vector<std::uint8_t> locally_removed_;
  std::vector<std::uint8_t> protected_core_;
  std::vector<std::uint8_t> round_removed_;
  ReductionKernelMetrics metrics_;
};

}  // namespace morseframes
