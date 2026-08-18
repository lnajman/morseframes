#pragma once

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstddef>
#include <functional>
#include <future>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "morseframes/complex_view.hpp"
#include "morseframes/reduction_kernel.hpp"

#ifndef MORSE_ENABLE_SEQUENCE_BUILDER_CHECKS
#define MORSE_ENABLE_SEQUENCE_BUILDER_CHECKS 0
#endif

namespace morseframes {

inline constexpr bool kValidateSequenceBuilder =
    MORSE_ENABLE_SEQUENCE_BUILDER_CHECKS != 0;

enum class MorseStepType {
  Critical,
  RegularPair,
};

enum class FloodingScheme {
  Maximal,
  Minimal,
  MinMax,
  MaxMin,
};

struct MorseStep {
  MorseStepType type = MorseStepType::Critical;
  SimplexId sigma = kInvalidSimplex;
  SimplexId tau = kInvalidSimplex;
  LevelId level = 0;
};

struct MorseSequenceBuildMetrics {
  std::uint64_t init_nanoseconds = 0;
  std::uint64_t candidate_seed_nanoseconds = 0;
  std::uint64_t candidate_loop_nanoseconds = 0;
  std::uint64_t emit_nanoseconds = 0;
  std::uint64_t callback_nanoseconds = 0;
  std::uint64_t replay_nanoseconds = 0;
  std::uint64_t reduction_kernel_facet_nanoseconds = 0;
  std::uint64_t reduction_kernel_essential_nanoseconds = 0;
  std::uint64_t reduction_kernel_core_nanoseconds = 0;
  std::uint64_t reduction_kernel_local_reduction_nanoseconds = 0;
  std::uint64_t reduction_kernel_aggregation_nanoseconds = 0;
  std::uint64_t reduction_kernel_merge_nanoseconds = 0;
  std::uint64_t reduction_kernel_closure_nanoseconds = 0;
  std::uint64_t reduction_kernel_setup_nanoseconds = 0;
  std::uint64_t reduction_kernel_level_wall_nanoseconds = 0;
  std::uint64_t reduction_kernel_replay_nanoseconds = 0;
  std::uint64_t process_lower_stars_builder_init_nanoseconds = 0;
  std::uint64_t process_lower_stars_setup_nanoseconds = 0;
  std::uint64_t process_lower_stars_local_wall_nanoseconds = 0;
  std::uint64_t process_lower_stars_replay_nanoseconds = 0;
  std::uint64_t process_lower_stars_cumulative_task_nanoseconds = 0;
  std::uint64_t process_lower_stars_min_task_nanoseconds = 0;
  std::uint64_t process_lower_stars_max_task_nanoseconds = 0;
  std::size_t candidate_pushes = 0;
  std::size_t candidate_pops = 0;
  std::size_t stale_candidate_skips = 0;
  std::size_t level_mismatch_skips = 0;
  std::size_t regular_pairs = 0;
  std::size_t criticals = 0;
  std::size_t reduction_kernel_levels = 0;
  std::size_t reduction_kernel_rounds = 0;
  std::size_t reduction_kernel_facet_kernels = 0;
  std::size_t reduction_kernel_reductions = 0;
  std::size_t reduction_kernel_perforations = 0;
  std::size_t reduction_kernel_parallel_batches = 0;
  std::size_t reduction_kernel_max_parallel_facets = 0;
  std::size_t reduction_kernel_parallel_level_batches = 0;
  std::size_t reduction_kernel_max_parallel_levels = 0;
  std::size_t reduction_kernel_executor_workers = 1;
  std::size_t reduction_kernel_facet_discovery_parallel_tasks = 0;
  std::size_t reduction_kernel_essential_parallel_tasks = 0;
  std::size_t reduction_kernel_aggregation_rounds = 0;
  std::size_t reduction_kernel_aggregation_parallel_tasks = 0;
  std::size_t reduction_kernel_facet_discovery_coboundary_visits = 0;
  std::size_t reduction_kernel_incidence_cell_visits = 0;
  std::size_t reduction_kernel_facet_cell_visits = 0;
  std::size_t reduction_kernel_local_candidate_visits = 0;
  std::size_t reduction_kernel_local_coboundary_visits = 0;
  std::size_t reduction_kernel_local_membership_tests = 0;
  std::size_t reduction_kernel_inline_cell_overflows = 0;
  std::size_t reduction_kernel_inline_event_overflows = 0;
  std::size_t process_lower_stars_count = 0;
  std::size_t process_lower_stars_max_star_size = 0;
  std::size_t process_lower_stars_executor_workers = 1;
  std::size_t process_lower_stars_setup_parallel_tasks = 0;
  std::size_t process_lower_stars_parallel_tasks = 0;
  std::size_t process_lower_stars_min_task_load = 0;
  std::size_t process_lower_stars_max_task_load = 0;
};

class MorseSequence {
 public:
  explicit MorseSequence(std::size_t num_simplices)
      : num_simplices_(num_simplices) {
    steps_.reserve(num_simplices);
    critical_simplices_.reserve(num_simplices);
  }

  void add_critical(SimplexId sigma, LevelId level) {
    MorseStep step;
    step.type = MorseStepType::Critical;
    step.sigma = sigma;
    step.level = level;
    steps_.push_back(step);

    const auto critical_id = static_cast<std::int32_t>(critical_simplices_.size());
    critical_simplices_.push_back(sigma);
    if (!critical_index_of_simplex_.empty()) {
      critical_index_of_simplex_[sigma] = critical_id;
    }
  }

  void add_regular_pair(SimplexId sigma, SimplexId tau, LevelId level) {
    MorseStep step;
    step.type = MorseStepType::RegularPair;
    step.sigma = sigma;
    step.tau = tau;
    step.level = level;
    steps_.push_back(step);
  }

  const std::vector<MorseStep>& steps() const { return steps_; }
  const std::vector<SimplexId>& critical_simplices() const { return critical_simplices_; }
  const std::vector<std::int32_t>& critical_index_of_simplex() const {
    ensure_critical_index_map();
    return critical_index_of_simplex_;
  }

  std::int32_t critical_index(SimplexId simplex) const {
    ensure_critical_index_map();
    return critical_index_of_simplex_[simplex];
  }

  bool is_critical(SimplexId simplex) const { return critical_index(simplex) >= 0; }

 private:
  void ensure_critical_index_map() const {
    if (critical_index_of_simplex_.size() == num_simplices_) {
      return;
    }
    critical_index_of_simplex_.assign(num_simplices_, -1);
    for (std::size_t index = 0; index < critical_simplices_.size(); ++index) {
      critical_index_of_simplex_[critical_simplices_[index]] =
          static_cast<std::int32_t>(index);
    }
  }

  std::size_t num_simplices_ = 0;
  std::vector<MorseStep> steps_;
  std::vector<SimplexId> critical_simplices_;
  mutable std::vector<std::int32_t> critical_index_of_simplex_;
};

template <class ComplexView = FilteredSimplicialComplex>
class FSequenceBuilder {
  static_assert(is_complex_view_v<ComplexView>,
                "FSequenceBuilder requires a type satisfying the Morse complex-view API.");

 private:
  using SequenceClock = std::chrono::steady_clock;

  static constexpr std::size_t kInvalidSimplexRank =
      std::numeric_limits<std::size_t>::max();

  static std::uint64_t elapsed_nanoseconds(SequenceClock::time_point start,
                                           SequenceClock::time_point stop) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count());
  }

  SequenceClock::time_point profile_start() const {
    return sequence_metrics_ == nullptr ? SequenceClock::time_point{}
                                        : SequenceClock::now();
  }

  void profile_add(std::uint64_t MorseSequenceBuildMetrics::* target,
                   SequenceClock::time_point start) const {
    if (sequence_metrics_ != nullptr) {
      sequence_metrics_->*target += elapsed_nanoseconds(start, SequenceClock::now());
    }
  }

  std::size_t simplex_rank(SimplexId simplex) const {
    if (simplex >= simplex_order_rank_.size() ||
        simplex_order_rank_[simplex] == kInvalidSimplexRank) {
      throw std::logic_error("Simplex is missing from the filtration order.");
    }
    return simplex_order_rank_[simplex];
  }

  bool simplex_key_less(SimplexId lhs, SimplexId rhs) const {
    return simplex_rank(lhs) < simplex_rank(rhs);
  }

  std::uint16_t simplex_dimension(SimplexId simplex) const {
    if (simplex >= simplex_dimensions_.size()) {
      throw std::logic_error("Simplex is missing from the dimension table.");
    }
    return simplex_dimensions_[simplex];
  }

  LevelId simplex_level(SimplexId simplex) const {
    if (simplex >= simplex_levels_.size()) {
      throw std::logic_error("Simplex is missing from the level table.");
    }
    return simplex_levels_[simplex];
  }

  struct MinSimplexPriority {
    const FSequenceBuilder* builder = nullptr;

    bool operator()(SimplexId lhs, SimplexId rhs) const {
      return builder->simplex_order_rank_[rhs] < builder->simplex_order_rank_[lhs];
    }
  };

  struct MaxSimplexPriority {
    const FSequenceBuilder* builder = nullptr;

    bool operator()(SimplexId lhs, SimplexId rhs) const {
      return builder->simplex_order_rank_[lhs] < builder->simplex_order_rank_[rhs];
    }
  };

  struct MaxDimensionMinSimplexPriority {
    const FSequenceBuilder* builder = nullptr;

    bool operator()(SimplexId lhs, SimplexId rhs) const {
      const auto lhs_dimension = builder->simplex_dimensions_[lhs];
      const auto rhs_dimension = builder->simplex_dimensions_[rhs];
      if (lhs_dimension != rhs_dimension) {
        return lhs_dimension < rhs_dimension;
      }
      return builder->simplex_order_rank_[rhs] < builder->simplex_order_rank_[lhs];
    }
  };

  struct MinDimensionMinSimplexPriority {
    const FSequenceBuilder* builder = nullptr;

    bool operator()(SimplexId lhs, SimplexId rhs) const {
      const auto lhs_dimension = builder->simplex_dimensions_[lhs];
      const auto rhs_dimension = builder->simplex_dimensions_[rhs];
      if (lhs_dimension != rhs_dimension) {
        return lhs_dimension > rhs_dimension;
      }
      return builder->simplex_order_rank_[rhs] < builder->simplex_order_rank_[lhs];
    }
  };

  struct MaxDimensionMaxSimplexPriority {
    const FSequenceBuilder* builder = nullptr;

    bool operator()(SimplexId lhs, SimplexId rhs) const {
      const auto lhs_dimension = builder->simplex_dimensions_[lhs];
      const auto rhs_dimension = builder->simplex_dimensions_[rhs];
      if (lhs_dimension != rhs_dimension) {
        return lhs_dimension < rhs_dimension;
      }
      return builder->simplex_order_rank_[lhs] < builder->simplex_order_rank_[rhs];
    }
  };

  struct PlateauFillableCandidate {
    SimplexId simplex = kInvalidSimplex;
    std::size_t unlock_count = 0;
  };

  struct PlateauFillablePriority {
    const FSequenceBuilder* builder = nullptr;

    bool operator()(const PlateauFillableCandidate& lhs,
                    const PlateauFillableCandidate& rhs) const {
      if (lhs.unlock_count != rhs.unlock_count) {
        return lhs.unlock_count < rhs.unlock_count;
      }
      const auto lhs_dimension = builder->simplex_dimensions_[lhs.simplex];
      const auto rhs_dimension = builder->simplex_dimensions_[rhs.simplex];
      if (lhs_dimension != rhs_dimension) {
        return lhs_dimension > rhs_dimension;
      }
      return builder->simplex_order_rank_[rhs.simplex] <
             builder->simplex_order_rank_[lhs.simplex];
    }
  };

 public:
  explicit FSequenceBuilder(const ComplexView& complex,
                            MorseSequenceBuildMetrics* sequence_metrics = nullptr)
      : complex_(complex),
        simplex_order_rank_(complex.size(), kInvalidSimplexRank),
        simplex_levels_(complex.size(), 0),
        simplex_dimensions_(complex.size(), 0),
        sequence_metrics_(sequence_metrics) {
    const auto& order = complex_.filtration_order();
    if (order.size() != complex_.size()) {
      throw std::logic_error("Filtration order size does not match complex size.");
    }
    for (std::size_t rank = 0; rank < order.size(); ++rank) {
      const SimplexId simplex = order[rank];
      if (simplex >= simplex_order_rank_.size()) {
        throw std::logic_error("Filtration order contains an invalid simplex id.");
      }
      if (simplex_order_rank_[simplex] != kInvalidSimplexRank) {
        throw std::logic_error("Filtration order contains a duplicate simplex id.");
      }
      simplex_order_rank_[simplex] = rank;
      simplex_levels_[simplex] = complex_.level(simplex);
      simplex_dimensions_[simplex] = complex_.dimension(simplex);
    }
  }

  MorseSequence build_saturated() const {
    return build_saturated_with_step_callback([](const MorseSequence&, const MorseStep&) {});
  }

  MorseSequence build_plateau_greedy() const {
    return build_plateau_greedy_with_step_callback([](const MorseSequence&, const MorseStep&) {});
  }

  MorseSequence build_coreduction() const {
    return build_coreduction_with_step_callback([](const MorseSequence&, const MorseStep&) {});
  }

  MorseSequence build_same_level_reduction() const {
    return build_same_level_reduction_with_step_callback(
        [](const MorseSequence&, const MorseStep&) {});
  }

  // Forward Max(S,F)-style seed-and-expand construction.  This is an
  // F-sequence builder, but not necessarily a flooding sequence: after a
  // critical seed is inserted it consumes all available same-level
  // coreduction-like pairs before choosing the next seed.
  MorseSequence build_f_max() const {
    return build_f_max_with_step_callback([](const MorseSequence&, const MorseStep&) {});
  }

  // Simplicial ProcessLowerStars: partition by the unique maximal vertex and
  // run the forward one-missing-face expansion independently in each star.
  MorseSequence build_process_lower_stars() const {
    return build_process_lower_stars_with_step_callback(
        [](const MorseSequence&, const MorseStep&) {});
  }

  MorseSequence build_process_lower_stars_parallel(
      std::size_t max_workers = 0) const {
    return build_process_lower_stars_parallel_with_step_callback(
        [](const MorseSequence&, const MorseStep&) {}, max_workers);
  }

  // Decreasing Min(S,F)-style dual construction.  Events are removed from the
  // high end using same-level reduction-like pairs, then replayed in reverse to
  // produce the increasing Morse sequence consumed by the rest of the pipeline.
  MorseSequence build_f_min() const {
    return build_f_min_with_step_callback([](const MorseSequence&, const MorseStep&) {});
  }

  // The named flooding builders are filtration-monotone: they exhaust one
  // filtration value before moving to the next.
  MorseSequence build_flooding_max() const {
    return build_flooding_max_with_step_callback([](const MorseSequence&, const MorseStep&) {});
  }

  MorseSequence build_flooding_min() const {
    return build_flooding_min_with_step_callback([](const MorseSequence&, const MorseStep&) {});
  }

  MorseSequence build_flooding_reduction_kernel() const {
    return build_flooding_reduction_kernel_with_step_callback(
        [](const MorseSequence&, const MorseStep&) {});
  }

  MorseSequence build_flooding_reduction_kernel_parallel(
      std::size_t max_workers = 0) const {
    return build_flooding_reduction_kernel_parallel_with_step_callback(
        [](const MorseSequence&, const MorseStep&) {}, max_workers);
  }

  MorseSequence build_flooding_minmax() const {
    return build_flooding_minmax_with_step_callback([](const MorseSequence&, const MorseStep&) {});
  }

  MorseSequence build_flooding_maxmin() const {
    return build_flooding_maxmin_with_step_callback([](const MorseSequence&, const MorseStep&) {});
  }

  template <typename StepCallback>
  MorseSequence build_saturated_with_step_callback(StepCallback&& on_step) const {
    const std::size_t n = complex_.size();
    MorseSequence sequence(n);
    auto&& callback = on_step;

    std::vector<std::uint8_t> inserted(n, 0);
    std::vector<std::uint32_t> missing_count(n, 0);
    std::vector<SimplexId> missing_xor(n, 0);
    std::vector<std::uint32_t> remaining_by_level(complex_.num_levels(), 0);

    for (SimplexId simplex = 0; simplex < n; ++simplex) {
      missing_count[simplex] = static_cast<std::uint32_t>(complex_.boundary(simplex).size());
      for (SimplexId face : complex_.boundary(simplex)) {
        missing_xor[simplex] ^= face;
      }
      ++remaining_by_level[simplex_level(simplex)];
    }

    auto is_fillable = [&](SimplexId simplex) {
      return !inserted[simplex] && missing_count[simplex] == 0;
    };

    auto is_pairable = [&](SimplexId tau, LevelId level) {
      if (inserted[tau] || simplex_dimension(tau) == 0 || missing_count[tau] != 1) {
        return false;
      }
      const SimplexId sigma = missing_xor[tau];
      return sigma < n && !inserted[sigma] && simplex_level(sigma) == level;
    };

    for (LevelId level = 0; level < complex_.num_levels(); ++level) {
      const auto& bucket = complex_.simplices_of_level(level);
      std::priority_queue<SimplexId, std::vector<SimplexId>, MinSimplexPriority>
          pair_candidates(MinSimplexPriority{this});
      std::priority_queue<SimplexId, std::vector<SimplexId>, MinSimplexPriority>
          fillable_candidates(MinSimplexPriority{this});

      auto enqueue_current_level_candidate = [&](SimplexId simplex) {
        if (simplex_level(simplex) != level || inserted[simplex]) {
          return;
        }
        if (is_pairable(simplex, level)) {
          pair_candidates.push(simplex);
        }
        if (is_fillable(simplex)) {
          fillable_candidates.push(simplex);
        }
      };

      for (SimplexId simplex : bucket) {
        enqueue_current_level_candidate(simplex);
      }

      auto insert_simplex = [&](SimplexId simplex) {
        if (inserted[simplex]) {
          throw std::logic_error("Tried to insert a simplex twice.");
        }
        inserted[simplex] = 1;
        --remaining_by_level[simplex_level(simplex)];

        for (SimplexId coface : complex_.coboundary(simplex)) {
          if (inserted[coface]) {
            continue;
          }
          if (missing_count[coface] == 0) {
            throw std::logic_error("Missing-face count underflow.");
          }
          --missing_count[coface];
          missing_xor[coface] ^= simplex;
          enqueue_current_level_candidate(coface);
        }
      };

      while (remaining_by_level[level] > 0) {
        bool inserted_pair = false;

        while (!pair_candidates.empty()) {
          const SimplexId tau = pair_candidates.top();
          pair_candidates.pop();
          if (!is_pairable(tau, level)) {
            continue;
          }

          const SimplexId sigma = missing_xor[tau];
          sequence.add_regular_pair(sigma, tau, level);
          callback(sequence, sequence.steps().back());
          insert_simplex(sigma);
          insert_simplex(tau);
          inserted_pair = true;
          break;
        }

        if (inserted_pair) {
          continue;
        }

        SimplexId fillable = kInvalidSimplex;
        while (!fillable_candidates.empty()) {
          const SimplexId simplex = fillable_candidates.top();
          fillable_candidates.pop();
          if (is_fillable(simplex)) {
            fillable = simplex;
            break;
          }
        }

        if (fillable == kInvalidSimplex) {
          throw std::logic_error("No valid F-sequence step found.");
        }

        sequence.add_critical(fillable, level);
        callback(sequence, sequence.steps().back());
        insert_simplex(fillable);
      }
    }

    return sequence;
  }

  template <typename StepCallback>
  MorseSequence build_plateau_greedy_with_step_callback(StepCallback&& on_step) const {
    const std::size_t n = complex_.size();
    MorseSequence sequence(n);
    auto&& callback = on_step;

    std::vector<std::uint8_t> inserted(n, 0);
    std::vector<std::uint32_t> missing_count(n, 0);
    std::vector<SimplexId> missing_xor(n, 0);
    std::vector<std::uint32_t> remaining_by_level(complex_.num_levels(), 0);

    for (SimplexId simplex = 0; simplex < n; ++simplex) {
      missing_count[simplex] = static_cast<std::uint32_t>(complex_.boundary(simplex).size());
      for (SimplexId face : complex_.boundary(simplex)) {
        missing_xor[simplex] ^= face;
      }
      ++remaining_by_level[simplex_level(simplex)];
    }

    auto is_fillable = [&](SimplexId simplex) {
      return !inserted[simplex] && missing_count[simplex] == 0;
    };

    auto is_pairable = [&](SimplexId tau, LevelId level) {
      if (inserted[tau] || simplex_dimension(tau) == 0 || missing_count[tau] != 1) {
        return false;
      }
      const SimplexId sigma = missing_xor[tau];
      return sigma < n && !inserted[sigma] && simplex_level(sigma) == level;
    };

    for (LevelId level = 0; level < complex_.num_levels(); ++level) {
      const auto& bucket = complex_.simplices_of_level(level);
      std::priority_queue<SimplexId,
                          std::vector<SimplexId>,
                          MaxDimensionMinSimplexPriority>
          pair_candidates(MaxDimensionMinSimplexPriority{this});
      std::priority_queue<PlateauFillableCandidate,
                          std::vector<PlateauFillableCandidate>,
                          PlateauFillablePriority>
          fillable_candidates(PlateauFillablePriority{this});

      auto count_unlocks = [&](SimplexId simplex) {
        std::size_t unlocks = 0;
        for (SimplexId coface : complex_.coboundary(simplex)) {
          if (inserted[coface] || simplex_level(coface) != level || missing_count[coface] != 2) {
            continue;
          }
          const SimplexId other_missing = missing_xor[coface] ^ simplex;
          if (other_missing < n && !inserted[other_missing] &&
              simplex_level(other_missing) == level) {
            ++unlocks;
          }
        }
        return unlocks;
      };

      auto enqueue_pair_candidate = [&](SimplexId simplex) {
        if (simplex_level(simplex) != level || !is_pairable(simplex, level)) {
          return;
        }
        pair_candidates.push(simplex);
      };

      auto enqueue_fillable_candidate = [&](SimplexId simplex) {
        if (simplex_level(simplex) != level || !is_fillable(simplex)) {
          return;
        }
        fillable_candidates.push(PlateauFillableCandidate{simplex, count_unlocks(simplex)});
      };

      auto enqueue_current_level_candidate = [&](SimplexId simplex) {
        enqueue_pair_candidate(simplex);
        enqueue_fillable_candidate(simplex);
      };

      auto enqueue_missing_faces_for_unlocks = [&](SimplexId simplex) {
        if (simplex_level(simplex) != level || inserted[simplex] || missing_count[simplex] != 2) {
          return;
        }
        for (SimplexId face : complex_.boundary(simplex)) {
          if (!inserted[face]) {
            enqueue_fillable_candidate(face);
          }
        }
      };

      for (SimplexId simplex : bucket) {
        enqueue_current_level_candidate(simplex);
      }

      auto insert_simplex = [&](SimplexId simplex) {
        if (inserted[simplex]) {
          throw std::logic_error("Tried to insert a simplex twice.");
        }
        inserted[simplex] = 1;
        --remaining_by_level[simplex_level(simplex)];

        for (SimplexId coface : complex_.coboundary(simplex)) {
          if (inserted[coface]) {
            continue;
          }
          if (missing_count[coface] == 0) {
            throw std::logic_error("Missing-face count underflow.");
          }
          --missing_count[coface];
          missing_xor[coface] ^= simplex;
          enqueue_current_level_candidate(coface);
          enqueue_missing_faces_for_unlocks(coface);
        }
      };

      while (remaining_by_level[level] > 0) {
        bool inserted_pair = false;

        while (!pair_candidates.empty()) {
          const SimplexId tau = pair_candidates.top();
          pair_candidates.pop();
          if (!is_pairable(tau, level)) {
            continue;
          }

          const SimplexId sigma = missing_xor[tau];
          sequence.add_regular_pair(sigma, tau, level);
          callback(sequence, sequence.steps().back());
          insert_simplex(sigma);
          insert_simplex(tau);
          inserted_pair = true;
          break;
        }

        if (inserted_pair) {
          continue;
        }

        SimplexId fillable = kInvalidSimplex;
        while (!fillable_candidates.empty()) {
          const PlateauFillableCandidate candidate = fillable_candidates.top();
          fillable_candidates.pop();
          if (!is_fillable(candidate.simplex)) {
            continue;
          }
          const std::size_t current_unlocks = count_unlocks(candidate.simplex);
          if (current_unlocks != candidate.unlock_count) {
            fillable_candidates.push(
                PlateauFillableCandidate{candidate.simplex, current_unlocks});
            continue;
          }
          fillable = candidate.simplex;
          break;
        }

        if (fillable == kInvalidSimplex) {
          throw std::logic_error("No valid F-sequence step found.");
        }

        sequence.add_critical(fillable, level);
        callback(sequence, sequence.steps().back());
        insert_simplex(fillable);
      }
    }

    return sequence;
  }

  template <typename StepCallback>
  MorseSequence build_coreduction_with_step_callback(StepCallback&& on_step) const {
    return build_same_level_reduction_with_step_callback(
        std::forward<StepCallback>(on_step));
  }

  template <typename StepCallback>
  MorseSequence build_same_level_reduction_with_step_callback(StepCallback&& on_step) const {
    const std::size_t n = complex_.size();
    MorseSequence sequence(n);
    auto&& callback = on_step;

    struct CollapsePair {
      SimplexId sigma = kInvalidSimplex;
      SimplexId tau = kInvalidSimplex;
    };

    std::vector<std::uint8_t> inserted(n, 0);
    std::vector<std::uint8_t> active(n, 0);
    std::vector<std::uint32_t> active_coface_count(n, 0);

    auto emit_critical = [&](SimplexId sigma, LevelId level) {
      if (inserted[sigma]) {
        throw std::logic_error("Tried to insert a simplex twice.");
      }
      if constexpr (kValidateSequenceBuilder) {
        for (SimplexId face : complex_.boundary(sigma)) {
          if (!inserted[face]) {
            throw std::logic_error(
                "Same-level reduction critical has a missing boundary face.");
          }
        }
      }
      sequence.add_critical(sigma, level);
      callback(sequence, sequence.steps().back());
      inserted[sigma] = 1;
    };

    auto emit_pair = [&](SimplexId sigma, SimplexId tau, LevelId level) {
      if (inserted[sigma] || inserted[tau]) {
        throw std::logic_error("Tried to insert a regular pair twice.");
      }
      if constexpr (kValidateSequenceBuilder) {
        for (SimplexId face : complex_.boundary(tau)) {
          if (face != sigma && !inserted[face]) {
            throw std::logic_error("Same-level reduction pair has a missing boundary face.");
          }
        }
      }
      sequence.add_regular_pair(sigma, tau, level);
      callback(sequence, sequence.steps().back());
      inserted[sigma] = 1;
      inserted[tau] = 1;
    };

    for (LevelId level = 0; level < complex_.num_levels(); ++level) {
      const auto& bucket = complex_.simplices_of_level(level);
      std::priority_queue<SimplexId, std::vector<SimplexId>, MinSimplexPriority>
          free_faces(MinSimplexPriority{this});
      std::vector<CollapsePair> collapse_pairs;
      collapse_pairs.reserve(bucket.size() / 2);

      for (SimplexId simplex : bucket) {
        active[simplex] = 1;
        active_coface_count[simplex] = 0;
      }

      auto count_active_same_level_cofaces = [&](SimplexId simplex) {
        std::uint32_t count = 0;
        for (SimplexId coface : complex_.coboundary(simplex)) {
          if (simplex_level(coface) == level && active[coface]) {
            ++count;
          }
        }
        return count;
      };

      for (SimplexId simplex : bucket) {
        active_coface_count[simplex] = count_active_same_level_cofaces(simplex);
        if (active_coface_count[simplex] == 1) {
          free_faces.push(simplex);
        }
      }

      auto unique_active_same_level_coface = [&](SimplexId simplex) {
        SimplexId result = kInvalidSimplex;
        for (SimplexId coface : complex_.coboundary(simplex)) {
          if (simplex_level(coface) != level || !active[coface]) {
            continue;
          }
          if (result != kInvalidSimplex) {
            return kInvalidSimplex;
          }
          result = coface;
        }
        return result;
      };

      auto decrement_active_faces_of = [&](SimplexId simplex) {
        for (SimplexId face : complex_.boundary(simplex)) {
          if (simplex_level(face) != level || !active[face]) {
            continue;
          }
          if (active_coface_count[face] == 0) {
            throw std::logic_error("Active coface count underflow.");
          }
          --active_coface_count[face];
          if (active_coface_count[face] == 1) {
            free_faces.push(face);
          }
        }
      };

      while (!free_faces.empty()) {
        const SimplexId sigma = free_faces.top();
        free_faces.pop();
        if (!active[sigma] || active_coface_count[sigma] != 1) {
          continue;
        }

        const SimplexId tau = unique_active_same_level_coface(sigma);
        if (tau == kInvalidSimplex) {
          continue;
        }

        collapse_pairs.push_back(CollapsePair{sigma, tau});
        active[sigma] = 0;
        active[tau] = 0;
        decrement_active_faces_of(tau);
        decrement_active_faces_of(sigma);
      }

      for (SimplexId simplex : bucket) {
        if (!active[simplex]) {
          continue;
        }
        active[simplex] = 0;
        emit_critical(simplex, level);
      }

      for (std::size_t index = collapse_pairs.size(); index > 0; --index) {
        const auto& pair = collapse_pairs[index - 1];
        emit_pair(pair.sigma, pair.tau, level);
      }
    }

    return sequence;
  }

  template <typename StepCallback>
  MorseSequence build_f_max_with_step_callback(StepCallback&& on_step) const {
    const std::size_t n = complex_.size();
    const auto init_start = profile_start();
    MorseSequence sequence(n);
    auto&& callback = on_step;

    std::vector<std::uint8_t> inserted(n, 0);
    std::vector<std::uint32_t> remaining_boundary_count(n, 0);
    std::vector<SimplexId> remaining_boundary_xor(n, 0);
    std::priority_queue<SimplexId, std::vector<SimplexId>, MinSimplexPriority>
        coreduction_candidates(MinSimplexPriority{this});

    for (SimplexId simplex = 0; simplex < n; ++simplex) {
      const auto& boundary = complex_.boundary(simplex);
      remaining_boundary_count[simplex] = static_cast<std::uint32_t>(boundary.size());
      for (SimplexId face : boundary) {
        remaining_boundary_xor[simplex] ^= face;
      }
    }
    profile_add(&MorseSequenceBuildMetrics::init_nanoseconds, init_start);

    auto enqueue_coreduction_candidate = [&](SimplexId tau) {
      if (remaining_boundary_count[tau] != 1) {
        return;
      }
      const SimplexId sigma = remaining_boundary_xor[tau];
      if (sigma >= n || inserted[sigma]) {
        return;
      }
      if (simplex_level(sigma) != simplex_level(tau)) {
        if (sequence_metrics_ != nullptr) {
          ++sequence_metrics_->level_mismatch_skips;
        }
        return;
      }
      coreduction_candidates.push(tau);
      if (sequence_metrics_ != nullptr) {
        ++sequence_metrics_->candidate_pushes;
      }
    };

    const auto seed_start = profile_start();
    for (SimplexId simplex : complex_.filtration_order()) {
      enqueue_coreduction_candidate(simplex);
    }
    profile_add(&MorseSequenceBuildMetrics::candidate_seed_nanoseconds, seed_start);

    auto decrement_boundary_count = [&](SimplexId simplex) {
      for (SimplexId coface : complex_.coboundary(simplex)) {
        if (inserted[coface]) {
          continue;
        }
        if (remaining_boundary_count[coface] == 0) {
          throw std::logic_error("F-Max boundary count underflow.");
        }
        --remaining_boundary_count[coface];
        remaining_boundary_xor[coface] ^= simplex;
        enqueue_coreduction_candidate(coface);
      }
    };

    auto emit_critical = [&](SimplexId simplex) {
      const auto emit_start = profile_start();
      if (inserted[simplex]) {
        throw std::logic_error("Tried to insert an F-Max critical twice.");
      }
      if constexpr (kValidateSequenceBuilder) {
        for (SimplexId face : complex_.boundary(simplex)) {
          if (!inserted[face]) {
            throw std::logic_error("F-Max critical has a missing boundary face.");
          }
        }
      }
      sequence.add_critical(simplex, simplex_level(simplex));
      const auto callback_start = profile_start();
      profile_add(&MorseSequenceBuildMetrics::emit_nanoseconds, emit_start);
      callback(sequence, sequence.steps().back());
      profile_add(&MorseSequenceBuildMetrics::callback_nanoseconds, callback_start);
      const auto emit_resume = profile_start();
      inserted[simplex] = 1;
      decrement_boundary_count(simplex);
      if (sequence_metrics_ != nullptr) {
        ++sequence_metrics_->criticals;
      }
      profile_add(&MorseSequenceBuildMetrics::emit_nanoseconds, emit_resume);
    };

    auto emit_pair = [&](SimplexId sigma, SimplexId tau) {
      const auto emit_start = profile_start();
      if (inserted[sigma] || inserted[tau]) {
        throw std::logic_error("Tried to insert an F-Max pair twice.");
      }
      if (simplex_level(sigma) != simplex_level(tau)) {
        throw std::logic_error("F-Max pair crosses filtration levels.");
      }
      if constexpr (kValidateSequenceBuilder) {
        for (SimplexId face : complex_.boundary(tau)) {
          if (face != sigma && !inserted[face]) {
            throw std::logic_error("F-Max pair has a missing boundary face.");
          }
        }
      }
      sequence.add_regular_pair(sigma, tau, simplex_level(tau));
      const auto callback_start = profile_start();
      profile_add(&MorseSequenceBuildMetrics::emit_nanoseconds, emit_start);
      callback(sequence, sequence.steps().back());
      profile_add(&MorseSequenceBuildMetrics::callback_nanoseconds, callback_start);
      const auto emit_resume = profile_start();
      inserted[sigma] = 1;
      inserted[tau] = 1;
      decrement_boundary_count(sigma);
      decrement_boundary_count(tau);
      if (sequence_metrics_ != nullptr) {
        ++sequence_metrics_->regular_pairs;
      }
      profile_add(&MorseSequenceBuildMetrics::emit_nanoseconds, emit_resume);
    };

    std::size_t order_index = 0;
    const auto& order = complex_.filtration_order();
    while (order_index < order.size()) {
      while (!coreduction_candidates.empty()) {
        const auto candidate_start = profile_start();
        const SimplexId tau = coreduction_candidates.top();
        coreduction_candidates.pop();
        if (sequence_metrics_ != nullptr) {
          ++sequence_metrics_->candidate_pops;
        }
        if (inserted[tau] || remaining_boundary_count[tau] != 1) {
          if (sequence_metrics_ != nullptr) {
            ++sequence_metrics_->stale_candidate_skips;
          }
          profile_add(&MorseSequenceBuildMetrics::candidate_loop_nanoseconds,
                      candidate_start);
          continue;
        }
        const SimplexId sigma = remaining_boundary_xor[tau];
        if (sigma >= n || inserted[sigma]) {
          if (sequence_metrics_ != nullptr) {
            ++sequence_metrics_->stale_candidate_skips;
          }
          profile_add(&MorseSequenceBuildMetrics::candidate_loop_nanoseconds,
                      candidate_start);
          continue;
        }
        if (simplex_level(sigma) != simplex_level(tau)) {
          if (sequence_metrics_ != nullptr) {
            ++sequence_metrics_->level_mismatch_skips;
          }
          profile_add(&MorseSequenceBuildMetrics::candidate_loop_nanoseconds,
                      candidate_start);
          continue;
        }
        profile_add(&MorseSequenceBuildMetrics::candidate_loop_nanoseconds,
                    candidate_start);
        emit_pair(sigma, tau);
      }

      const auto scan_start = profile_start();
      while (order_index < order.size() && inserted[order[order_index]]) {
        ++order_index;
      }
      profile_add(&MorseSequenceBuildMetrics::candidate_loop_nanoseconds, scan_start);
      if (order_index < order.size()) {
        emit_critical(order[order_index]);
      }
    }

    return sequence;
  }

  template <typename StepCallback>
  MorseSequence build_process_lower_stars_with_step_callback(
      StepCallback&& on_step) const {
    return build_process_lower_stars_with_execution_options(
        std::forward<StepCallback>(on_step), 1);
  }

  template <typename StepCallback>
  MorseSequence build_process_lower_stars_parallel_with_step_callback(
      StepCallback&& on_step, std::size_t max_workers = 0) const {
    return build_process_lower_stars_with_execution_options(
        std::forward<StepCallback>(on_step), max_workers);
  }

  template <typename StepCallback>
  MorseSequence build_process_lower_stars_with_execution_options(
      StepCallback&& on_step, std::size_t max_workers) const {
    const auto setup_start = profile_start();
    const std::size_t n = complex_.size();
    MorseSequence sequence(n);
    auto&& callback = on_step;

    std::vector<SimplexId> vertex_order;
    vertex_order.reserve(n);
    std::unordered_map<VertexId, SimplexId> vertex_simplex;
    for (SimplexId simplex : complex_.filtration_order()) {
      if (simplex_dimension(simplex) != 0) {
        continue;
      }
      const auto& vertices = complex_.vertices(simplex);
      if (vertices.size() != 1) {
        throw std::invalid_argument(
            "ProcessLowerStars requires zero-cells with one vertex.");
      }
      if (!vertex_simplex.emplace(vertices[0], simplex).second) {
        throw std::invalid_argument(
            "ProcessLowerStars requires unique vertex identifiers.");
      }
      if (!vertex_order.empty() &&
          complex_.filtration(vertex_order.back()) ==
              complex_.filtration(simplex)) {
        throw std::invalid_argument(
            "ProcessLowerStars requires injective vertex filtration values.");
      }
      vertex_order.push_back(simplex);
    }
    if (vertex_order.empty()) {
      throw std::invalid_argument(
          "ProcessLowerStars requires at least one zero-cell.");
    }

    std::unordered_map<VertexId, std::size_t> vertex_rank;
    vertex_rank.reserve(vertex_order.size());
    for (std::size_t rank = 0; rank < vertex_order.size(); ++rank) {
      vertex_rank.emplace(complex_.vertices(vertex_order[rank])[0], rank);
    }

    BoundedTaskExecutor executor(max_workers);
    const std::size_t worker_count = executor.worker_count();
    std::vector<SimplexId> owner(n, kInvalidSimplex);
    std::vector<std::vector<SimplexId>> owned(n);
    std::vector<std::vector<std::size_t>> robins_key(n);
    auto build_owner_and_key = [&](SimplexId simplex) {
      const auto& vertices = complex_.vertices(simplex);
      if (vertices.empty()) {
        throw std::invalid_argument(
            "ProcessLowerStars does not support the empty simplex.");
      }
      auto& key = robins_key[simplex];
      key.reserve(vertices.size());
      SimplexId simplex_owner = kInvalidSimplex;
      std::size_t owner_rank = 0;
      for (VertexId vertex : vertices) {
        const auto rank_it = vertex_rank.find(vertex);
        const auto simplex_it = vertex_simplex.find(vertex);
        if (rank_it == vertex_rank.end() || simplex_it == vertex_simplex.end()) {
          throw std::invalid_argument(
              "ProcessLowerStars found a cell with an unknown vertex.");
        }
        key.push_back(rank_it->second);
        if (simplex_owner == kInvalidSimplex || rank_it->second > owner_rank) {
          simplex_owner = simplex_it->second;
          owner_rank = rank_it->second;
        }
      }
      std::sort(key.begin(), key.end(), std::greater<std::size_t>());
      if (complex_.level(simplex) != complex_.level(simplex_owner)) {
        throw std::invalid_argument(
            "ProcessLowerStars requires the max-vertex lower-star extension.");
      }
      owner[simplex] = simplex_owner;
    };

    constexpr std::size_t kParallelSetupThreshold = 512;
    if (worker_count > 1 && n >= kParallelSetupThreshold) {
      const std::size_t task_count = std::min(worker_count, n);
      const std::size_t chunk_size = (n + task_count - 1) / task_count;
      std::vector<std::future<void>> futures;
      futures.reserve(task_count);
      for (std::size_t first = 0; first < n; first += chunk_size) {
        const std::size_t last = std::min(n, first + chunk_size);
        futures.push_back(executor.submit([first, last, &build_owner_and_key]() {
          for (SimplexId simplex = first; simplex < last; ++simplex) {
            build_owner_and_key(simplex);
          }
        }));
      }
      if (sequence_metrics_ != nullptr) {
        sequence_metrics_->process_lower_stars_setup_parallel_tasks =
            futures.size();
      }
      for (auto& future : futures) {
        executor.get(future);
      }
    } else {
      for (SimplexId simplex = 0; simplex < n; ++simplex) {
        build_owner_and_key(simplex);
      }
    }
    for (SimplexId simplex = 0; simplex < n; ++simplex) {
      owned[owner[simplex]].push_back(simplex);
    }

    struct RobinsMinPriority {
      const std::vector<std::vector<std::size_t>>* keys = nullptr;

      bool operator()(SimplexId lhs, SimplexId rhs) const {
        const auto& lhs_key = (*keys)[lhs];
        const auto& rhs_key = (*keys)[rhs];
        if (lhs_key != rhs_key) {
          return std::lexicographical_compare(
              rhs_key.begin(), rhs_key.end(), lhs_key.begin(), lhs_key.end());
        }
        return rhs < lhs;
      }
    };

    struct LowerStarEvent {
      MorseStepType type = MorseStepType::Critical;
      SimplexId sigma = kInvalidSimplex;
      SimplexId tau = kInvalidSimplex;
    };
    std::vector<std::vector<LowerStarEvent>> events_by_star(vertex_order.size());

    auto process_lower_star = [&](std::size_t star_rank) {
      const SimplexId star_vertex = vertex_order[star_rank];
      const auto& lower_star = owned[star_vertex];
      auto& events = events_by_star[star_rank];
      events.reserve(lower_star.size());

      std::unordered_map<SimplexId, std::size_t> local_index;
      local_index.reserve(lower_star.size());
      for (std::size_t index = 0; index < lower_star.size(); ++index) {
        local_index.emplace(lower_star[index], index);
      }
      std::vector<std::uint8_t> classified(lower_star.size(), 0);
      std::vector<std::uint32_t> local_boundary_count(lower_star.size(), 0);
      std::vector<SimplexId> local_boundary_xor(lower_star.size(), 0);

      RobinsMinPriority priority{&robins_key};
      std::priority_queue<SimplexId, std::vector<SimplexId>, RobinsMinPriority>
          pair_candidates(priority);
      std::priority_queue<SimplexId, std::vector<SimplexId>, RobinsMinPriority>
          zero_candidates(priority);
      std::size_t remaining = lower_star.size();

      auto enqueue = [&](SimplexId simplex) {
        const std::size_t index = local_index.at(simplex);
        if (classified[index]) {
          return;
        }
        if (local_boundary_count[index] == 1) {
          pair_candidates.push(simplex);
        } else if (local_boundary_count[index] == 0) {
          zero_candidates.push(simplex);
        }
      };

      for (SimplexId simplex : lower_star) {
        const std::size_t index = local_index.at(simplex);
        for (SimplexId face : complex_.boundary(simplex)) {
          if (owner[face] == star_vertex) {
            ++local_boundary_count[index];
            local_boundary_xor[index] ^= face;
          }
        }
        enqueue(simplex);
      }

      auto mark_classified = [&](SimplexId simplex) {
        const std::size_t index = local_index.at(simplex);
        if (classified[index]) {
          throw std::logic_error(
              "ProcessLowerStars classified a simplex twice.");
        }
        classified[index] = 1;
        --remaining;
        for (SimplexId coface : complex_.coboundary(simplex)) {
          if (owner[coface] != star_vertex) {
            continue;
          }
          const std::size_t coface_index = local_index.at(coface);
          if (classified[coface_index]) {
            continue;
          }
          if (local_boundary_count[coface_index] == 0) {
            throw std::logic_error(
                "ProcessLowerStars local boundary count underflow.");
          }
          --local_boundary_count[coface_index];
          local_boundary_xor[coface_index] ^= simplex;
          enqueue(coface);
        }
      };

      while (remaining > 0) {
        bool paired = false;
        while (!pair_candidates.empty()) {
          const SimplexId tau = pair_candidates.top();
          pair_candidates.pop();
          const std::size_t tau_index = local_index.at(tau);
          if (classified[tau_index] || local_boundary_count[tau_index] != 1) {
            continue;
          }
          const SimplexId sigma = local_boundary_xor[tau_index];
          const auto sigma_it = local_index.find(sigma);
          if (sigma >= n || sigma_it == local_index.end() ||
              classified[sigma_it->second]) {
            continue;
          }
          events.push_back(
              LowerStarEvent{MorseStepType::RegularPair, sigma, tau});
          mark_classified(sigma);
          mark_classified(tau);
          paired = true;
          break;
        }
        if (paired) {
          continue;
        }

        SimplexId critical = kInvalidSimplex;
        while (!zero_candidates.empty()) {
          const SimplexId candidate = zero_candidates.top();
          zero_candidates.pop();
          const std::size_t candidate_index = local_index.at(candidate);
          if (!classified[candidate_index] &&
              local_boundary_count[candidate_index] == 0) {
            critical = candidate;
            break;
          }
        }
        if (critical == kInvalidSimplex) {
          throw std::logic_error(
              "ProcessLowerStars found no local expansion or critical step.");
        }
        events.push_back(
            LowerStarEvent{MorseStepType::Critical, critical, kInvalidSimplex});
        mark_classified(critical);
      }
    };

    if (sequence_metrics_ != nullptr) {
      sequence_metrics_->process_lower_stars_count = vertex_order.size();
      for (SimplexId star_vertex : vertex_order) {
        sequence_metrics_->process_lower_stars_max_star_size = std::max(
            sequence_metrics_->process_lower_stars_max_star_size,
            owned[star_vertex].size());
      }
      sequence_metrics_->process_lower_stars_executor_workers = worker_count;
    }
    profile_add(&MorseSequenceBuildMetrics::process_lower_stars_setup_nanoseconds,
                setup_start);
    const auto local_start = profile_start();
    if (worker_count <= 1 || vertex_order.size() <= 1) {
      for (std::size_t star_rank = 0; star_rank < vertex_order.size(); ++star_rank) {
        process_lower_star(star_rank);
      }
      if (sequence_metrics_ != nullptr) {
        sequence_metrics_->process_lower_stars_min_task_load = n;
        sequence_metrics_->process_lower_stars_max_task_load = n;
        const auto task_nanoseconds =
            elapsed_nanoseconds(local_start, SequenceClock::now());
        sequence_metrics_->process_lower_stars_cumulative_task_nanoseconds =
            task_nanoseconds;
        sequence_metrics_->process_lower_stars_min_task_nanoseconds =
            task_nanoseconds;
        sequence_metrics_->process_lower_stars_max_task_nanoseconds =
            task_nanoseconds;
      }
    } else {
      const std::size_t task_count = std::min(worker_count, vertex_order.size());
      std::vector<std::size_t> star_ranks(vertex_order.size());
      std::iota(star_ranks.begin(), star_ranks.end(), 0);
      std::sort(star_ranks.begin(), star_ranks.end(),
                [&](std::size_t lhs, std::size_t rhs) {
                  const std::size_t lhs_size = owned[vertex_order[lhs]].size();
                  const std::size_t rhs_size = owned[vertex_order[rhs]].size();
                  return lhs_size != rhs_size ? lhs_size > rhs_size : lhs < rhs;
                });

      std::vector<std::vector<std::size_t>> task_stars(task_count);
      std::vector<std::size_t> task_loads(task_count, 0);
      for (std::size_t star_rank : star_ranks) {
        const auto lightest =
            std::min_element(task_loads.begin(), task_loads.end());
        const std::size_t task_index =
            static_cast<std::size_t>(lightest - task_loads.begin());
        task_stars[task_index].push_back(star_rank);
        task_loads[task_index] += owned[vertex_order[star_rank]].size();
      }

      std::vector<std::future<void>> futures;
      futures.reserve(task_count);
      const bool measure_tasks = sequence_metrics_ != nullptr;
      std::vector<std::uint64_t> task_nanoseconds(
          measure_tasks ? task_count : 0, 0);
      for (std::size_t task_index = 0; task_index < task_count; ++task_index) {
        const auto& task = task_stars[task_index];
        futures.push_back(executor.submit([task, task_index, measure_tasks,
                                           &process_lower_star,
                                           &task_nanoseconds]() {
          const auto task_start =
              measure_tasks ? SequenceClock::now() : SequenceClock::time_point{};
          for (std::size_t star_rank : task) {
            process_lower_star(star_rank);
          }
          if (measure_tasks) {
            task_nanoseconds[task_index] =
                elapsed_nanoseconds(task_start, SequenceClock::now());
          }
        }));
      }
      if (sequence_metrics_ != nullptr) {
        sequence_metrics_->process_lower_stars_parallel_tasks = futures.size();
        sequence_metrics_->process_lower_stars_min_task_load =
            *std::min_element(task_loads.begin(), task_loads.end());
        sequence_metrics_->process_lower_stars_max_task_load =
            *std::max_element(task_loads.begin(), task_loads.end());
      }
      for (auto& future : futures) {
        executor.get(future);
      }
      if (sequence_metrics_ != nullptr) {
        sequence_metrics_->process_lower_stars_cumulative_task_nanoseconds =
            std::accumulate(task_nanoseconds.begin(), task_nanoseconds.end(),
                            std::uint64_t{0});
        sequence_metrics_->process_lower_stars_min_task_nanoseconds =
            *std::min_element(task_nanoseconds.begin(), task_nanoseconds.end());
        sequence_metrics_->process_lower_stars_max_task_nanoseconds =
            *std::max_element(task_nanoseconds.begin(), task_nanoseconds.end());
      }
    }
    profile_add(
        &MorseSequenceBuildMetrics::process_lower_stars_local_wall_nanoseconds,
        local_start);

    const auto replay_start = profile_start();
    for (const auto& events : events_by_star) {
      for (const LowerStarEvent& event : events) {
        if (event.type == MorseStepType::Critical) {
          sequence.add_critical(event.sigma, complex_.level(event.sigma));
          if (sequence_metrics_ != nullptr) {
            ++sequence_metrics_->criticals;
          }
        } else {
          sequence.add_regular_pair(event.sigma, event.tau,
                                    complex_.level(event.tau));
          if (sequence_metrics_ != nullptr) {
            ++sequence_metrics_->regular_pairs;
          }
        }
        callback(sequence, sequence.steps().back());
      }
    }
    profile_add(&MorseSequenceBuildMetrics::process_lower_stars_replay_nanoseconds,
                replay_start);

    return sequence;
  }

  template <typename StepCallback>
  MorseSequence build_f_min_with_step_callback(StepCallback&& on_step) const {
    const std::size_t n = complex_.size();
    const auto init_start = profile_start();
    MorseSequence sequence(n);
    auto&& callback = on_step;

    struct Event {
      MorseStepType type = MorseStepType::Critical;
      SimplexId sigma = kInvalidSimplex;
      SimplexId tau = kInvalidSimplex;
      LevelId level = 0;
    };

    std::vector<std::uint8_t> removed(n, 0);
    std::vector<std::uint32_t> remaining_coboundary_count(n, 0);
    std::vector<SimplexId> remaining_coboundary_xor(n, 0);
    std::priority_queue<SimplexId, std::vector<SimplexId>, MaxSimplexPriority>
        reduction_candidates(MaxSimplexPriority{this});
    std::vector<Event> decreasing_events;
    decreasing_events.reserve(n);

    for (SimplexId simplex = 0; simplex < n; ++simplex) {
      const auto& coboundary = complex_.coboundary(simplex);
      remaining_coboundary_count[simplex] =
          static_cast<std::uint32_t>(coboundary.size());
      for (SimplexId coface : coboundary) {
        remaining_coboundary_xor[simplex] ^= coface;
      }
    }
    profile_add(&MorseSequenceBuildMetrics::init_nanoseconds, init_start);

    const auto& increasing_order = complex_.filtration_order();
    auto enqueue_reduction_candidate = [&](SimplexId sigma) {
      if (remaining_coboundary_count[sigma] != 1) {
        return;
      }
      const SimplexId tau = remaining_coboundary_xor[sigma];
      if (tau >= n || removed[tau]) {
        return;
      }
      if (simplex_level(sigma) != simplex_level(tau)) {
        if (sequence_metrics_ != nullptr) {
          ++sequence_metrics_->level_mismatch_skips;
        }
        return;
      }
      reduction_candidates.push(sigma);
      if (sequence_metrics_ != nullptr) {
        ++sequence_metrics_->candidate_pushes;
      }
    };

    const auto seed_start = profile_start();
    for (std::size_t index = increasing_order.size(); index > 0; --index) {
      const SimplexId simplex = increasing_order[index - 1];
      enqueue_reduction_candidate(simplex);
    }
    profile_add(&MorseSequenceBuildMetrics::candidate_seed_nanoseconds, seed_start);

    auto decrement_coboundary_count = [&](SimplexId simplex) {
      for (SimplexId face : complex_.boundary(simplex)) {
        if (removed[face]) {
          continue;
        }
        if (remaining_coboundary_count[face] == 0) {
          throw std::logic_error("F-Min coboundary count underflow.");
        }
        --remaining_coboundary_count[face];
        remaining_coboundary_xor[face] ^= simplex;
        enqueue_reduction_candidate(face);
      }
    };

    auto remove_pair = [&](SimplexId sigma, SimplexId tau) {
      const auto emit_start = profile_start();
      if (removed[sigma] || removed[tau]) {
        throw std::logic_error("Tried to remove an F-Min pair twice.");
      }
      if (simplex_level(sigma) != simplex_level(tau)) {
        throw std::logic_error("F-Min pair crosses filtration levels.");
      }
      decreasing_events.push_back(
          Event{MorseStepType::RegularPair, sigma, tau, simplex_level(tau)});
      removed[sigma] = 1;
      removed[tau] = 1;
      decrement_coboundary_count(sigma);
      decrement_coboundary_count(tau);
      if (sequence_metrics_ != nullptr) {
        ++sequence_metrics_->regular_pairs;
      }
      profile_add(&MorseSequenceBuildMetrics::emit_nanoseconds, emit_start);
    };

    auto remove_critical = [&](SimplexId simplex) {
      const auto emit_start = profile_start();
      if (removed[simplex]) {
        throw std::logic_error("Tried to remove an F-Min critical twice.");
      }
      decreasing_events.push_back(
          Event{MorseStepType::Critical, simplex, kInvalidSimplex, simplex_level(simplex)});
      removed[simplex] = 1;
      decrement_coboundary_count(simplex);
      if (sequence_metrics_ != nullptr) {
        ++sequence_metrics_->criticals;
      }
      profile_add(&MorseSequenceBuildMetrics::emit_nanoseconds, emit_start);
    };

    std::size_t reverse_index = increasing_order.size();
    while (reverse_index > 0) {
      while (!reduction_candidates.empty()) {
        const auto candidate_start = profile_start();
        const SimplexId sigma = reduction_candidates.top();
        reduction_candidates.pop();
        if (sequence_metrics_ != nullptr) {
          ++sequence_metrics_->candidate_pops;
        }
        if (removed[sigma] || remaining_coboundary_count[sigma] != 1) {
          if (sequence_metrics_ != nullptr) {
            ++sequence_metrics_->stale_candidate_skips;
          }
          profile_add(&MorseSequenceBuildMetrics::candidate_loop_nanoseconds,
                      candidate_start);
          continue;
        }
        const SimplexId tau = remaining_coboundary_xor[sigma];
        if (tau >= n || removed[tau]) {
          if (sequence_metrics_ != nullptr) {
            ++sequence_metrics_->stale_candidate_skips;
          }
          profile_add(&MorseSequenceBuildMetrics::candidate_loop_nanoseconds,
                      candidate_start);
          continue;
        }
        if (simplex_level(sigma) != simplex_level(tau)) {
          if (sequence_metrics_ != nullptr) {
            ++sequence_metrics_->level_mismatch_skips;
          }
          profile_add(&MorseSequenceBuildMetrics::candidate_loop_nanoseconds,
                      candidate_start);
          continue;
        }
        profile_add(&MorseSequenceBuildMetrics::candidate_loop_nanoseconds,
                    candidate_start);
        remove_pair(sigma, tau);
      }

      const auto scan_start = profile_start();
      while (reverse_index > 0 && removed[increasing_order[reverse_index - 1]]) {
        --reverse_index;
      }
      profile_add(&MorseSequenceBuildMetrics::candidate_loop_nanoseconds, scan_start);
      if (reverse_index > 0) {
        remove_critical(increasing_order[reverse_index - 1]);
      }
    }

    auto replay_segment_start = profile_start();
    for (std::size_t index = decreasing_events.size(); index > 0; --index) {
      const Event& event = decreasing_events[index - 1];
      if (event.type == MorseStepType::Critical) {
        sequence.add_critical(event.sigma, event.level);
      } else {
        sequence.add_regular_pair(event.sigma, event.tau, event.level);
      }
      const auto callback_start = profile_start();
      profile_add(&MorseSequenceBuildMetrics::replay_nanoseconds, replay_segment_start);
      callback(sequence, sequence.steps().back());
      profile_add(&MorseSequenceBuildMetrics::callback_nanoseconds, callback_start);
      replay_segment_start = profile_start();
    }

    return sequence;
  }

  template <typename StepCallback>
  MorseSequence build_flooding_max_with_step_callback(StepCallback&& on_step) const {
    return build_flooding_with_step_callback(FloodingScheme::Maximal,
                                             std::forward<StepCallback>(on_step));
  }
  template <typename StepCallback>
  MorseSequence build_flooding_reduction_kernel_with_step_callback(
      StepCallback&& on_step) const {
    return build_flooding_reduction_kernel_with_execution_options(
        ReductionKernelExecutionOptions{},
        std::forward<StepCallback>(on_step));
  }

  template <typename StepCallback>
  MorseSequence build_flooding_reduction_kernel_parallel_with_step_callback(
      StepCallback&& on_step, std::size_t max_workers = 0) const {
    ReductionKernelExecutionOptions options;
    options.policy = ReductionKernelExecutionPolicy::Parallel;
    options.max_workers = max_workers;
    return build_flooding_reduction_kernel_with_execution_options(
        options, std::forward<StepCallback>(on_step));
  }

  template <typename StepCallback>
  MorseSequence build_flooding_reduction_kernel_with_execution_options(
      ReductionKernelExecutionOptions options, StepCallback&& on_step) const {
    options.collect_metrics = sequence_metrics_ != nullptr;
    const auto setup_start = profile_start();
    const std::size_t n = complex_.size();
    MorseSequence sequence(n);
    auto&& callback = on_step;
    const std::size_t num_levels = complex_.num_levels();
    std::shared_ptr<BoundedTaskExecutor> executor;
    if (options.policy == ReductionKernelExecutionPolicy::Parallel) {
      executor = std::make_shared<BoundedTaskExecutor>(options.max_workers);
    }
    const std::size_t workers =
        executor == nullptr ? 1 : executor->worker_count();
    const std::size_t level_workers = std::min(num_levels, workers);
    ReductionKernelWorkspace<ComplexView> workspace(complex_, options, executor);
    std::vector<std::size_t> event_offsets(num_levels + 1, 0);
    for (LevelId level = 0; level < num_levels; ++level) {
      event_offsets[level + 1] =
          event_offsets[level] + complex_.simplices_of_level(level).size();
    }
    const std::size_t event_capacity = event_offsets.back();
    if (event_capacity >
        std::numeric_limits<std::size_t>::max() /
            sizeof(ReductionKernelEvent)) {
      throw std::length_error("Reduction-kernel event arena is too large.");
    }
    std::unique_ptr<unsigned char[]> event_arena;
    if (event_capacity > 0) {
      event_arena.reset(
          new unsigned char[event_capacity * sizeof(ReductionKernelEvent)]);
    }
    auto* level_events =
        reinterpret_cast<ReductionKernelEvent*>(event_arena.get());
    std::vector<std::size_t> level_event_counts(num_levels, 0);
    std::vector<ReductionKernelMetrics> level_metrics(
        options.collect_metrics ? num_levels : 0);
    ReductionKernelMetrics kernel_metrics;
    kernel_metrics.executor_workers = workers;
    profile_add(&MorseSequenceBuildMetrics::reduction_kernel_setup_nanoseconds,
                setup_start);

    const auto level_start = profile_start();
    if (level_workers == 1 || num_levels <= 1) {
      for (LevelId level = 0; level < num_levels; ++level) {
        if (options.collect_metrics) {
          level_metrics[level] = workspace.compute_level_isolated_into(
              level, 0, level_events + event_offsets[level],
              event_offsets[level + 1] - event_offsets[level],
              level_event_counts[level]);
        } else {
          workspace.compute_level_isolated_into_unprofiled(
              level, 0, level_events + event_offsets[level],
              event_offsets[level + 1] - event_offsets[level],
              level_event_counts[level]);
        }
      }
    } else {
      // Levels own disjoint workspace entries. Long-lived tasks dynamically
      // claim levels so the executor sees only one task per worker and balance
      // follows actual kernel cost rather than a simplex-count proxy. Nested
      // facet tasks are disabled on this path; a single large plateau still
      // uses the intra-level parallel algorithm.
      const std::size_t task_count = std::min(level_workers, num_levels);
      std::atomic<LevelId> next_level{0};
      if (options.collect_metrics) {
        ++kernel_metrics.parallel_level_batches;
        kernel_metrics.max_parallel_levels = task_count;
      }
      std::vector<std::future<void>> futures;
      futures.reserve(task_count);
      for (std::size_t task = 0; task < task_count; ++task) {
        futures.push_back(executor->submit(
            [task, &next_level, num_levels, &workspace, &event_offsets,
             level_events, &level_event_counts, &level_metrics,
             collect_metrics = options.collect_metrics]() {
              while (true) {
                const LevelId level = next_level.fetch_add(
                    1, std::memory_order_relaxed);
                if (level >= num_levels) {
                  return;
                }
                if (collect_metrics) {
                  level_metrics[level] = workspace.compute_level_isolated_into(
                      level, task, level_events + event_offsets[level],
                      event_offsets[level + 1] - event_offsets[level],
                      level_event_counts[level], false);
                } else {
                  workspace.compute_level_isolated_into_unprofiled(
                      level, task, level_events + event_offsets[level],
                      event_offsets[level + 1] - event_offsets[level],
                      level_event_counts[level], false);
                }
              }
            }));
      }
      for (auto& future : futures) {
        executor->get(future);
      }
    }
    profile_add(
        &MorseSequenceBuildMetrics::reduction_kernel_level_wall_nanoseconds,
        level_start);

    const auto replay_start = profile_start();
    for (LevelId level = 0; level < num_levels; ++level) {
      if (options.collect_metrics) {
        ReductionKernelWorkspace<ComplexView>::accumulate_metrics(
            kernel_metrics, level_metrics[level]);
      }
      const std::size_t first = event_offsets[level];
      for (std::size_t index = level_event_counts[level]; index > 0; --index) {
        const auto& event = level_events[first + index - 1];
        if (event.is_perforation()) {
          sequence.add_critical(event.sigma, level);
          if (sequence_metrics_ != nullptr) {
            ++sequence_metrics_->criticals;
          }
        } else {
          sequence.add_regular_pair(event.sigma, event.tau, level);
          if (sequence_metrics_ != nullptr) {
            ++sequence_metrics_->regular_pairs;
          }
        }
        callback(sequence, sequence.steps().back());
      }
    }
    profile_add(
        &MorseSequenceBuildMetrics::reduction_kernel_replay_nanoseconds,
        replay_start);

    if (sequence_metrics_ != nullptr) {
      sequence_metrics_->reduction_kernel_facet_nanoseconds =
          kernel_metrics.facet_nanoseconds;
      sequence_metrics_->reduction_kernel_essential_nanoseconds =
          kernel_metrics.essential_nanoseconds;
      sequence_metrics_->reduction_kernel_core_nanoseconds =
          kernel_metrics.core_nanoseconds;
      sequence_metrics_->reduction_kernel_local_reduction_nanoseconds =
          kernel_metrics.local_reduction_nanoseconds;
      sequence_metrics_->reduction_kernel_aggregation_nanoseconds =
          kernel_metrics.aggregation_nanoseconds;
      sequence_metrics_->reduction_kernel_merge_nanoseconds =
          kernel_metrics.merge_nanoseconds;
      sequence_metrics_->reduction_kernel_closure_nanoseconds =
          kernel_metrics.closure_nanoseconds;
      sequence_metrics_->reduction_kernel_levels = kernel_metrics.levels;
      sequence_metrics_->reduction_kernel_rounds = kernel_metrics.kernel_rounds;
      sequence_metrics_->reduction_kernel_facet_kernels =
          kernel_metrics.facet_kernels;
      sequence_metrics_->reduction_kernel_reductions =
          kernel_metrics.reductions;
      sequence_metrics_->reduction_kernel_perforations =
          kernel_metrics.perforations;
      sequence_metrics_->reduction_kernel_parallel_batches =
          kernel_metrics.parallel_batches;
      sequence_metrics_->reduction_kernel_max_parallel_facets =
          kernel_metrics.max_parallel_facets;
      sequence_metrics_->reduction_kernel_parallel_level_batches =
          kernel_metrics.parallel_level_batches;
      sequence_metrics_->reduction_kernel_max_parallel_levels =
          kernel_metrics.max_parallel_levels;
      sequence_metrics_->reduction_kernel_executor_workers =
          kernel_metrics.executor_workers;
      sequence_metrics_->reduction_kernel_facet_discovery_parallel_tasks =
          kernel_metrics.facet_discovery_parallel_tasks;
      sequence_metrics_->reduction_kernel_essential_parallel_tasks =
          kernel_metrics.essential_parallel_tasks;
      sequence_metrics_->reduction_kernel_aggregation_rounds =
          kernel_metrics.aggregation_rounds;
      sequence_metrics_->reduction_kernel_aggregation_parallel_tasks =
          kernel_metrics.aggregation_parallel_tasks;
      sequence_metrics_->reduction_kernel_facet_discovery_coboundary_visits =
          kernel_metrics.facet_discovery_coboundary_visits;
      sequence_metrics_->reduction_kernel_incidence_cell_visits =
          kernel_metrics.incidence_cell_visits;
      sequence_metrics_->reduction_kernel_facet_cell_visits =
          kernel_metrics.facet_cell_visits;
      sequence_metrics_->reduction_kernel_local_candidate_visits =
          kernel_metrics.local_candidate_visits;
      sequence_metrics_->reduction_kernel_local_coboundary_visits =
          kernel_metrics.local_coboundary_visits;
      sequence_metrics_->reduction_kernel_local_membership_tests =
          kernel_metrics.local_membership_tests;
      sequence_metrics_->reduction_kernel_inline_cell_overflows =
          kernel_metrics.inline_cell_overflows;
      sequence_metrics_->reduction_kernel_inline_event_overflows =
          kernel_metrics.inline_event_overflows;
    }

    return sequence;
  }

  template <typename StepCallback>
  MorseSequence build_flooding_min_with_step_callback(StepCallback&& on_step) const {
    return build_flooding_with_step_callback(FloodingScheme::Minimal,
                                             std::forward<StepCallback>(on_step));
  }

  template <typename StepCallback>
  MorseSequence build_flooding_minmax_with_step_callback(StepCallback&& on_step) const {
    return build_flooding_with_step_callback(FloodingScheme::MinMax,
                                             std::forward<StepCallback>(on_step));
  }

  template <typename StepCallback>
  MorseSequence build_flooding_maxmin_with_step_callback(StepCallback&& on_step) const {
    return build_flooding_with_step_callback(FloodingScheme::MaxMin,
                                             std::forward<StepCallback>(on_step));
  }

  template <typename StepCallback>
  MorseSequence build_flooding_with_step_callback(FloodingScheme scheme,
                                                  StepCallback&& on_step) const {
    const std::size_t n = complex_.size();
    MorseSequence sequence(n);
    auto&& callback = on_step;

    struct LevelEvent {
      MorseStepType type = MorseStepType::Critical;
      SimplexId sigma = kInvalidSimplex;
      SimplexId tau = kInvalidSimplex;
    };

    std::vector<std::uint8_t> inserted(n, 0);
    std::vector<std::uint8_t> active(n, 0);
    std::vector<std::uint32_t> active_boundary_count(n, 0);
    std::vector<std::uint32_t> active_coboundary_count(n, 0);
    std::vector<SimplexId> active_boundary_xor(n, 0);
    std::vector<SimplexId> active_coboundary_xor(n, 0);

    auto emit_event = [&](const LevelEvent& event, LevelId level) {
      if (event.type == MorseStepType::Critical) {
        if (inserted[event.sigma]) {
          throw std::logic_error("Tried to insert a flooding critical twice.");
        }
        if constexpr (kValidateSequenceBuilder) {
          for (SimplexId face : complex_.boundary(event.sigma)) {
            if (!inserted[face]) {
              throw std::logic_error("Flooding critical has a missing boundary face.");
            }
          }
        }
        sequence.add_critical(event.sigma, level);
        callback(sequence, sequence.steps().back());
        inserted[event.sigma] = 1;
        return;
      }

      if (inserted[event.sigma] || inserted[event.tau]) {
        throw std::logic_error("Tried to insert a flooding regular pair twice.");
      }
      if constexpr (kValidateSequenceBuilder) {
        for (SimplexId face : complex_.boundary(event.tau)) {
          if (face != event.sigma && !inserted[face]) {
            throw std::logic_error("Flooding regular pair has a missing boundary face.");
          }
        }
      }
      sequence.add_regular_pair(event.sigma, event.tau, level);
      callback(sequence, sequence.steps().back());
      inserted[event.sigma] = 1;
      inserted[event.tau] = 1;
    };

    for (LevelId level = 0; level < complex_.num_levels(); ++level) {
      const auto& bucket = complex_.simplices_of_level(level);
      std::priority_queue<
          SimplexId,
          std::vector<SimplexId>,
          MaxDimensionMaxSimplexPriority>
          coreduction_candidates(MaxDimensionMaxSimplexPriority{this});
      std::priority_queue<
          SimplexId,
          std::vector<SimplexId>,
          MinDimensionMinSimplexPriority>
          reduction_candidates(MinDimensionMinSimplexPriority{this});
      std::priority_queue<
          SimplexId,
          std::vector<SimplexId>,
          MaxDimensionMaxSimplexPriority>
          coperforation_candidates(MaxDimensionMaxSimplexPriority{this});
      std::priority_queue<
          SimplexId,
          std::vector<SimplexId>,
          MinDimensionMinSimplexPriority>
          perforation_candidates(MinDimensionMinSimplexPriority{this});
      std::vector<LevelEvent> left_events;
      std::vector<LevelEvent> right_events;
      left_events.reserve(bucket.size());
      right_events.reserve(bucket.size());
      std::size_t remaining = bucket.size();

      for (SimplexId simplex : bucket) {
        active[simplex] = 1;
        active_boundary_count[simplex] = 0;
        active_coboundary_count[simplex] = 0;
        active_boundary_xor[simplex] = 0;
        active_coboundary_xor[simplex] = 0;
      }

      for (SimplexId simplex : bucket) {
        for (SimplexId face : complex_.boundary(simplex)) {
          if (simplex_level(face) == level && active[face]) {
            ++active_boundary_count[simplex];
            active_boundary_xor[simplex] ^= face;
          }
        }
        for (SimplexId coface : complex_.coboundary(simplex)) {
          if (simplex_level(coface) == level && active[coface]) {
            ++active_coboundary_count[simplex];
            active_coboundary_xor[simplex] ^= coface;
          }
        }
      }

      auto is_coreduction = [&](SimplexId tau) {
        return active[tau] && active_boundary_count[tau] == 1 &&
               active_boundary_xor[tau] < n && active[active_boundary_xor[tau]];
      };
      auto is_reduction = [&](SimplexId sigma) {
        return active[sigma] && active_coboundary_count[sigma] == 1 &&
               active_coboundary_xor[sigma] < n && active[active_coboundary_xor[sigma]];
      };
      auto is_coperforation = [&](SimplexId simplex) {
        return active[simplex] && active_boundary_count[simplex] == 0;
      };
      auto is_perforation = [&](SimplexId simplex) {
        return active[simplex] && active_coboundary_count[simplex] == 0;
      };

      auto enqueue_if_candidate = [&](SimplexId simplex) {
        if (!active[simplex]) {
          return;
        }
        if (is_coreduction(simplex)) {
          coreduction_candidates.push(simplex);
        }
        if (is_reduction(simplex)) {
          reduction_candidates.push(simplex);
        }
        if (is_coperforation(simplex)) {
          coperforation_candidates.push(simplex);
        }
        if (is_perforation(simplex)) {
          perforation_candidates.push(simplex);
        }
      };

      for (SimplexId simplex : bucket) {
        enqueue_if_candidate(simplex);
      }

      auto remove_simplex = [&](SimplexId simplex) {
        if (!active[simplex]) {
          throw std::logic_error("Tried to remove an inactive flooding simplex.");
        }
        active[simplex] = 0;
        --remaining;

        for (SimplexId coface : complex_.coboundary(simplex)) {
          if (simplex_level(coface) != level || !active[coface]) {
            continue;
          }
          if (active_boundary_count[coface] == 0) {
            throw std::logic_error("Active boundary count underflow.");
          }
          --active_boundary_count[coface];
          active_boundary_xor[coface] ^= simplex;
          enqueue_if_candidate(coface);
        }

        for (SimplexId face : complex_.boundary(simplex)) {
          if (simplex_level(face) != level || !active[face]) {
            continue;
          }
          if (active_coboundary_count[face] == 0) {
            throw std::logic_error("Active coboundary count underflow.");
          }
          --active_coboundary_count[face];
          active_coboundary_xor[face] ^= simplex;
          enqueue_if_candidate(face);
        }
      };

      auto take_coreduction = [&]() {
        while (!coreduction_candidates.empty()) {
          const SimplexId tau = coreduction_candidates.top();
          coreduction_candidates.pop();
          if (!is_coreduction(tau)) {
            continue;
          }
          const SimplexId sigma = active_boundary_xor[tau];
          left_events.push_back(LevelEvent{MorseStepType::RegularPair, sigma, tau});
          remove_simplex(sigma);
          remove_simplex(tau);
          return true;
        }
        return false;
      };

      auto take_reduction = [&]() {
        while (!reduction_candidates.empty()) {
          const SimplexId sigma = reduction_candidates.top();
          reduction_candidates.pop();
          if (!is_reduction(sigma)) {
            continue;
          }
          const SimplexId tau = active_coboundary_xor[sigma];
          right_events.push_back(LevelEvent{MorseStepType::RegularPair, sigma, tau});
          remove_simplex(sigma);
          remove_simplex(tau);
          return true;
        }
        return false;
      };

      auto take_coperforation = [&]() {
        while (!coperforation_candidates.empty()) {
          const SimplexId simplex = coperforation_candidates.top();
          coperforation_candidates.pop();
          if (!is_coperforation(simplex)) {
            continue;
          }
          left_events.push_back(LevelEvent{MorseStepType::Critical, simplex, kInvalidSimplex});
          remove_simplex(simplex);
          return true;
        }
        return false;
      };

      auto take_perforation = [&]() {
        while (!perforation_candidates.empty()) {
          const SimplexId simplex = perforation_candidates.top();
          perforation_candidates.pop();
          if (!is_perforation(simplex)) {
            continue;
          }
          right_events.push_back(LevelEvent{MorseStepType::Critical, simplex, kInvalidSimplex});
          remove_simplex(simplex);
          return true;
        }
        return false;
      };

      while (remaining > 0) {
        bool removed = false;
        if (scheme == FloodingScheme::Maximal) {
          removed = take_coreduction() || take_coperforation();
        } else if (scheme == FloodingScheme::Minimal) {
          removed = take_reduction() || take_perforation();
        } else if (scheme == FloodingScheme::MinMax) {
          removed = take_reduction() || take_coreduction() ||
                    take_perforation() || take_coperforation();
        } else {
          removed = take_coreduction() || take_reduction() ||
                    take_coperforation() || take_perforation();
        }
        if (!removed) {
          throw std::logic_error("No valid flooding operation found.");
        }
      }

      for (const LevelEvent& event : left_events) {
        emit_event(event, level);
      }
      for (std::size_t index = right_events.size(); index > 0; --index) {
        emit_event(right_events[index - 1], level);
      }
    }

    return sequence;
  }

 private:
  const ComplexView& complex_;
  std::vector<std::size_t> simplex_order_rank_;
  std::vector<LevelId> simplex_levels_;
  std::vector<std::uint16_t> simplex_dimensions_;
  MorseSequenceBuildMetrics* sequence_metrics_ = nullptr;
};

template <class ComplexView>
FSequenceBuilder(const ComplexView&) -> FSequenceBuilder<ComplexView>;

}  // namespace morseframes
