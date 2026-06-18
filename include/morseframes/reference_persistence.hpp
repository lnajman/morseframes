#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "morseframes/annotation.hpp"
#include "morseframes/field_annotation_store.hpp"
#include "morseframes/filtered_complex.hpp"
#include "morseframes/inverse_annotation_store.hpp"
#include "morseframes/morse_sequence.hpp"
#include "morseframes/working_sets.hpp"

namespace morseframes {

struct PersistencePair {
  SimplexId birth = kInvalidSimplex;
  SimplexId death = kInvalidSimplex;
  std::uint16_t dimension = 0;
  double birth_value = 0.0;
  double death_value = 0.0;
};

struct EssentialInterval {
  SimplexId birth = kInvalidSimplex;
  std::uint16_t dimension = 0;
  double birth_value = 0.0;
};

struct PersistenceDiagram {
  std::vector<PersistencePair> finite_pairs;
  std::vector<EssentialInterval> essential;
};

struct MorseReferenceReductionMetrics {
  std::size_t working_set_size = 0;
  std::size_t critical_count = 0;
  std::uint64_t reducer_setup_nanoseconds = 0;
  std::uint64_t reducer_compute_nanoseconds = 0;
  std::size_t boundary_plan_face_scans = 0;
  std::size_t boundary_annotation_candidate_criticals = 0;
  std::size_t boundary_annotation_zero_skipped_criticals = 0;
  std::size_t boundary_annotation_zero_skipped_faces = 0;
  std::size_t boundary_annotation_xors = 0;
  std::size_t boundary_annotation_total_input_size = 0;
  std::size_t boundary_annotation_total_output_size = 0;
  std::size_t boundary_annotation_max_size = 0;
  std::size_t boundary_annotation_max_output_size = 0;
  std::size_t pivot_eliminations = 0;
  std::size_t finite_pairs = 0;
  std::size_t essential_intervals = 0;
  InverseAnnotationStoreMetrics inverse_store;
};

struct MorseReferenceReductionResult {
  PersistenceDiagram diagram;
  MorseReferenceReductionMetrics metrics;
};

struct MorseReferenceFrame {
  MorseSequence sequence;
  std::vector<Annotation> references;
};

struct MorseFieldReferenceFrame {
  MorseSequence sequence;
  std::vector<FieldAnnotation> references;
  std::uint32_t modulus = 2;
};

enum class ReferenceFrameReleasePolicy {
  Eager,
  Deferred,
};

struct MorseReferenceFrameMetrics {
  std::uint64_t remaining_cofaces_nanoseconds = 0;
  std::uint64_t sequence_total_nanoseconds = 0;
  std::uint64_t sequence_core_nanoseconds = 0;
  std::uint64_t reference_update_nanoseconds = 0;
  std::uint64_t reduction_plan_nanoseconds = 0;
  std::uint64_t release_cleanup_nanoseconds = 0;
  std::uint64_t working_set_pack_nanoseconds = 0;
  std::uint64_t local_index_nanoseconds = 0;
  std::uint64_t sequence_init_nanoseconds = 0;
  std::uint64_t sequence_candidate_seed_nanoseconds = 0;
  std::uint64_t sequence_candidate_loop_nanoseconds = 0;
  std::uint64_t sequence_emit_nanoseconds = 0;
  std::uint64_t sequence_callback_nanoseconds = 0;
  std::uint64_t sequence_replay_nanoseconds = 0;
  std::size_t sequence_candidate_pushes = 0;
  std::size_t sequence_candidate_pops = 0;
  std::size_t sequence_stale_candidate_skips = 0;
  std::size_t sequence_level_mismatch_skips = 0;
  std::size_t sequence_regular_pairs = 0;
  std::size_t sequence_criticals = 0;
  std::size_t final_live_nonempty_annotations = 0;
  std::size_t final_live_total_annotation_size = 0;
  std::size_t peak_live_nonempty_annotations = 0;
  std::size_t peak_live_total_annotation_size = 0;
  std::size_t released_annotations = 0;
  std::size_t released_total_annotation_size = 0;
};

struct ReferenceReductionCandidate {
  CriticalId critical_id = 0;
  SimplexId simplex = kInvalidSimplex;
  std::size_t first_boundary_face = 0;
  std::size_t boundary_face_count = 0;
};

struct ReferenceReductionFace {
  SimplexId face = kInvalidSimplex;
  std::size_t boundary_index = 0;
  std::size_t local_index = std::numeric_limits<std::size_t>::max();
};

struct ReferenceReductionPlan {
  std::vector<SimplexId> working_set;
  std::vector<ReferenceReductionCandidate> boundary_candidates;
  std::vector<ReferenceReductionFace> boundary_annotation_faces;
  std::vector<CriticalId> zero_boundary_critical_ids;
  std::size_t boundary_face_scans = 0;
  std::size_t zero_boundary_skipped_faces = 0;
};

struct MorseReferenceReductionInput {
  MorseSequence sequence;
  ReferenceReductionPlan reduction_plan;
  std::vector<Annotation> annotations;
  MorseReferenceFrameMetrics frame_metrics;
};

struct MorseReferenceProfile {
  std::size_t num_simplices = 0;
  std::size_t num_levels = 0;
  std::size_t num_critical_simplices = 0;
  std::size_t num_regular_pairs = 0;
  MorseReferenceFrameMetrics frame_metrics;
  MorseReferenceReductionMetrics reduction_metrics;
  std::size_t estimated_reducer_work = 0;
};

inline CriticalId checked_reduction_plan_critical_id(std::size_t value) {
  if (value > std::numeric_limits<CriticalId>::max()) {
    throw std::overflow_error("Too many critical simplices for CriticalId.");
  }
  return static_cast<CriticalId>(value);
}

inline void assign_reduction_plan_local_indices(std::size_t universe_size,
                                                ReferenceReductionPlan& plan) {
  constexpr std::size_t invalid = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> simplex_to_local(universe_size, invalid);
  for (std::size_t local = 0; local < plan.working_set.size(); ++local) {
    const SimplexId simplex = plan.working_set[local];
    if (simplex >= universe_size) {
      throw std::out_of_range("Working set contains an invalid simplex id.");
    }
    if (simplex_to_local[simplex] != invalid) {
      throw std::invalid_argument("Working set contains a duplicate simplex id.");
    }
    simplex_to_local[simplex] = local;
  }
  for (ReferenceReductionFace& face : plan.boundary_annotation_faces) {
    if (face.face >= universe_size || simplex_to_local[face.face] == invalid) {
      throw std::logic_error("Boundary annotation face is outside the working set.");
    }
    face.local_index = simplex_to_local[face.face];
  }
}

template <class ComplexView>
ReferenceReductionPlan build_reference_reduction_plan(
    const ComplexView& complex,
    const MorseSequence& sequence,
    const std::vector<Annotation>& references) {
  static_assert(is_complex_view_v<ComplexView>,
                "build_reference_reduction_plan requires a Morse complex-view type.");
  std::vector<std::uint8_t> present(complex.size(), 0);
  ReferenceReductionPlan plan;

  const auto& critical_simplices = sequence.critical_simplices();
  plan.boundary_candidates.reserve(critical_simplices.size());
  plan.zero_boundary_critical_ids.reserve(critical_simplices.size());

  for (std::size_t index = 0; index < critical_simplices.size(); ++index) {
    const CriticalId critical_id = checked_reduction_plan_critical_id(index);
    const SimplexId critical = critical_simplices[index];
    present[critical] = 1;

    bool may_have_nonzero_boundary = false;
    const auto& boundary = complex.boundary(critical);
    plan.boundary_face_scans += boundary.size();
    const std::size_t first_boundary_face = plan.boundary_annotation_faces.size();
    for (std::size_t boundary_index = 0; boundary_index < boundary.size(); ++boundary_index) {
      const SimplexId face = boundary[boundary_index];
      present[face] = 1;
      if (!references[face].empty()) {
        may_have_nonzero_boundary = true;
        plan.boundary_annotation_faces.push_back(
            ReferenceReductionFace{face, boundary_index});
      }
    }

    // Empty annotations stay empty: pivot updates only touch annotations that
    // already contain the pivot label.
    if (may_have_nonzero_boundary) {
      plan.boundary_candidates.push_back(ReferenceReductionCandidate{
          critical_id,
          critical,
          first_boundary_face,
          plan.boundary_annotation_faces.size() - first_boundary_face});
    } else {
      plan.zero_boundary_critical_ids.push_back(critical_id);
      plan.zero_boundary_skipped_faces += boundary.size();
    }
  }

  for (SimplexId simplex = 0; simplex < present.size(); ++simplex) {
    if (present[simplex]) {
      plan.working_set.push_back(simplex);
    }
  }
  assign_reduction_plan_local_indices(complex.size(), plan);
  return plan;
}

template <class ComplexView>
ReferenceReductionPlan build_reference_reduction_plan(
    const ComplexView& complex,
    const MorseSequence& sequence,
    const std::vector<FieldAnnotation>& references) {
  static_assert(is_complex_view_v<ComplexView>,
                "build_reference_reduction_plan requires a Morse complex-view type.");
  std::vector<std::uint8_t> present(complex.size(), 0);
  ReferenceReductionPlan plan;

  const auto& critical_simplices = sequence.critical_simplices();
  plan.boundary_candidates.reserve(critical_simplices.size());
  plan.zero_boundary_critical_ids.reserve(critical_simplices.size());

  for (std::size_t index = 0; index < critical_simplices.size(); ++index) {
    const CriticalId critical_id = checked_reduction_plan_critical_id(index);
    const SimplexId critical = critical_simplices[index];
    present[critical] = 1;

    bool may_have_nonzero_boundary = false;
    const auto& boundary = complex.boundary(critical);
    plan.boundary_face_scans += boundary.size();
    const std::size_t first_boundary_face = plan.boundary_annotation_faces.size();
    for (std::size_t boundary_index = 0; boundary_index < boundary.size(); ++boundary_index) {
      const SimplexId face = boundary[boundary_index];
      present[face] = 1;
      if (!references[face].empty()) {
        may_have_nonzero_boundary = true;
        plan.boundary_annotation_faces.push_back(
            ReferenceReductionFace{face, boundary_index});
      }
    }

    if (may_have_nonzero_boundary) {
      plan.boundary_candidates.push_back(ReferenceReductionCandidate{
          critical_id,
          critical,
          first_boundary_face,
          plan.boundary_annotation_faces.size() - first_boundary_face});
    } else {
      plan.zero_boundary_critical_ids.push_back(critical_id);
      plan.zero_boundary_skipped_faces += boundary.size();
    }
  }

  for (SimplexId simplex = 0; simplex < present.size(); ++simplex) {
    if (present[simplex]) {
      plan.working_set.push_back(simplex);
    }
  }
  assign_reduction_plan_local_indices(complex.size(), plan);
  return plan;
}

template <class ComplexView>
MorseReferenceProfile profile_morse_reference_reduction_input(
    const ComplexView& complex,
    const MorseReferenceReductionInput& input) {
  static_assert(is_complex_view_v<ComplexView>,
                "profile_morse_reference_reduction_input requires a Morse complex-view type.");
  MorseReferenceProfile profile;
  profile.num_simplices = complex.size();
  profile.num_levels = complex.num_levels();
  profile.num_critical_simplices = input.sequence.critical_simplices().size();
  profile.frame_metrics = input.frame_metrics;

  for (const MorseStep& step : input.sequence.steps()) {
    if (step.type == MorseStepType::RegularPair) {
      ++profile.num_regular_pairs;
    }
  }

  auto& metrics = profile.reduction_metrics;
  metrics.working_set_size = input.reduction_plan.working_set.size();
  metrics.critical_count = profile.num_critical_simplices;
  metrics.boundary_plan_face_scans = input.reduction_plan.boundary_face_scans;
  metrics.boundary_annotation_candidate_criticals =
      input.reduction_plan.boundary_candidates.size();
  metrics.boundary_annotation_zero_skipped_criticals =
      input.reduction_plan.zero_boundary_critical_ids.size();
  metrics.boundary_annotation_zero_skipped_faces =
      input.reduction_plan.zero_boundary_skipped_faces;

  for (const Annotation& annotation : input.annotations) {
    if (!annotation.empty()) {
      ++metrics.inverse_store.initial_nonempty_annotations;
      metrics.inverse_store.initial_total_annotation_size += annotation.size();
      metrics.inverse_store.initial_inverse_list_entries += annotation.size();
      if (annotation.size() > metrics.inverse_store.initial_max_annotation_size) {
        metrics.inverse_store.initial_max_annotation_size = annotation.size();
      }
    }
  }

  for (const auto& candidate : input.reduction_plan.boundary_candidates) {
    Annotation boundary_annotation;
    const std::size_t last_boundary_face =
        candidate.first_boundary_face + candidate.boundary_face_count;
    for (std::size_t face_index = candidate.first_boundary_face;
         face_index < last_boundary_face;
         ++face_index) {
      const std::size_t index =
          input.reduction_plan.boundary_annotation_faces[face_index].local_index;
      const Annotation& face_annotation = input.annotations.at(index);
      if (face_annotation.empty()) {
        continue;
      }
      ++metrics.boundary_annotation_xors;
      metrics.boundary_annotation_total_input_size += face_annotation.size();
      xor_annotations_in_place(boundary_annotation, face_annotation);
    }
    if (boundary_annotation.size() > metrics.boundary_annotation_max_size) {
      metrics.boundary_annotation_max_size = boundary_annotation.size();
    }
    metrics.boundary_annotation_total_output_size += boundary_annotation.size();
    if (boundary_annotation.size() > metrics.boundary_annotation_max_output_size) {
      metrics.boundary_annotation_max_output_size = boundary_annotation.size();
    }
  }

  profile.estimated_reducer_work =
      metrics.boundary_plan_face_scans +
      metrics.inverse_store.initial_total_annotation_size +
      metrics.boundary_annotation_xors +
      metrics.boundary_annotation_total_input_size +
      metrics.boundary_annotation_total_output_size;
  return profile;
}

template <class ComplexView = FilteredSimplicialComplex>
class MorseReferenceComputer {
  static_assert(is_complex_view_v<ComplexView>,
                "MorseReferenceComputer requires a Morse complex-view type.");

 public:
  MorseReferenceComputer(const ComplexView& complex, const MorseSequence& sequence)
      : complex_(complex), sequence_(sequence) {}

  std::vector<Annotation> compute_full_references() const {
    std::vector<Annotation> references(complex_.size());

    for (const MorseStep& step : sequence_.steps()) {
      if (step.type == MorseStepType::Critical) {
        const auto critical_id = checked_critical_id(sequence_.critical_index(step.sigma));
        references[step.sigma] = Annotation{critical_id};
        continue;
      }

      references[step.tau].clear();
      Annotation lower_reference;
      for (SimplexId face : complex_.boundary(step.tau)) {
        if (face != step.sigma) {
          xor_annotations_in_place(lower_reference, references[face]);
        }
      }
      references[step.sigma] = std::move(lower_reference);
    }

    return references;
  }

 private:
  static CriticalId checked_critical_id(std::int32_t value) {
    if (value < 0) {
      throw std::logic_error("Expected a critical simplex.");
    }
    return static_cast<CriticalId>(value);
  }

  const ComplexView& complex_;
  const MorseSequence& sequence_;
};

template <class ComplexView = FilteredSimplicialComplex>
class MorseReferenceFrameBuilder {
  static_assert(is_complex_view_v<ComplexView>,
                "MorseReferenceFrameBuilder requires a Morse complex-view type.");

 public:
  explicit MorseReferenceFrameBuilder(const ComplexView& complex,
                                      bool collect_frame_timing = false,
                                      ReferenceFrameReleasePolicy release_policy =
                                          ReferenceFrameReleasePolicy::Eager)
      : complex_(complex),
        collect_frame_timing_(collect_frame_timing),
        release_policy_(release_policy) {}

  MorseReferenceFrame build_saturated() const {
    std::vector<Annotation> references(complex_.size());
    Annotation reference_update_scratch;
    auto sequence = FSequenceBuilder(complex_).build_saturated_with_step_callback(
        [&](const MorseSequence& sequence, const MorseStep& step) {
          update_reference_for_step(sequence, step, references, reference_update_scratch);
        });
    return MorseReferenceFrame{std::move(sequence), std::move(references)};
  }

  MorseReferenceFrame build_plateau_greedy() const {
    std::vector<Annotation> references(complex_.size());
    Annotation reference_update_scratch;
    auto sequence = FSequenceBuilder(complex_).build_plateau_greedy_with_step_callback(
        [&](const MorseSequence& sequence, const MorseStep& step) {
          update_reference_for_step(sequence, step, references, reference_update_scratch);
        });
    return MorseReferenceFrame{std::move(sequence), std::move(references)};
  }

  MorseReferenceFrame build_coreduction() const {
    return build_same_level_reduction();
  }

  MorseReferenceFrame build_same_level_reduction() const {
    std::vector<Annotation> references(complex_.size());
    Annotation reference_update_scratch;
    auto sequence = FSequenceBuilder(complex_).build_same_level_reduction_with_step_callback(
        [&](const MorseSequence& sequence, const MorseStep& step) {
          update_reference_for_step(sequence, step, references, reference_update_scratch);
        });
    return MorseReferenceFrame{std::move(sequence), std::move(references)};
  }

  MorseReferenceFrame build_f_max() const {
    std::vector<Annotation> references(complex_.size());
    Annotation reference_update_scratch;
    auto sequence = FSequenceBuilder(complex_).build_f_max_with_step_callback(
        [&](const MorseSequence& sequence, const MorseStep& step) {
          update_reference_for_step(sequence, step, references, reference_update_scratch);
        });
    return MorseReferenceFrame{std::move(sequence), std::move(references)};
  }

  MorseReferenceFrame build_f_min() const {
    std::vector<Annotation> references(complex_.size());
    Annotation reference_update_scratch;
    auto sequence = FSequenceBuilder(complex_).build_f_min_with_step_callback(
        [&](const MorseSequence& sequence, const MorseStep& step) {
          update_reference_for_step(sequence, step, references, reference_update_scratch);
        });
    return MorseReferenceFrame{std::move(sequence), std::move(references)};
  }

  MorseReferenceFrame build_flooding_max() const {
    std::vector<Annotation> references(complex_.size());
    Annotation reference_update_scratch;
    auto sequence = FSequenceBuilder(complex_).build_flooding_max_with_step_callback(
        [&](const MorseSequence& sequence, const MorseStep& step) {
          update_reference_for_step(sequence, step, references, reference_update_scratch);
        });
    return MorseReferenceFrame{std::move(sequence), std::move(references)};
  }

  MorseReferenceFrame build_flooding_min() const {
    std::vector<Annotation> references(complex_.size());
    Annotation reference_update_scratch;
    auto sequence = FSequenceBuilder(complex_).build_flooding_min_with_step_callback(
        [&](const MorseSequence& sequence, const MorseStep& step) {
          update_reference_for_step(sequence, step, references, reference_update_scratch);
        });
    return MorseReferenceFrame{std::move(sequence), std::move(references)};
  }

  MorseReferenceFrame build_flooding_minmax() const {
    std::vector<Annotation> references(complex_.size());
    Annotation reference_update_scratch;
    auto sequence = FSequenceBuilder(complex_).build_flooding_minmax_with_step_callback(
        [&](const MorseSequence& sequence, const MorseStep& step) {
          update_reference_for_step(sequence, step, references, reference_update_scratch);
        });
    return MorseReferenceFrame{std::move(sequence), std::move(references)};
  }

  MorseReferenceFrame build_flooding_maxmin() const {
    std::vector<Annotation> references(complex_.size());
    Annotation reference_update_scratch;
    auto sequence = FSequenceBuilder(complex_).build_flooding_maxmin_with_step_callback(
        [&](const MorseSequence& sequence, const MorseStep& step) {
          update_reference_for_step(sequence, step, references, reference_update_scratch);
        });
    return MorseReferenceFrame{std::move(sequence), std::move(references)};
  }

  MorseReferenceReductionInput build_saturated_reduction_input() const {
    return build_reduction_input_with([&](auto&& step_callback, auto* sequence_metrics) {
      return FSequenceBuilder(complex_, sequence_metrics).build_saturated_with_step_callback(
          std::forward<decltype(step_callback)>(step_callback));
    });
  }

  MorseReferenceReductionInput build_plateau_greedy_reduction_input() const {
    return build_reduction_input_with([&](auto&& step_callback, auto* sequence_metrics) {
      return FSequenceBuilder(complex_, sequence_metrics).build_plateau_greedy_with_step_callback(
          std::forward<decltype(step_callback)>(step_callback));
    });
  }

  MorseReferenceReductionInput build_coreduction_reduction_input() const {
    return build_same_level_reduction_reduction_input();
  }

  MorseReferenceReductionInput build_same_level_reduction_reduction_input() const {
    return build_reduction_input_with([&](auto&& step_callback, auto* sequence_metrics) {
      return FSequenceBuilder(complex_, sequence_metrics)
          .build_same_level_reduction_with_step_callback(
          std::forward<decltype(step_callback)>(step_callback));
    });
  }

  MorseReferenceReductionInput build_f_max_reduction_input() const {
    return build_reduction_input_with([&](auto&& step_callback, auto* sequence_metrics) {
      return FSequenceBuilder(complex_, sequence_metrics).build_f_max_with_step_callback(
          std::forward<decltype(step_callback)>(step_callback));
    });
  }

  MorseReferenceReductionInput build_f_min_reduction_input() const {
    return build_reduction_input_with([&](auto&& step_callback, auto* sequence_metrics) {
      return FSequenceBuilder(complex_, sequence_metrics).build_f_min_with_step_callback(
          std::forward<decltype(step_callback)>(step_callback));
    });
  }

  MorseReferenceReductionInput build_flooding_max_reduction_input() const {
    return build_reduction_input_with([&](auto&& step_callback, auto* sequence_metrics) {
      return FSequenceBuilder(complex_, sequence_metrics).build_flooding_max_with_step_callback(
          std::forward<decltype(step_callback)>(step_callback));
    });
  }

  MorseReferenceReductionInput build_flooding_min_reduction_input() const {
    return build_reduction_input_with([&](auto&& step_callback, auto* sequence_metrics) {
      return FSequenceBuilder(complex_, sequence_metrics).build_flooding_min_with_step_callback(
          std::forward<decltype(step_callback)>(step_callback));
    });
  }

  MorseReferenceReductionInput build_flooding_minmax_reduction_input() const {
    return build_reduction_input_with([&](auto&& step_callback, auto* sequence_metrics) {
      return FSequenceBuilder(complex_, sequence_metrics).build_flooding_minmax_with_step_callback(
          std::forward<decltype(step_callback)>(step_callback));
    });
  }

  MorseReferenceReductionInput build_flooding_maxmin_reduction_input() const {
    return build_reduction_input_with([&](auto&& step_callback, auto* sequence_metrics) {
      return FSequenceBuilder(complex_, sequence_metrics).build_flooding_maxmin_with_step_callback(
          std::forward<decltype(step_callback)>(step_callback));
    });
  }

 private:
  using Clock = std::chrono::steady_clock;

  static std::uint64_t elapsed_nanoseconds(Clock::time_point start,
                                           Clock::time_point stop) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count());
  }

  static void accumulate_nanoseconds(std::uint64_t& total,
                                     Clock::time_point start,
                                     Clock::time_point stop) {
    total += elapsed_nanoseconds(start, stop);
  }

  static CriticalId checked_critical_id(std::int32_t value) {
    if (value < 0) {
      throw std::logic_error("Expected a critical simplex.");
    }
    return static_cast<CriticalId>(value);
  }

  static CriticalId current_callback_critical_id(const MorseSequence& sequence) {
    const auto& critical_simplices = sequence.critical_simplices();
    if (critical_simplices.empty()) {
      throw std::logic_error("Critical callback has no critical simplex.");
    }
    return checked_reduction_plan_critical_id(critical_simplices.size() - 1);
  }

  void reserve_reduction_plan_storage(ReferenceReductionPlan& plan) const {
    const std::size_t simplex_count = complex_.size();
    plan.working_set.reserve(simplex_count);
    plan.boundary_candidates.reserve(simplex_count);
    plan.boundary_annotation_faces.reserve(simplex_count);
    plan.zero_boundary_critical_ids.reserve(simplex_count);
  }

  void update_reference_for_step(const MorseSequence& sequence,
                                 const MorseStep& step,
                                 std::vector<Annotation>& references,
                                 Annotation& scratch) const {
    if (step.type == MorseStepType::Critical) {
      references[step.sigma] =
          Annotation{current_callback_critical_id(sequence)};
      return;
    }

    references[step.tau].clear();
    Annotation lower_reference;
    for (SimplexId face : complex_.boundary(step.tau)) {
      if (face != step.sigma && !references[face].empty()) {
        xor_annotations_in_place(lower_reference, references[face], scratch);
      }
    }
    references[step.sigma] = std::move(lower_reference);
  }

  void update_reference_for_step(const MorseSequence& sequence,
                                 const MorseStep& step,
                                 std::vector<Annotation>& references) const {
    Annotation scratch;
    update_reference_for_step(sequence, step, references, scratch);
  }

  template <typename ReplaceReference>
  void update_reference_for_step(const MorseSequence& sequence,
                                 const MorseStep& step,
                                 const std::vector<Annotation>& references,
                                 Annotation& scratch,
                                 ReplaceReference&& replace_reference) const {
    if (step.type == MorseStepType::Critical) {
      replace_reference(
          step.sigma, Annotation{current_callback_critical_id(sequence)});
      return;
    }

    if (!references[step.tau].empty()) {
      replace_reference(step.tau, Annotation{});
    }
    Annotation lower_reference;
    for (SimplexId face : complex_.boundary(step.tau)) {
      if (face != step.sigma && !references[face].empty()) {
        xor_annotations_in_place(lower_reference, references[face], scratch);
      }
    }
    replace_reference(step.sigma, std::move(lower_reference));
  }

  template <typename ReplaceReference>
  void update_reference_for_step(const MorseSequence& sequence,
                                 const MorseStep& step,
                                 const std::vector<Annotation>& references,
                                 ReplaceReference&& replace_reference) const {
    Annotation scratch;
    update_reference_for_step(sequence, step, references, scratch, replace_reference);
  }

  template <typename MarkPresent>
  void update_reduction_plan_for_step(const MorseSequence& sequence,
                                      const MorseStep& step,
                                      const std::vector<Annotation>& references,
                                      MarkPresent&& mark_present,
                                      ReferenceReductionPlan& plan) const {
    if (step.type != MorseStepType::Critical) {
      return;
    }

    const CriticalId critical_id = current_callback_critical_id(sequence);
    mark_present(step.sigma);

    bool may_have_nonzero_boundary = false;
    const auto& boundary = complex_.boundary(step.sigma);
    plan.boundary_face_scans += boundary.size();
    const std::size_t first_boundary_face = plan.boundary_annotation_faces.size();
    for (std::size_t boundary_index = 0; boundary_index < boundary.size(); ++boundary_index) {
      const SimplexId face = boundary[boundary_index];
      mark_present(face);
      if (!references[face].empty()) {
        may_have_nonzero_boundary = true;
        plan.boundary_annotation_faces.push_back(
            ReferenceReductionFace{face, boundary_index});
      }
    }

    if (may_have_nonzero_boundary) {
      plan.boundary_candidates.push_back(ReferenceReductionCandidate{
          critical_id,
          step.sigma,
          first_boundary_face,
          plan.boundary_annotation_faces.size() - first_boundary_face});
    } else {
      plan.zero_boundary_critical_ids.push_back(critical_id);
      plan.zero_boundary_skipped_faces += boundary.size();
    }
  }

  void update_reduction_plan_for_step(const MorseSequence& sequence,
                                      const MorseStep& step,
                                      const std::vector<Annotation>& references,
                                      std::vector<std::uint8_t>& present,
                                      ReferenceReductionPlan& plan) const {
    auto mark_present = [&](SimplexId simplex) {
      present[simplex] = 1;
    };
    update_reduction_plan_for_step(sequence, step, references, mark_present, plan);
  }

  template <typename BuildSequence>
  MorseReferenceReductionInput build_reduction_input_with(BuildSequence&& build_sequence) const {
    std::vector<Annotation> references(complex_.size());
    std::vector<std::uint8_t> present(complex_.size(), 0);
    std::vector<std::uint32_t> remaining_cofaces(complex_.size(), 0);
    ReferenceReductionPlan plan;
    MorseReferenceFrameMetrics frame_metrics;
    Annotation reference_update_scratch;
    MorseSequenceBuildMetrics sequence_build_metrics;
    reserve_reduction_plan_storage(plan);
    const bool eager_release = release_policy_ == ReferenceFrameReleasePolicy::Eager;
    const bool track_live_reference_metrics = collect_frame_timing_ || eager_release;
    auto* sequence_metrics = collect_frame_timing_ ? &sequence_build_metrics : nullptr;

    const auto remaining_start =
        collect_frame_timing_ ? Clock::now() : Clock::time_point{};
    if (eager_release) {
      for (SimplexId simplex = 0; simplex < complex_.size(); ++simplex) {
        remaining_cofaces[simplex] =
            static_cast<std::uint32_t>(complex_.coboundary(simplex).size());
      }
    }
    if (collect_frame_timing_) {
      const auto remaining_done = Clock::now();
      frame_metrics.remaining_cofaces_nanoseconds =
          elapsed_nanoseconds(remaining_start, remaining_done);
    }

    auto replace_reference = [&](SimplexId simplex, Annotation annotation) {
      if (!track_live_reference_metrics) {
        references[simplex] = std::move(annotation);
        return;
      }

      if (!references[simplex].empty()) {
        --frame_metrics.final_live_nonempty_annotations;
        frame_metrics.final_live_total_annotation_size -= references[simplex].size();
      }

      const std::size_t new_size = annotation.size();
      references[simplex] = std::move(annotation);
      if (new_size == 0) {
        return;
      }

      ++frame_metrics.final_live_nonempty_annotations;
      frame_metrics.final_live_total_annotation_size += new_size;
      if (frame_metrics.final_live_nonempty_annotations >
          frame_metrics.peak_live_nonempty_annotations) {
        frame_metrics.peak_live_nonempty_annotations =
            frame_metrics.final_live_nonempty_annotations;
      }
      if (frame_metrics.final_live_total_annotation_size >
          frame_metrics.peak_live_total_annotation_size) {
        frame_metrics.peak_live_total_annotation_size =
            frame_metrics.final_live_total_annotation_size;
      }
    };

    auto move_reference = [&](SimplexId simplex) {
      Annotation annotation = std::move(references[simplex]);
      if (track_live_reference_metrics && !annotation.empty()) {
        --frame_metrics.final_live_nonempty_annotations;
        frame_metrics.final_live_total_annotation_size -= annotation.size();
      }
      return annotation;
    };

    auto release_if_done = [&](SimplexId simplex) {
      if (remaining_cofaces[simplex] == 0 &&
          !present[simplex] &&
          !references[simplex].empty()) {
        Annotation released = move_reference(simplex);
        ++frame_metrics.released_annotations;
        frame_metrics.released_total_annotation_size += released.size();
      }
    };

    auto mark_coface_processed = [&](SimplexId coface) {
      if (!eager_release) {
        return;
      }
      for (SimplexId face : complex_.boundary(coface)) {
        if (remaining_cofaces[face] == 0) {
          throw std::logic_error("Remaining coface count underflow.");
        }
        --remaining_cofaces[face];
        release_if_done(face);
      }
      release_if_done(coface);
    };

    auto mark_present = [&](SimplexId simplex) {
      present[simplex] = 1;
    };

    auto sequence = [&]() {
      if (!collect_frame_timing_) {
        return build_sequence([&](const MorseSequence& sequence, const MorseStep& step) {
          update_reference_for_step(
              sequence, step, references, reference_update_scratch, replace_reference);
          update_reduction_plan_for_step(sequence,
                                         step,
                                         references,
                                         mark_present,
                                         plan);
          if (step.type == MorseStepType::Critical) {
            mark_coface_processed(step.sigma);
          } else {
            mark_coface_processed(step.sigma);
            mark_coface_processed(step.tau);
          }
        }, sequence_metrics);
      }

      const auto sequence_start = Clock::now();
      auto timed_sequence = build_sequence(
          [&](const MorseSequence& sequence, const MorseStep& step) {
            const auto reference_start = Clock::now();
            update_reference_for_step(
                sequence, step, references, reference_update_scratch, replace_reference);
            const auto reference_done = Clock::now();
            accumulate_nanoseconds(frame_metrics.reference_update_nanoseconds,
                                   reference_start,
                                   reference_done);

            const auto plan_start = Clock::now();
            update_reduction_plan_for_step(sequence,
                                           step,
                                           references,
                                           mark_present,
                                           plan);
            const auto plan_done = Clock::now();
            accumulate_nanoseconds(frame_metrics.reduction_plan_nanoseconds,
                                   plan_start,
                                   plan_done);

            const auto release_start = Clock::now();
            if (step.type == MorseStepType::Critical) {
              mark_coface_processed(step.sigma);
            } else {
              mark_coface_processed(step.sigma);
              mark_coface_processed(step.tau);
            }
            const auto release_done = Clock::now();
            accumulate_nanoseconds(frame_metrics.release_cleanup_nanoseconds,
                                   release_start,
                                   release_done);
          },
          sequence_metrics);
      const auto sequence_done = Clock::now();
      frame_metrics.sequence_total_nanoseconds =
          elapsed_nanoseconds(sequence_start, sequence_done);
      const std::uint64_t callback_nanoseconds =
          frame_metrics.reference_update_nanoseconds +
          frame_metrics.reduction_plan_nanoseconds +
          frame_metrics.release_cleanup_nanoseconds;
      frame_metrics.sequence_core_nanoseconds =
          frame_metrics.sequence_total_nanoseconds > callback_nanoseconds
              ? frame_metrics.sequence_total_nanoseconds - callback_nanoseconds
              : 0;
      return timed_sequence;
    }();

    if (collect_frame_timing_) {
      frame_metrics.sequence_init_nanoseconds =
          sequence_build_metrics.init_nanoseconds;
      frame_metrics.sequence_candidate_seed_nanoseconds =
          sequence_build_metrics.candidate_seed_nanoseconds;
      frame_metrics.sequence_candidate_loop_nanoseconds =
          sequence_build_metrics.candidate_loop_nanoseconds;
      frame_metrics.sequence_emit_nanoseconds =
          sequence_build_metrics.emit_nanoseconds;
      frame_metrics.sequence_callback_nanoseconds =
          sequence_build_metrics.callback_nanoseconds;
      frame_metrics.sequence_replay_nanoseconds =
          sequence_build_metrics.replay_nanoseconds;
      frame_metrics.sequence_candidate_pushes =
          sequence_build_metrics.candidate_pushes;
      frame_metrics.sequence_candidate_pops =
          sequence_build_metrics.candidate_pops;
      frame_metrics.sequence_stale_candidate_skips =
          sequence_build_metrics.stale_candidate_skips;
      frame_metrics.sequence_level_mismatch_skips =
          sequence_build_metrics.level_mismatch_skips;
      frame_metrics.sequence_regular_pairs =
          sequence_build_metrics.regular_pairs;
      frame_metrics.sequence_criticals =
          sequence_build_metrics.criticals;
    }

    const auto pack_start =
        collect_frame_timing_ ? Clock::now() : Clock::time_point{};
    std::vector<Annotation> annotations;
    for (SimplexId simplex = 0; simplex < present.size(); ++simplex) {
      if (present[simplex]) {
        plan.working_set.push_back(simplex);
      }
    }
    annotations.reserve(plan.working_set.size());
    for (SimplexId simplex : plan.working_set) {
      annotations.push_back(move_reference(simplex));
    }
    if (collect_frame_timing_) {
      const auto pack_done = Clock::now();
      frame_metrics.working_set_pack_nanoseconds =
          elapsed_nanoseconds(pack_start, pack_done);
    }

    const auto local_index_start =
        collect_frame_timing_ ? Clock::now() : Clock::time_point{};
    assign_reduction_plan_local_indices(complex_.size(), plan);
    if (collect_frame_timing_) {
      const auto local_index_done = Clock::now();
      frame_metrics.local_index_nanoseconds =
          elapsed_nanoseconds(local_index_start, local_index_done);
    }

    return MorseReferenceReductionInput{
        std::move(sequence), std::move(plan), std::move(annotations), frame_metrics};
  }

  const ComplexView& complex_;
  bool collect_frame_timing_ = false;
  ReferenceFrameReleasePolicy release_policy_ = ReferenceFrameReleasePolicy::Eager;
};

template <class ComplexView>
MorseReferenceComputer(const ComplexView&, const MorseSequence&)
    -> MorseReferenceComputer<ComplexView>;

template <class ComplexView>
MorseReferenceFrameBuilder(const ComplexView&,
                           bool = false,
                           ReferenceFrameReleasePolicy = ReferenceFrameReleasePolicy::Eager)
    -> MorseReferenceFrameBuilder<ComplexView>;

template <class ComplexView = FilteredSimplicialComplex>
class MorseFieldReferenceComputer {
  static_assert(is_complex_view_v<ComplexView>,
                "MorseFieldReferenceComputer requires a Morse complex-view type.");

 public:
  MorseFieldReferenceComputer(const ComplexView& complex,
                              const MorseSequence& sequence,
                              std::uint32_t modulus)
      : complex_(complex), sequence_(sequence), modulus_(modulus) {
    validate_prime_field_characteristic(modulus_);
  }

  std::vector<FieldAnnotation> compute_full_references() const {
    std::vector<FieldAnnotation> references(complex_.size());

    for (const MorseStep& step : sequence_.steps()) {
      update_reference_for_step(step, references);
    }

    return references;
  }

 private:
  static CriticalId checked_critical_id(std::int32_t value) {
    if (value < 0) {
      throw std::logic_error("Expected a critical simplex.");
    }
    return static_cast<CriticalId>(value);
  }

  void update_reference_for_step(const MorseStep& step,
                                 std::vector<FieldAnnotation>& references) const {
    if (step.type == MorseStepType::Critical) {
      references[step.sigma] = FieldAnnotation{
          FieldAnnotationEntry{checked_critical_id(sequence_.critical_index(step.sigma)), 1}};
      return;
    }

    references[step.tau].clear();
    FieldAnnotation lower_reference;
    std::uint32_t paired_face_coefficient = 0;
    bool found_paired_face = false;
    const auto& boundary = complex_.boundary(step.tau);
    for (std::size_t removed_index = 0; removed_index < boundary.size(); ++removed_index) {
      const SimplexId face = boundary[removed_index];
      const std::uint32_t coefficient = boundary_coefficient(removed_index, modulus_);
      if (face == step.sigma) {
        paired_face_coefficient = coefficient;
        found_paired_face = true;
        continue;
      }
      add_scaled_field_annotation_in_place(
          lower_reference, references[face], coefficient, modulus_);
    }

    if (!found_paired_face) {
      throw std::logic_error("Regular pair is not a face/coface pair.");
    }

    const std::uint32_t scale =
        (modulus_ - modp_inverse(paired_face_coefficient, modulus_)) % modulus_;
    scale_field_annotation_in_place(lower_reference, scale, modulus_);
    references[step.sigma] = std::move(lower_reference);
  }

  const ComplexView& complex_;
  const MorseSequence& sequence_;
  std::uint32_t modulus_ = 2;
};

template <class ComplexView>
MorseFieldReferenceComputer(const ComplexView&, const MorseSequence&, std::uint32_t)
    -> MorseFieldReferenceComputer<ComplexView>;

template <class ComplexView = FilteredSimplicialComplex>
class MorseReferencePersistenceReducer {
  static_assert(is_complex_view_v<ComplexView>,
                "MorseReferencePersistenceReducer requires a Morse complex-view type.");

 public:
  MorseReferencePersistenceReducer(const ComplexView& complex,
                                   const MorseSequence& sequence,
                                   const std::vector<Annotation>& references)
      : complex_(complex),
        sequence_(sequence),
        reduction_plan_(build_reference_reduction_plan(complex, sequence, references)),
        annotations_(references,
                     reduction_plan_.working_set,
                     sequence.critical_simplices().size()) {}

  MorseReferencePersistenceReducer(const ComplexView& complex,
                                   const MorseSequence& sequence,
                                   ReferenceReductionPlan reduction_plan,
                                   std::vector<Annotation> annotations)
      : complex_(complex),
        sequence_(sequence),
        reduction_plan_(std::move(reduction_plan)),
        annotations_(std::move(annotations),
                     reduction_plan_.working_set,
                     complex.size(),
                     sequence.critical_simplices().size()) {}

  PersistenceDiagram compute() {
    return compute_with_metrics().diagram;
  }

  MorseReferenceReductionResult compute_with_metrics() {
    const auto& critical_simplices = sequence_.critical_simplices();
    std::vector<std::uint8_t> active(critical_simplices.size(), 0);
    std::vector<std::uint8_t> killed(critical_simplices.size(), 0);
    MorseReferenceReductionResult result;
    result.diagram.finite_pairs.reserve(reduction_plan_.boundary_candidates.size());
    result.diagram.essential.reserve(critical_simplices.size());
    result.metrics.working_set_size = annotations_.size();
    result.metrics.critical_count = critical_simplices.size();
    result.metrics.boundary_plan_face_scans = reduction_plan_.boundary_face_scans;
    result.metrics.boundary_annotation_candidate_criticals =
        reduction_plan_.boundary_candidates.size();
    result.metrics.boundary_annotation_zero_skipped_criticals =
        reduction_plan_.zero_boundary_critical_ids.size();
    result.metrics.boundary_annotation_zero_skipped_faces =
        reduction_plan_.zero_boundary_skipped_faces;

    for (CriticalId critical_id : reduction_plan_.zero_boundary_critical_ids) {
      active[critical_id] = 1;
    }

    Annotation boundary_annotation;
    for (const auto& candidate : reduction_plan_.boundary_candidates) {
      const CriticalId sigma_critical_id = candidate.critical_id;
      const SimplexId sigma = candidate.simplex;
      boundary_annotation.clear();
      const std::size_t last_boundary_face =
          candidate.first_boundary_face + candidate.boundary_face_count;
      for (std::size_t face_index = candidate.first_boundary_face;
           face_index < last_boundary_face;
           ++face_index) {
        const auto& face =
            reduction_plan_.boundary_annotation_faces[face_index];
        const auto& face_annotation =
            annotations_.annotation_by_local_unchecked(face.local_index);
        if (face_annotation.empty()) {
          continue;
        }
        ++result.metrics.boundary_annotation_xors;
        result.metrics.boundary_annotation_total_input_size += face_annotation.size();
        xor_annotations_in_place(boundary_annotation, face_annotation);
      }
      if (boundary_annotation.size() > result.metrics.boundary_annotation_max_size) {
        result.metrics.boundary_annotation_max_size = boundary_annotation.size();
      }
      result.metrics.boundary_annotation_total_output_size += boundary_annotation.size();
      if (boundary_annotation.size() > result.metrics.boundary_annotation_max_output_size) {
        result.metrics.boundary_annotation_max_output_size = boundary_annotation.size();
      }

      if (boundary_annotation.empty()) {
        active[sigma_critical_id] = 1;
        continue;
      }

      const CriticalId pivot = boundary_annotation.back();
      const SimplexId birth = critical_simplices[pivot];
      const std::uint16_t dimension = complex_.dimension(birth);

      result.diagram.finite_pairs.push_back(PersistencePair{
          birth,
          sigma,
          dimension,
          complex_.filtration(birth),
          complex_.filtration(sigma),
      });

      killed[pivot] = 1;
      ++result.metrics.pivot_eliminations;

      remove_negative_label(sigma_critical_id);
      eliminate_pivot(pivot, boundary_annotation);
    }

    for (CriticalId id = 0; id < critical_simplices.size(); ++id) {
      if (!active[id] || killed[id]) {
        continue;
      }
      const SimplexId birth = critical_simplices[id];
      result.diagram.essential.push_back(EssentialInterval{
          birth,
          complex_.dimension(birth),
          complex_.filtration(birth),
      });
    }

    result.metrics.finite_pairs = result.diagram.finite_pairs.size();
    result.metrics.essential_intervals = result.diagram.essential.size();
    result.metrics.inverse_store = annotations_.metrics();
    return result;
  }

 private:
  static CriticalId checked_critical_id(std::int32_t value) {
    if (value < 0) {
      throw std::logic_error("Expected a critical simplex.");
    }
    return static_cast<CriticalId>(value);
  }

  void remove_negative_label(CriticalId label) {
    annotations_.remove_label_from_all(label);
  }

  void eliminate_pivot(CriticalId pivot, const Annotation& boundary_annotation) {
    annotations_.xor_into_all_containing_known_pivot(pivot, boundary_annotation);
  }

  const ComplexView& complex_;
  const MorseSequence& sequence_;
  ReferenceReductionPlan reduction_plan_;
  InverseAnnotationStore annotations_;
};

template <class ComplexView>
MorseReferencePersistenceReducer(const ComplexView&,
                                 const MorseSequence&,
                                 const std::vector<Annotation>&)
    -> MorseReferencePersistenceReducer<ComplexView>;

template <class ComplexView>
MorseReferencePersistenceReducer(const ComplexView&,
                                 const MorseSequence&,
                                 ReferenceReductionPlan,
                                 std::vector<Annotation>)
    -> MorseReferencePersistenceReducer<ComplexView>;

template <class ComplexView = FilteredSimplicialComplex>
class MorseFieldReferencePersistenceReducer {
  static_assert(is_complex_view_v<ComplexView>,
                "MorseFieldReferencePersistenceReducer requires a Morse complex-view type.");

 public:
  MorseFieldReferencePersistenceReducer(const ComplexView& complex,
                                        const MorseSequence& sequence,
                                        std::vector<FieldAnnotation> references,
                                        std::uint32_t modulus)
      : complex_(complex),
        sequence_(sequence),
        references_(std::move(references)),
        modulus_(modulus) {
    validate_prime_field_characteristic(modulus_);
    if (references_.size() != complex_.size()) {
      throw std::invalid_argument("Reference table size does not match the complex size.");
    }
  }

  PersistenceDiagram compute() {
    const auto& critical_simplices = sequence_.critical_simplices();
    std::vector<std::uint8_t> active(critical_simplices.size(), 0);
    std::vector<std::uint8_t> killed(critical_simplices.size(), 0);
    PersistenceDiagram diagram;

    for (CriticalId sigma_critical_id = 0;
         sigma_critical_id < critical_simplices.size();
         ++sigma_critical_id) {
      const SimplexId sigma = critical_simplices[sigma_critical_id];
      FieldAnnotation boundary_annotation;
      const auto& boundary = complex_.boundary(sigma);
      for (std::size_t removed_index = 0; removed_index < boundary.size(); ++removed_index) {
        add_scaled_field_annotation_in_place(
            boundary_annotation,
            references_[boundary[removed_index]],
            boundary_coefficient(removed_index, modulus_),
            modulus_);
      }

      if (boundary_annotation.empty()) {
        active[sigma_critical_id] = 1;
        continue;
      }

      const CriticalId pivot = boundary_annotation.back().label;
      const std::uint32_t pivot_coefficient = boundary_annotation.back().coefficient;
      const SimplexId birth = critical_simplices.at(pivot);

      diagram.finite_pairs.push_back(PersistencePair{
          birth,
          sigma,
          complex_.dimension(birth),
          complex_.filtration(birth),
          complex_.filtration(sigma),
      });

      killed[pivot] = 1;

      for (auto& annotation : references_) {
        remove_label_from_field_annotation(annotation, sigma_critical_id);
      }

      const std::uint32_t inverse_pivot = modp_inverse(pivot_coefficient, modulus_);
      for (auto& annotation : references_) {
        const std::uint32_t coefficient =
            field_annotation_coefficient(annotation, pivot);
        if (coefficient == 0) {
          continue;
        }
        const std::uint32_t scale =
            (modulus_ - modp_multiply(coefficient, inverse_pivot, modulus_)) % modulus_;
        add_scaled_field_annotation_in_place(
            annotation, boundary_annotation, scale, modulus_);
      }
    }

    for (CriticalId id = 0; id < critical_simplices.size(); ++id) {
      if (!active[id] || killed[id]) {
        continue;
      }
      const SimplexId birth = critical_simplices[id];
      diagram.essential.push_back(EssentialInterval{
          birth,
          complex_.dimension(birth),
          complex_.filtration(birth),
      });
    }

    return diagram;
  }

 private:
  const ComplexView& complex_;
  const MorseSequence& sequence_;
  std::vector<FieldAnnotation> references_;
  std::uint32_t modulus_ = 2;
};

template <class ComplexView>
MorseFieldReferencePersistenceReducer(const ComplexView&,
                                      const MorseSequence&,
                                      std::vector<FieldAnnotation>,
                                      std::uint32_t)
    -> MorseFieldReferencePersistenceReducer<ComplexView>;

template <class ComplexView = FilteredSimplicialComplex>
class MorseCompactFieldReferencePersistenceReducer {
  static_assert(is_complex_view_v<ComplexView>,
                "MorseCompactFieldReferencePersistenceReducer requires a Morse "
                "complex-view type.");

 public:
  MorseCompactFieldReferencePersistenceReducer(const ComplexView& complex,
                                               const MorseSequence& sequence,
                                               ReferenceReductionPlan reduction_plan,
                                               std::vector<FieldAnnotation> annotations,
                                               std::uint32_t modulus)
      : complex_(complex),
        sequence_(sequence),
        reduction_plan_(std::move(reduction_plan)),
        annotations_(std::move(annotations),
                     reduction_plan_.working_set,
                     complex.size(),
                     sequence.critical_simplices().size(),
                     modulus),
        modulus_(modulus) {
    validate_prime_field_characteristic(modulus_);
  }

  PersistenceDiagram compute() {
    return compute_with_metrics().diagram;
  }

  MorseReferenceReductionResult compute_with_metrics() {
    const auto& critical_simplices = sequence_.critical_simplices();
    std::vector<std::uint8_t> active(critical_simplices.size(), 0);
    std::vector<std::uint8_t> killed(critical_simplices.size(), 0);
    MorseReferenceReductionResult result;
    result.metrics.working_set_size = annotations_.size();
    result.metrics.critical_count = critical_simplices.size();
    result.metrics.boundary_plan_face_scans = reduction_plan_.boundary_face_scans;
    result.metrics.boundary_annotation_candidate_criticals =
        reduction_plan_.boundary_candidates.size();
    result.metrics.boundary_annotation_zero_skipped_criticals =
        reduction_plan_.zero_boundary_critical_ids.size();
    result.metrics.boundary_annotation_zero_skipped_faces =
        reduction_plan_.zero_boundary_skipped_faces;

    for (CriticalId critical_id : reduction_plan_.zero_boundary_critical_ids) {
      active[critical_id] = 1;
    }

    for (const auto& candidate : reduction_plan_.boundary_candidates) {
      const CriticalId sigma_critical_id = candidate.critical_id;
      const SimplexId sigma = candidate.simplex;
      FieldAnnotation boundary_annotation;
      const std::size_t last_boundary_face =
          candidate.first_boundary_face + candidate.boundary_face_count;
      for (std::size_t face_index = candidate.first_boundary_face;
           face_index < last_boundary_face;
           ++face_index) {
        const auto& face =
            reduction_plan_.boundary_annotation_faces[face_index];
        const auto& face_annotation =
            annotations_.annotation_by_local_unchecked(face.local_index);
        if (face_annotation.empty()) {
          continue;
        }
        ++result.metrics.boundary_annotation_xors;
        result.metrics.boundary_annotation_total_input_size += face_annotation.size();
        add_scaled_field_annotation_in_place(boundary_annotation,
                                             face_annotation,
                                             boundary_coefficient(face.boundary_index, modulus_),
                                             modulus_);
      }
      if (boundary_annotation.size() > result.metrics.boundary_annotation_max_size) {
        result.metrics.boundary_annotation_max_size = boundary_annotation.size();
      }
      result.metrics.boundary_annotation_total_output_size += boundary_annotation.size();
      if (boundary_annotation.size() > result.metrics.boundary_annotation_max_output_size) {
        result.metrics.boundary_annotation_max_output_size = boundary_annotation.size();
      }

      if (boundary_annotation.empty()) {
        active[sigma_critical_id] = 1;
        continue;
      }

      const CriticalId pivot = boundary_annotation.back().label;
      const std::uint32_t pivot_coefficient = boundary_annotation.back().coefficient;
      const SimplexId birth = critical_simplices.at(pivot);
      const std::uint16_t dimension = complex_.dimension(birth);

      result.diagram.finite_pairs.push_back(PersistencePair{
          birth,
          sigma,
          dimension,
          complex_.filtration(birth),
          complex_.filtration(sigma),
      });

      killed[pivot] = 1;
      ++result.metrics.pivot_eliminations;

      annotations_.remove_label_from_all(sigma_critical_id);
      annotations_.eliminate_pivot(pivot, boundary_annotation, pivot_coefficient);
    }

    for (CriticalId id = 0; id < critical_simplices.size(); ++id) {
      if (!active[id] || killed[id]) {
        continue;
      }
      const SimplexId birth = critical_simplices[id];
      result.diagram.essential.push_back(EssentialInterval{
          birth,
          complex_.dimension(birth),
          complex_.filtration(birth),
      });
    }

    result.metrics.finite_pairs = result.diagram.finite_pairs.size();
    result.metrics.essential_intervals = result.diagram.essential.size();
    result.metrics.inverse_store = annotations_.metrics();
    return result;
  }

 private:
  const ComplexView& complex_;
  const MorseSequence& sequence_;
  ReferenceReductionPlan reduction_plan_;
  FieldAnnotationStore annotations_;
  std::uint32_t modulus_ = 2;
};

template <class ComplexView>
MorseCompactFieldReferencePersistenceReducer(const ComplexView&,
                                             const MorseSequence&,
                                             ReferenceReductionPlan,
                                             std::vector<FieldAnnotation>,
                                             std::uint32_t)
    -> MorseCompactFieldReferencePersistenceReducer<ComplexView>;

template <class ComplexView>
inline PersistenceDiagram compute_morse_reference_prime_field_persistence(
    const ComplexView& complex,
    const MorseSequence& sequence,
    std::uint32_t modulus) {
  static_assert(is_complex_view_v<ComplexView>,
                "compute_morse_reference_prime_field_persistence requires a Morse "
                "complex-view type.");
  auto references =
      MorseFieldReferenceComputer(complex, sequence, modulus).compute_full_references();
  auto reduction_plan = build_reference_reduction_plan(complex, sequence, references);
  std::vector<FieldAnnotation> annotations;
  annotations.reserve(reduction_plan.working_set.size());
  for (SimplexId simplex : reduction_plan.working_set) {
    annotations.push_back(std::move(references.at(simplex)));
  }
  MorseCompactFieldReferencePersistenceReducer reducer(
      complex, sequence, std::move(reduction_plan), std::move(annotations), modulus);
  return reducer.compute();
}

template <class ComplexView>
inline PersistenceDiagram compute_morse_reference_persistence(
    const ComplexView& complex, const MorseSequence& sequence) {
  static_assert(is_complex_view_v<ComplexView>,
                "compute_morse_reference_persistence requires a Morse complex-view type.");
  MorseReferenceComputer reference_computer(complex, sequence);
  auto references = reference_computer.compute_full_references();
  MorseReferencePersistenceReducer reducer(complex, sequence, references);
  return reducer.compute();
}

template <class ComplexView>
inline PersistenceDiagram compute_fused_morse_reference_persistence(
    const ComplexView& complex) {
  static_assert(is_complex_view_v<ComplexView>,
                "compute_fused_morse_reference_persistence requires a Morse complex-view type.");
  auto input = MorseReferenceFrameBuilder(complex).build_saturated_reduction_input();
  MorseReferencePersistenceReducer reducer(complex,
                                           input.sequence,
                                           std::move(input.reduction_plan),
                                           std::move(input.annotations));
  return reducer.compute();
}

inline std::vector<PersistencePair> off_diagonal_pairs(const PersistenceDiagram& diagram) {
  std::vector<PersistencePair> result;
  for (const auto& pair : diagram.finite_pairs) {
    if (pair.birth_value < pair.death_value) {
      result.push_back(pair);
    }
  }
  return result;
}

}  // namespace morseframes
