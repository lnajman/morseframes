#pragma once

#include <stdexcept>
#include <string>
#include <utility>

#include "gudhi/Morse_persistence/internal/reference_persistence.h"

namespace Gudhi { namespace morse_persistence { namespace internal {

enum class MorseSequenceStrategy {
  SameLevelReduction,
  FMax,
  FMin,
  PlateauGreedy,
};

inline const char* morse_sequence_strategy_name(MorseSequenceStrategy strategy) {
  switch (strategy) {
    case MorseSequenceStrategy::SameLevelReduction:
      return "same-level-reduction";
    case MorseSequenceStrategy::FMax:
      return "f-max";
    case MorseSequenceStrategy::FMin:
      return "f-min";
    case MorseSequenceStrategy::PlateauGreedy:
      return "plateau-greedy";
  }
  throw std::invalid_argument("Unknown Morse sequence strategy.");
}

inline MorseSequenceStrategy morse_sequence_strategy_from_name(const std::string& strategy) {
  if (strategy == "same-level-reduction") {
    return MorseSequenceStrategy::SameLevelReduction;
  }
  if (strategy == "f-max") {
    return MorseSequenceStrategy::FMax;
  }
  if (strategy == "f-min") {
    return MorseSequenceStrategy::FMin;
  }
  if (strategy == "plateau-greedy") {
    return MorseSequenceStrategy::PlateauGreedy;
  }
  throw std::invalid_argument("Unknown Morse sequence strategy: " + strategy);
}

template <class ComplexView>
MorseReferenceFrame build_morse_reference_frame(
    const ComplexView& complex,
    MorseSequenceStrategy strategy = MorseSequenceStrategy::FMax) {
  static_assert(is_complex_view_v<ComplexView>,
                "build_morse_reference_frame requires a Morse complex-view type.");

  MorseReferenceFrameBuilder builder(complex);
  switch (strategy) {
    case MorseSequenceStrategy::SameLevelReduction:
      return builder.build_same_level_reduction();
    case MorseSequenceStrategy::FMax:
      return builder.build_f_max();
    case MorseSequenceStrategy::FMin:
      return builder.build_f_min();
    case MorseSequenceStrategy::PlateauGreedy:
      return builder.build_plateau_greedy();
  }
  throw std::invalid_argument("Unknown Morse sequence strategy.");
}

template <class ComplexView>
MorseReferenceReductionInput build_morse_reference_reduction_input(
    const ComplexView& complex,
    MorseSequenceStrategy strategy = MorseSequenceStrategy::FMax,
    bool collect_frame_timing = false,
    ReferenceFrameReleasePolicy release_policy = ReferenceFrameReleasePolicy::Eager) {
  static_assert(is_complex_view_v<ComplexView>,
                "build_morse_reference_reduction_input requires a Morse complex-view type.");

  MorseReferenceFrameBuilder builder(complex, collect_frame_timing, release_policy);
  switch (strategy) {
    case MorseSequenceStrategy::SameLevelReduction:
      return builder.build_same_level_reduction_reduction_input();
    case MorseSequenceStrategy::FMax:
      return builder.build_f_max_reduction_input();
    case MorseSequenceStrategy::FMin:
      return builder.build_f_min_reduction_input();
    case MorseSequenceStrategy::PlateauGreedy:
      return builder.build_plateau_greedy_reduction_input();
  }
  throw std::invalid_argument("Unknown Morse sequence strategy.");
}

template <class ComplexView>
MorseReferenceReductionResult compute_morse_reference_persistence_with_metrics(
    const ComplexView& complex,
    MorseSequenceStrategy strategy = MorseSequenceStrategy::FMax) {
  static_assert(is_complex_view_v<ComplexView>,
                "compute_morse_reference_persistence_with_metrics requires a Morse "
                "complex-view type.");

  auto input = build_morse_reference_reduction_input(complex, strategy);
  MorseReferencePersistenceReducer reducer(complex,
                                           input.sequence,
                                           std::move(input.reduction_plan),
                                           std::move(input.annotations));
  return reducer.compute_with_metrics();
}

template <class ComplexView>
PersistenceDiagram compute_morse_reference_persistence(
    const ComplexView& complex,
    MorseSequenceStrategy strategy) {
  return compute_morse_reference_persistence_with_metrics(complex, strategy).diagram;
}

}  // namespace internal
}  // namespace morse_persistence
}  // namespace Gudhi
