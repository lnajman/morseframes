#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "morseframes/annotation.hpp"
#include "morseframes/complex_view.hpp"
#include "morseframes/filtered_complex.hpp"
#include "morseframes/morse_sequence.hpp"
#include "morseframes/reference_persistence.hpp"

namespace morseframes {

inline void require_debug_condition(bool condition, const char* message) {
  if (!condition) {
    throw std::logic_error(message);
  }
}

template <class ComplexView>
inline bool boundary_contains(const ComplexView& complex,
                              SimplexId simplex,
                              SimplexId face) {
  static_assert(is_complex_view_v<ComplexView>,
                "boundary_contains requires a Morse complex-view type.");
  const auto& boundary = complex.boundary(simplex);
  return std::find(boundary.begin(), boundary.end(), face) != boundary.end();
}

template <class ComplexView>
inline void validate_morse_sequence(const ComplexView& complex,
                                    const MorseSequence& sequence) {
  static_assert(is_complex_view_v<ComplexView>,
                "validate_morse_sequence requires a Morse complex-view type.");
  std::vector<bool> inserted(complex.size(), false);

  for (const MorseStep& step : sequence.steps()) {
    if (step.type == MorseStepType::Critical) {
      for (SimplexId face : complex.boundary(step.sigma)) {
        require_debug_condition(inserted[face], "Critical filling has a missing boundary face.");
      }
      require_debug_condition(!inserted[step.sigma], "Critical simplex inserted twice.");
      inserted[step.sigma] = true;
      continue;
    }

    require_debug_condition(!inserted[step.sigma], "Regular lower simplex inserted twice.");
    require_debug_condition(!inserted[step.tau], "Regular upper simplex inserted twice.");
    require_debug_condition(complex.dimension(step.tau) == complex.dimension(step.sigma) + 1,
                            "Regular pair dimensions are invalid.");
    require_debug_condition(boundary_contains(complex, step.tau, step.sigma),
                            "Regular lower simplex is not a face of the upper simplex.");
    require_debug_condition(complex.level(step.sigma) == complex.level(step.tau),
                            "Regular pair crosses filtration levels.");

    for (SimplexId face : complex.boundary(step.tau)) {
      if (face != step.sigma) {
        require_debug_condition(inserted[face], "Regular pair has a missing boundary face.");
      }
    }

    inserted[step.sigma] = true;
    inserted[step.tau] = true;
  }

  for (bool was_inserted : inserted) {
    require_debug_condition(was_inserted, "Morse sequence did not insert every simplex.");
  }
}

template <class ComplexView>
inline void validate_reference_invariants(const ComplexView& complex,
                                          const MorseSequence& sequence,
                                          const std::vector<Annotation>& references) {
  static_assert(is_complex_view_v<ComplexView>,
                "validate_reference_invariants requires a Morse complex-view type.");
  require_debug_condition(references.size() == complex.size(), "Reference table has wrong size.");

  for (const MorseStep& step : sequence.steps()) {
    if (step.type != MorseStepType::RegularPair) {
      continue;
    }

    require_debug_condition(references[step.tau].empty(),
                            "Reference of regular upper simplex is not zero.");

    Annotation boundary_reference;
    for (SimplexId face : complex.boundary(step.tau)) {
      xor_annotations_in_place(boundary_reference, references[face]);
    }
    require_debug_condition(boundary_reference.empty(),
                            "Reference of regular upper boundary is not zero.");
  }
}

template <class ComplexView>
inline void validate_coreference_invariants(const ComplexView& complex,
                                            const MorseSequence& sequence,
                                            const std::vector<Annotation>& coreferences) {
  static_assert(is_complex_view_v<ComplexView>,
                "validate_coreference_invariants requires a Morse complex-view type.");
  require_debug_condition(coreferences.size() == complex.size(), "Coreference table has wrong size.");

  for (const MorseStep& step : sequence.steps()) {
    if (step.type != MorseStepType::RegularPair) {
      continue;
    }

    require_debug_condition(coreferences[step.sigma].empty(),
                            "Coreference of regular lower simplex is not zero.");

    Annotation coboundary_coreference;
    for (SimplexId coface : complex.coboundary(step.sigma)) {
      xor_annotations_in_place(coboundary_coreference, coreferences[coface]);
    }
    require_debug_condition(coboundary_coreference.empty(),
                            "Coreference of regular lower coboundary is not zero.");
  }
}

inline void validate_persistence_diagram(const PersistenceDiagram& diagram) {
  for (const auto& pair : diagram.finite_pairs) {
    require_debug_condition(pair.birth != kInvalidSimplex, "Persistence pair has invalid birth.");
    require_debug_condition(pair.death != kInvalidSimplex, "Persistence pair has invalid death.");
    require_debug_condition(pair.birth_value <= pair.death_value,
                            "Persistence pair has birth after death.");
  }

  for (const auto& interval : diagram.essential) {
    require_debug_condition(interval.birth != kInvalidSimplex,
                            "Essential interval has invalid birth.");
  }
}

}  // namespace morseframes
