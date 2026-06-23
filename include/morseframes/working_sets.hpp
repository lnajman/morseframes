#pragma once

#include <vector>

#include "morseframes/complex_view.hpp"
#include "morseframes/filtered_complex.hpp"
#include "morseframes/morse_sequence.hpp"

namespace morseframes {

template <class ComplexView>
inline std::vector<SimplexId> reference_working_set(const ComplexView& complex,
                                                    const MorseSequence& sequence) {
  static_assert(is_complex_view_v<ComplexView>,
                "reference_working_set requires a Morse complex-view type.");
  std::vector<bool> present(complex.size(), false);

  for (SimplexId critical : sequence.critical_simplices()) {
    present[critical] = true;
    for (SimplexId face : complex.boundary(critical)) {
      present[face] = true;
    }
  }

  std::vector<SimplexId> result;
  for (SimplexId simplex = 0; simplex < present.size(); ++simplex) {
    if (present[simplex]) {
      result.push_back(simplex);
    }
  }
  return result;
}

template <class ComplexView>
inline std::vector<SimplexId> coreference_working_set(const ComplexView& complex,
                                                      const MorseSequence& sequence) {
  static_assert(is_complex_view_v<ComplexView>,
                "coreference_working_set requires a Morse complex-view type.");
  std::vector<bool> present(complex.size(), false);

  for (SimplexId critical : sequence.critical_simplices()) {
    present[critical] = true;
    for (SimplexId coface : complex.coboundary(critical)) {
      present[coface] = true;
    }
  }

  std::vector<SimplexId> result;
  for (SimplexId simplex = 0; simplex < present.size(); ++simplex) {
    if (present[simplex]) {
      result.push_back(simplex);
    }
  }
  return result;
}

}  // namespace morseframes
