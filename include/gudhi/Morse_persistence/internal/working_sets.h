#pragma once

#include <vector>

#include "gudhi/Morse_persistence/internal/filtered_complex.h"
#include "gudhi/Morse_persistence/internal/morse_sequence.h"

namespace Gudhi { namespace morse_persistence { namespace internal {

inline std::vector<SimplexId> reference_working_set(const FilteredSimplicialComplex& complex,
                                                    const MorseSequence& sequence) {
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

inline std::vector<SimplexId> coreference_working_set(const FilteredSimplicialComplex& complex,
                                                      const MorseSequence& sequence) {
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

}  // namespace internal
}  // namespace morse_persistence
}  // namespace Gudhi
