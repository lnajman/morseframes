#pragma once

/** \file gudhi/Morse_persistence/diagram.h
 *  \brief Persistence diagram types returned by Morse persistence.
 */

#include <vector>

#include "gudhi/Morse_persistence/internal/reference_persistence.h"

namespace Gudhi {
namespace morse_persistence {

/** \brief Contiguous local simplex identifier used by the Morse reducer.
 *
 *  \ingroup morse_persistence
 */
using Simplex_id = internal::SimplexId;

/** \brief Vertex identifier used in canonical simplex vertex tuples.
 *
 *  \ingroup morse_persistence
 */
using Vertex_id = internal::VertexId;

/** \brief Index of a distinct filtration value.
 *
 *  \ingroup morse_persistence
 */
using Level_id = internal::LevelId;

/** \brief Sentinel used when no simplex id is available. */
constexpr Simplex_id invalid_simplex = internal::kInvalidSimplex;

/** \brief Finite persistence interval with birth and death simplex ids. */
using Persistence_pair = internal::PersistencePair;

/** \brief Essential persistence interval with a birth simplex id. */
using Essential_interval = internal::EssentialInterval;

/** \brief Persistence diagram over \f$\mathbb Z_2\f$.
 *
 *  \ingroup morse_persistence
 *
 *  Finite intervals may include zero-length intervals on plateaus. Use
 *  `off_diagonal_pairs()` when comparing public persistent features.
 */
using Persistence_diagram = internal::PersistenceDiagram;

/** \brief Return finite intervals whose birth and death filtration values differ. */
inline std::vector<Persistence_pair> off_diagonal_pairs(
    const Persistence_diagram& diagram) {
  return internal::off_diagonal_pairs(diagram);
}

}  // namespace morse_persistence
}  // namespace Gudhi
