#pragma once

/** \file gudhi/Morse_persistence.h
 *  \brief Umbrella header for the GUDHI Morse persistence API.
 *
 *  This header exposes a thin public layer over the current Morse persistence
 *  kernels. The main entry point is
 *  `Gudhi::morse_persistence::compute_morse_persistence`.
 *
 *  \defgroup morse_persistence Morse persistence
 *  \brief Same-level Morse sequence, reference-map, and persistence routines.
 *
 *  The module computes Morse persistence directly on a filtered
 *  `Gudhi::Simplex_tree` without refining plateaus into lower-star filtrations.
 *  Equal filtration values are treated as genuine levels where same-level
 *  face/coface pairings may be selected.
 *
 *  The first public GUDHI-facing API computes ordinary persistence over
 *  \f$\mathbb Z_2\f$. Prime-field kernels are kept in the prototype layer for
 *  experiments and can be promoted later without changing the Simplex-tree view
 *  abstraction.
 */

#include "gudhi/Morse_persistence/complex_view.h"
#include "gudhi/Morse_persistence/diagram.h"
#include "gudhi/Morse_persistence/morse_sequence.h"
#include "gudhi/Morse_persistence/persistence_reducer.h"
#include "gudhi/Morse_persistence/reference_map.h"
#include "gudhi/Morse_persistence/strategy.h"
