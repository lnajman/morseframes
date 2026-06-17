/*    Candidate documentation for the Gudhi Library - https://gudhi.inria.fr/
 *    The Gudhi Library is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Laurent Najman
 *
 *    Copyright (C) 2026
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_MORSE_PERSISTENCE_INTRO_MORSE_PERSISTENCE_H_
#define DOC_MORSE_PERSISTENCE_INTRO_MORSE_PERSISTENCE_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
namespace morse_persistence {

/** \defgroup morse_persistence Morse persistence
 * @{
 * \author    Laurent Najman
 *
 * \section morse_persistence_intro Introduction
 *
 * This module computes ordinary persistence from a filtered
 * `Gudhi::Simplex_tree` using same-level Morse reductions. It is designed for
 * filtrations with plateaus: equal filtration values are treated as genuine
 * levels where face/coface pairings may be selected directly, without refining
 * the input into a lower-star or simplex-wise filtration.
 *
 * The main entry point is @ref compute_morse_persistence. It builds a
 * same-level Morse sequence, constructs the associated reference map, and
 * reduces the resulting annotations to obtain a persistence diagram. The
 * result keeps a view of the input `Simplex_tree`, so local simplex ids in the
 * Morse sequence and the diagram can be mapped back to
 * `Simplex_tree::Simplex_handle` values while the input complex remains alive
 * and structurally unchanged.
 *
 * \section morse_persistence_scope Scope
 *
 * The first public version is intentionally narrow:
 *
 * \li input complexes are filtered `Gudhi::Simplex_tree` instances;
 * \li the computed invariant is ordinary persistence;
 * \li coefficients are in \f$\mathbb Z/2\mathbb Z\f$;
 * \li the exposed strategies are @ref Morse_sequence_strategy::SAME_LEVEL_REDUCTION,
 * @ref Morse_sequence_strategy::F_MAX, @ref Morse_sequence_strategy::F_MIN, and
 * @ref Morse_sequence_strategy::PLATEAU_GREEDY.
 *
 * The implementation keeps the Morse sequence and reference-map construction
 * fused, so the first reduction pass does not need to materialize a separate
 * Morse complex before persistence is computed.
 *
 * \section morse_persistence_example Example
 *
 * A minimal example is available in
 * `src/Morse_persistence/example/example_morse_persistence_from_simplex_tree.cpp`.
 *
 * @}
 */

}  // namespace morse_persistence
}  // namespace Gudhi

#endif  // DOC_MORSE_PERSISTENCE_INTRO_MORSE_PERSISTENCE_H_
