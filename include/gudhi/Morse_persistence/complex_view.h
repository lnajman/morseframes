#pragma once

/** \file gudhi/Morse_persistence/complex_view.h
 *  \brief Read-only views used to run Morse persistence on existing complexes.
 */

#include "morseframes/simplex_tree_builder.hpp"

namespace Gudhi {
namespace morse_persistence {

/** \brief Read-only adapter from a GUDHI-style `Simplex_tree` to dense ids.
 *
 *  \ingroup morse_persistence
 *
 *  The view stores local contiguous ids, boundary/coboundary lists, filtration
 *  levels, and a map back to `Simplex_tree::Simplex_handle`. The input
 *  Simplex-tree must outlive the view.
 */
template <class SimplexTree>
using Simplex_tree_view = morseframes::SimplexTreeComplexView<SimplexTree>;

/** \brief Owning filtered complex used for validation and imports.
 *
 *  \ingroup morse_persistence
 */
using Filtered_simplicial_complex = morseframes::FilteredSimplicialComplex;

/** \brief Build a read-only Morse view over an existing Simplex-tree. */
template <class SimplexTree>
Simplex_tree_view<SimplexTree> make_simplex_tree_view(const SimplexTree& simplex_tree) {
  return Simplex_tree_view<SimplexTree>(simplex_tree);
}

/** \brief Copy a Simplex-tree into an owning filtered complex.
 *
 *  The direct view path is the preferred GUDHI-facing entry point. This import
 *  helper remains useful for validation and file-based experiments.
 */
template <class SimplexTree>
Filtered_simplicial_complex filtered_complex_from_simplex_tree(
    const SimplexTree& simplex_tree) {
  return morseframes::filtered_complex_from_simplex_tree(simplex_tree);
}

}  // namespace morse_persistence
}  // namespace Gudhi
