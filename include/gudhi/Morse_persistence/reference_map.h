#pragma once

/** \file gudhi/Morse_persistence/reference_map.h
 *  \brief Fused Morse sequence and reference-map construction.
 */

#include <type_traits>
#include <utility>
#include <vector>

#include "gudhi/Morse_persistence/complex_view.h"
#include "gudhi/Morse_persistence/diagram.h"
#include "gudhi/Morse_persistence/morse_sequence.h"
#include "gudhi/Morse_persistence/strategy.h"
#include "gudhi/Morse_persistence/internal/reference_persistence.h"

namespace Gudhi {
namespace morse_persistence {

/** \brief Sparse annotation over critical simplex ids.
 *
 *  \ingroup morse_persistence
 */
using Annotation = morseframes::Annotation;

/** \brief Morse sequence together with the full reference map.
 *
 *  \ingroup morse_persistence
 */
using Morse_reference_frame = morseframes::MorseReferenceFrame;

/** \brief Metrics collected while constructing a Morse reference frame.
 *
 *  \ingroup morse_persistence
 */
using Morse_reference_frame_metrics = morseframes::MorseReferenceFrameMetrics;

/** \brief Compact reducer input built while computing a Morse sequence.
 *
 *  \ingroup morse_persistence
 */
using Morse_reference_reduction_input = morseframes::MorseReferenceReductionInput;

/** \brief Frame result for a Simplex-tree input.
 *
 *  \ingroup morse_persistence
 *
 *  The object owns the temporary view so the local simplex ids in the frame can
 *  still be mapped back to Simplex-tree handles. The input Simplex-tree must
 *  outlive this result.
 */
template <class SimplexTree>
struct Simplex_tree_morse_frame {
  Simplex_tree_view<SimplexTree> view;
  Morse_reference_frame frame;

  /** \brief Return the computed Morse sequence. */
  const Morse_sequence& sequence() const { return frame.sequence; }

  /** \brief Return the full reference map indexed by local simplex id. */
  const std::vector<Annotation>& reference_map() const { return frame.references; }

  /** \brief Map a local simplex id back to a Simplex-tree handle. */
  auto simplex_tree_handle(Simplex_id simplex) const ->
      typename Simplex_tree_view<SimplexTree>::SimplexHandle {
    return view.handle(simplex);
  }
};

/** \brief Compute a Morse sequence and full reference map on a complex view.
 *
 *  \ingroup morse_persistence
 *
 *  \param[in] complex Read-only complex view satisfying the Morse complex-view API.
 *  \param[in] strategy Same-level Morse sequence construction strategy.
 *  \return The computed Morse sequence and full reference map.
 */
template <class ComplexView,
          typename std::enable_if<morseframes::is_complex_view_v<ComplexView>, int>::type = 0>
Morse_reference_frame compute_morse_sequence_and_reference_map(
    const ComplexView& complex,
    Morse_sequence_strategy strategy = Morse_sequence_strategy::F_MAX) {
  return morseframes::build_morse_reference_frame(complex, to_kernel_strategy(strategy));
}

/** \brief Compute a compact reduction input on a complex view.
 *
 *  \ingroup morse_persistence
 *
 *  This fused path avoids materializing reference annotations that the reducer
 *  will never inspect.
 *
 *  \param[in] complex Read-only complex view satisfying the Morse complex-view API.
 *  \param[in] strategy Same-level Morse sequence construction strategy.
 *  \return A compact reduction input containing the sequence, reduction plan,
 *  and live annotations required by the reducer.
 */
template <class ComplexView,
          typename std::enable_if<morseframes::is_complex_view_v<ComplexView>, int>::type = 0>
Morse_reference_reduction_input compute_morse_reference_reduction_input(
    const ComplexView& complex,
    Morse_sequence_strategy strategy = Morse_sequence_strategy::F_MAX) {
  return morseframes::build_morse_reference_reduction_input(complex, to_kernel_strategy(strategy));
}

/** \brief Compute a Morse sequence and reference map directly from a Simplex-tree.
 *
 *  \ingroup morse_persistence
 *
 *  \param[in] simplex_tree Filtered Simplex-tree. It must outlive the returned result.
 *  \param[in] strategy Same-level Morse sequence construction strategy.
 *  \return A result owning the temporary view and the computed frame.
 */
template <class SimplexTree,
          typename std::enable_if<!morseframes::is_complex_view_v<SimplexTree>, int>::type = 0>
Simplex_tree_morse_frame<SimplexTree> compute_morse_sequence_and_reference_map(
    const SimplexTree& simplex_tree,
    Morse_sequence_strategy strategy = Morse_sequence_strategy::F_MAX) {
  auto view = make_simplex_tree_view(simplex_tree);
  auto frame = compute_morse_sequence_and_reference_map(view, strategy);
  return Simplex_tree_morse_frame<SimplexTree>{std::move(view), std::move(frame)};
}

}  // namespace morse_persistence
}  // namespace Gudhi
