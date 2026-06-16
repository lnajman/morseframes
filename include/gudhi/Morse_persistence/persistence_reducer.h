#pragma once

/** \file gudhi/Morse_persistence/persistence_reducer.h
 *  \brief Public persistence computation entry points.
 */

#include <type_traits>
#include <utility>
#include <vector>

#include "gudhi/Morse_persistence/diagram.h"
#include "gudhi/Morse_persistence/reference_map.h"
#include "morseframes/simplex_tree_morse.hpp"

namespace Gudhi {
namespace morse_persistence {

/** \brief Metrics collected by the Morse-reference reducer.
 *
 *  \ingroup morse_persistence
 */
using Morse_reference_reduction_metrics = morseframes::MorseReferenceReductionMetrics;

/** \brief Persistence result with reducer metrics.
 *
 *  \ingroup morse_persistence
 */
using Morse_reference_reduction_result = morseframes::MorseReferenceReductionResult;

/** \brief Full result for direct Simplex-tree Morse persistence.
 *
 *  \ingroup morse_persistence
 *
 *  The result owns the temporary view so that local simplex ids in the sequence
 *  and diagram can be mapped back to Simplex-tree handles. The input
 *  Simplex-tree must outlive this result.
 */
template <class SimplexTree>
struct Simplex_tree_morse_result {
  Simplex_tree_view<SimplexTree> view;
  Morse_sequence sequence;
  Persistence_diagram diagram;
  Morse_reference_frame_metrics frame_metrics;
  Morse_reference_reduction_metrics reduction_metrics;

  /** \brief All finite intervals, including zero-length intervals on plateaus. */
  const std::vector<Persistence_pair>& finite_intervals() const {
    return diagram.finite_pairs;
  }

  /** \brief Finite intervals with birth value strictly smaller than death value. */
  std::vector<Persistence_pair> off_diagonal_intervals() const {
    return ::Gudhi::morse_persistence::off_diagonal_pairs(diagram);
  }

  /** \brief Essential intervals. */
  const std::vector<Essential_interval>& essential_intervals() const {
    return diagram.essential;
  }

  /** \brief Map a local simplex id back to a Simplex-tree handle. */
  auto simplex_tree_handle(Simplex_id simplex) const ->
      typename Simplex_tree_view<SimplexTree>::SimplexHandle {
    return view.handle(simplex);
  }
};

/** \brief Compute Morse persistence on a complex view.
 *
 *  \ingroup morse_persistence
 *
 *  \param[in] complex Read-only complex view satisfying the Morse complex-view API.
 *  \param[in] strategy Same-level Morse sequence construction strategy.
 *  \return Persistence diagram over \f$\mathbb Z_2\f$.
 */
template <class ComplexView,
          typename std::enable_if<morseframes::is_complex_view_v<ComplexView>, int>::type = 0>
Persistence_diagram compute_morse_persistence(
    const ComplexView& complex,
    Morse_sequence_strategy strategy = Morse_sequence_strategy::F_MAX) {
  return morseframes::compute_morse_reference_persistence(complex, to_kernel_strategy(strategy));
}

/** \brief Compute Morse persistence with reducer metrics on a complex view.
 *
 *  \ingroup morse_persistence
 *
 *  \param[in] complex Read-only complex view satisfying the Morse complex-view API.
 *  \param[in] strategy Same-level Morse sequence construction strategy.
 *  \return Persistence diagram over \f$\mathbb Z_2\f$ and reducer metrics.
 */
template <class ComplexView,
          typename std::enable_if<morseframes::is_complex_view_v<ComplexView>, int>::type = 0>
Morse_reference_reduction_result compute_morse_persistence_with_metrics(
    const ComplexView& complex,
    Morse_sequence_strategy strategy = Morse_sequence_strategy::F_MAX) {
  return morseframes::compute_morse_reference_persistence_with_metrics(
      complex,
      to_kernel_strategy(strategy));
}

/** \brief Reduce persistence from a precomputed Morse sequence/reference map.
 *
 *  \ingroup morse_persistence
 *
 *  \param[in] complex Read-only complex view used to compute the frame.
 *  \param[in] frame Precomputed Morse sequence and full reference map.
 *  \return Persistence diagram over \f$\mathbb Z_2\f$.
 */
template <class ComplexView,
          typename std::enable_if<morseframes::is_complex_view_v<ComplexView>, int>::type = 0>
Persistence_diagram compute_morse_persistence(
    const ComplexView& complex,
    const Morse_reference_frame& frame) {
  morseframes::MorseReferencePersistenceReducer reducer(
      complex,
      frame.sequence,
      frame.references);
  return reducer.compute();
}

/** \brief Reduce persistence from a precomputed Simplex-tree Morse frame.
 *
 *  \ingroup morse_persistence
 */
template <class SimplexTree>
Persistence_diagram compute_morse_persistence(
    const Simplex_tree_morse_frame<SimplexTree>& frame) {
  return compute_morse_persistence(frame.view, frame.frame);
}

/** \brief Reduce persistence from a precomputed Simplex-tree Morse frame.
 *
 *  \ingroup morse_persistence
 *
 *  The `simplex_tree` parameter documents the intended source object and lets
 *  call sites read naturally. The frame itself owns the view used by the
 *  reducer.
 */
template <class SimplexTree>
Persistence_diagram compute_morse_persistence(
    const SimplexTree& simplex_tree,
    const Simplex_tree_morse_frame<SimplexTree>& frame) {
  (void)simplex_tree;
  return compute_morse_persistence(frame);
}

/** \brief Compute Morse persistence directly from a Simplex-tree.
 *
 *  \ingroup morse_persistence
 *
 *  \param[in] simplex_tree Filtered Simplex-tree. It must outlive the returned result.
 *  \param[in] strategy Same-level Morse sequence construction strategy.
 *  \return Result containing the temporary view, Morse sequence, \f$\mathbb Z_2\f$
 *          diagram, and metrics.
 */
template <class SimplexTree,
          typename std::enable_if<!morseframes::is_complex_view_v<SimplexTree>, int>::type = 0>
Simplex_tree_morse_result<SimplexTree> compute_morse_persistence(
    const SimplexTree& simplex_tree,
    Morse_sequence_strategy strategy = Morse_sequence_strategy::F_MAX) {
  auto internal = morseframes::compute_simplex_tree_morse_reference_persistence(
      simplex_tree,
      to_kernel_strategy(strategy));
  return Simplex_tree_morse_result<SimplexTree>{
      std::move(internal.view),
      std::move(internal.sequence),
      std::move(internal.diagram),
      internal.frame_metrics,
      internal.reduction_metrics,
  };
}

}  // namespace morse_persistence
}  // namespace Gudhi
