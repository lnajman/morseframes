#pragma once

#include <utility>

#include "morseframes/morse_reference_api.hpp"
#include "morseframes/simplex_tree_builder.hpp"

namespace morseframes {

template <class SimplexTree>
using GudhiSimplexTreeComplexView = SimplexTreeComplexView<SimplexTree>;

template <class SimplexTree>
SimplexTreeComplexView<SimplexTree> make_simplex_tree_complex_view(
    const SimplexTree& simplex_tree,
    SimplexTreeFiltrationOrder order_policy =
        SimplexTreeFiltrationOrder::CanonicalLexicographic) {
  return SimplexTreeComplexView<SimplexTree>(simplex_tree, order_policy);
}

template <class SimplexTree>
struct SimplexTreeMorseReferenceResult {
  SimplexTreeComplexView<SimplexTree> view;
  MorseSequence sequence;
  PersistenceDiagram diagram;
  MorseReferenceFrameMetrics frame_metrics;
  MorseReferenceReductionMetrics reduction_metrics;
};

template <class SimplexTree>
SimplexTreeMorseReferenceResult<SimplexTree> compute_simplex_tree_morse_reference_persistence(
    const SimplexTree& simplex_tree,
    MorseSequenceStrategy strategy = MorseSequenceStrategy::Saturated,
    SimplexTreeFiltrationOrder order_policy =
        SimplexTreeFiltrationOrder::CanonicalLexicographic) {
  auto view = make_simplex_tree_complex_view(simplex_tree, order_policy);
  auto input = build_morse_reference_reduction_input(view, strategy);
  MorseReferencePersistenceReducer reducer(view,
                                           input.sequence,
                                           std::move(input.reduction_plan),
                                           std::move(input.annotations));
  auto reduction = reducer.compute_with_metrics();

  return SimplexTreeMorseReferenceResult<SimplexTree>{
      std::move(view),
      std::move(input.sequence),
      std::move(reduction.diagram),
      input.frame_metrics,
      reduction.metrics,
  };
}

}  // namespace morseframes
