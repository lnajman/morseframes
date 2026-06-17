#include <algorithm>
#include <cstdint>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gudhi/Morse_persistence.h>
#include <gudhi/Simplex_tree.h>

#include "morseframes/annotation.hpp"
#include "morseframes/morse_sequence.hpp"
#include "morseframes/reference_persistence.hpp"
#include "morseframes/simplex_tree_builder.hpp"
#include "morseframes/simplex_tree_morse.hpp"
#include "morseframes/standard_persistence.hpp"

namespace {

using GudhiSimplexTree = Gudhi::Simplex_tree<>;
namespace mp = Gudhi::morse_persistence;

constexpr double kEps = 1e-12;

void check(bool condition, const char* expression, const char* file, int line) {
  if (!condition) {
    throw std::logic_error(std::string("Check failed: ") + expression + " at " +
                           file + ":" + std::to_string(line));
  }
}

#define CHECK(expression) check(static_cast<bool>(expression), #expression, __FILE__, __LINE__)

static_assert(morseframes::is_complex_view_v<morseframes::FilteredSimplicialComplex>,
              "FilteredSimplicialComplex must satisfy the Morse complex-view API.");
static_assert(morseframes::is_complex_view_v<morseframes::SimplexTreeComplexView<GudhiSimplexTree>>,
              "SimplexTreeComplexView must satisfy the Morse complex-view API.");

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <= kEps;
}

void insert(GudhiSimplexTree& simplex_tree,
            std::vector<GudhiSimplexTree::Vertex_handle> vertices,
            double filtration) {
  std::sort(vertices.begin(), vertices.end());
  simplex_tree.insert_simplex_and_subfaces(vertices, filtration);
  simplex_tree.clear_filtration();
}


struct SimplexSpec {
  std::vector<GudhiSimplexTree::Vertex_handle> vertices;
  double filtration = 0.0;
};

GudhiSimplexTree make_simplex_tree(std::initializer_list<SimplexSpec> simplices) {
  GudhiSimplexTree simplex_tree;
  for (const auto& simplex : simplices) {
    insert(simplex_tree, simplex.vertices, simplex.filtration);
  }
  simplex_tree.initialize_filtration();
  return simplex_tree;
}

GudhiSimplexTree make_single_vertex() {
  return make_simplex_tree({
      {{0}, 0.0},
  });
}

GudhiSimplexTree make_increasing_edge() {
  return make_simplex_tree({
      {{0}, 0.0},
      {{1}, 0.0},
      {{0, 1}, 1.0},
  });
}

GudhiSimplexTree make_plateau_filled_triangle() {
  return make_simplex_tree({
      {{0, 1, 2}, 0.0},
  });
}

GudhiSimplexTree make_filtered_triangle_with_tail() {
  return make_simplex_tree({
      {{0, 1, 2}, 2.0},
      {{0, 1}, 1.0},
      {{0, 2}, 1.0},
      {{1, 2}, 1.0},
      {{0}, 0.0},
      {{1}, 0.0},
      {{2}, 0.0},
      {{2, 3}, 2.5},
      {{3}, 0.25},
  });
}

GudhiSimplexTree make_two_components_joined_late() {
  return make_simplex_tree({
      {{0}, 0.0},
      {{1}, 0.0},
      {{0, 1}, 1.0},
      {{2}, 0.25},
      {{3}, 0.25},
      {{2, 3}, 0.5},
      {{1, 2}, 2.0},
  });
}

GudhiSimplexTree make_cycle_killed_by_later_triangle() {
  return make_simplex_tree({
      {{0}, 0.0},
      {{1}, 0.0},
      {{2}, 0.0},
      {{0, 1}, 1.0},
      {{0, 2}, 1.0},
      {{1, 2}, 1.0},
      {{0, 1, 2}, 2.0},
  });
}

template <class ComplexView, class MorseSequence>
void validate_sequence_shape(const ComplexView& complex, const MorseSequence& sequence) {
  std::vector<bool> seen(complex.size(), false);
  std::size_t regular_pairs = 0;

  for (const auto& step : sequence.steps()) {
    using StepType = std::decay_t<decltype(step.type)>;
    CHECK(step.sigma < complex.size());
    CHECK(!seen[step.sigma]);

    if (step.type == StepType::Critical) {
      seen[step.sigma] = true;
      continue;
    }

    CHECK(step.tau < complex.size());
    CHECK(!seen[step.tau]);
    CHECK(complex.level(step.sigma) == complex.level(step.tau));

    const auto& boundary = complex.boundary(step.tau);
    CHECK(std::find(boundary.begin(), boundary.end(), step.sigma) != boundary.end());

    seen[step.sigma] = true;
    seen[step.tau] = true;
    ++regular_pairs;
  }

  CHECK(sequence.critical_simplices().size() + 2 * regular_pairs == complex.size());
  CHECK(std::all_of(seen.begin(), seen.end(), [](bool value) { return value; }));
}

template <class ComplexView>
void validate_reference_recurrence(const ComplexView& complex,
                                   const morseframes::MorseSequence& sequence,
                                   const std::vector<morseframes::Annotation>& references) {
  CHECK(references.size() == complex.size());

  for (const morseframes::MorseStep& step : sequence.steps()) {
    if (step.type == morseframes::MorseStepType::Critical) {
      const auto critical_id =
          static_cast<morseframes::CriticalId>(sequence.critical_index(step.sigma));
      CHECK((references.at(step.sigma) == morseframes::Annotation{critical_id}));
      continue;
    }

    CHECK(references.at(step.tau).empty());
    morseframes::Annotation expected;
    for (morseframes::SimplexId face : complex.boundary(step.tau)) {
      if (face != step.sigma) {
        morseframes::xor_annotations_in_place(expected, references.at(face));
      }
    }
    CHECK(references.at(step.sigma) == expected);
  }
}

template <class Diagram>
std::vector<std::tuple<std::uint16_t, double, double>> finite_barcode_key(
    const Diagram& diagram) {
  std::vector<std::tuple<std::uint16_t, double, double>> result;
  result.reserve(diagram.finite_pairs.size());
  for (const auto& interval : diagram.finite_pairs) {
    result.emplace_back(interval.dimension, interval.birth_value, interval.death_value);
  }
  std::sort(result.begin(), result.end());
  return result;
}

template <class Diagram>
auto off_diagonal_pairs_from_diagram(const Diagram& diagram) {
  using Pair = typename std::decay_t<decltype(diagram.finite_pairs)>::value_type;
  std::vector<Pair> result;
  for (const auto& pair : diagram.finite_pairs) {
    if (pair.birth_value < pair.death_value) {
      result.push_back(pair);
    }
  }
  return result;
}

template <class Diagram>
std::vector<std::tuple<std::uint16_t, double, double>> off_diagonal_barcode_key(
    const Diagram& diagram) {
  std::vector<std::tuple<std::uint16_t, double, double>> result;
  for (const auto& interval : off_diagonal_pairs_from_diagram(diagram)) {
    result.emplace_back(interval.dimension, interval.birth_value, interval.death_value);
  }
  std::sort(result.begin(), result.end());
  return result;
}

template <class Diagram>
std::vector<std::pair<std::uint16_t, double>> essential_barcode_key(
    const Diagram& diagram) {
  std::vector<std::pair<std::uint16_t, double>> result;
  result.reserve(diagram.essential.size());
  for (const auto& interval : diagram.essential) {
    result.emplace_back(interval.dimension, interval.birth_value);
  }
  std::sort(result.begin(), result.end());
  return result;
}

template <class Pairs>
std::size_t count_finite_dim(const Pairs& pairs, std::uint16_t dimension) {
  return static_cast<std::size_t>(
      std::count_if(pairs.begin(), pairs.end(), [dimension](const auto& pair) {
        return pair.dimension == dimension;
      }));
}

template <class Diagram>
std::size_t count_essential_dim(const Diagram& diagram, std::uint16_t dimension) {
  return static_cast<std::size_t>(
      std::count_if(diagram.essential.begin(),
                    diagram.essential.end(),
                    [dimension](const auto& interval) {
                      return interval.dimension == dimension;
                    }));
}

template <class Pair>
void assert_off_diagonal_interval(const Pair& pair,
                                  std::uint16_t dimension,
                                  double birth,
                                  double death) {
  CHECK(pair.dimension == dimension);
  CHECK(close(pair.birth_value, birth));
  CHECK(close(pair.death_value, death));
}

template <class ComplexView>
std::vector<std::tuple<int,
                       morseframes::LevelId,
                       std::vector<morseframes::VertexId>,
                       std::vector<morseframes::VertexId>>>
sequence_signature(const ComplexView& complex, const morseframes::MorseSequence& sequence) {
  std::vector<std::tuple<int,
                         morseframes::LevelId,
                         std::vector<morseframes::VertexId>,
                         std::vector<morseframes::VertexId>>>
      result;
  result.reserve(sequence.steps().size());
  for (const auto& step : sequence.steps()) {
    result.emplace_back(step.type == morseframes::MorseStepType::Critical ? 0 : 1,
                        step.level,
                        complex.vertices(step.sigma),
                        step.type == morseframes::MorseStepType::Critical
                            ? std::vector<morseframes::VertexId>{}
                            : complex.vertices(step.tau));
  }
  return result;
}

const std::vector<morseframes::MorseSequenceStrategy>& public_maintainer_strategies() {
  static const std::vector<morseframes::MorseSequenceStrategy> strategies = {
      morseframes::MorseSequenceStrategy::SameLevelReduction,
      morseframes::MorseSequenceStrategy::FMax,
      morseframes::MorseSequenceStrategy::FMin,
      morseframes::MorseSequenceStrategy::PlateauGreedy,
  };
  return strategies;
}

const std::vector<mp::Morse_sequence_strategy>& public_gudhi_strategies() {
  static const std::vector<mp::Morse_sequence_strategy> strategies = {
      mp::Morse_sequence_strategy::SAME_LEVEL_REDUCTION,
      mp::Morse_sequence_strategy::F_MAX,
      mp::Morse_sequence_strategy::F_MIN,
      mp::Morse_sequence_strategy::PLATEAU_GREEDY,
  };
  return strategies;
}

template <class ComplexView>
morseframes::PersistenceDiagram compute_diagram_from_frame(
    const ComplexView& complex,
    const morseframes::MorseReferenceFrame& frame) {
  return morseframes::MorseReferencePersistenceReducer(complex, frame.sequence, frame.references)
      .compute();
}

template <class Result>
void assert_handles_are_mappable(const Result& result) {
  for (const auto& interval : result.diagram.finite_pairs) {
    CHECK(result.view.handle(interval.birth) != GudhiSimplexTree::null_simplex());
    CHECK(result.view.handle(interval.death) != GudhiSimplexTree::null_simplex());
  }
  for (const auto& interval : result.diagram.essential) {
    CHECK(result.view.handle(interval.birth) != GudhiSimplexTree::null_simplex());
  }
}

void validate_strategy_on_simplex_tree(const GudhiSimplexTree& simplex_tree,
                                       morseframes::MorseSequenceStrategy strategy) {
  morseframes::SimplexTreeComplexView<GudhiSimplexTree> view(simplex_tree);
  const auto compact = morseframes::filtered_complex_from_simplex_tree(simplex_tree);
  const auto standard = morseframes::compute_standard_z2_persistence(compact);

  const auto direct_frame = morseframes::build_morse_reference_frame(view, strategy);
  validate_sequence_shape(view, direct_frame.sequence);
  validate_reference_recurrence(view, direct_frame.sequence, direct_frame.references);

  const auto direct_diagram = compute_diagram_from_frame(view, direct_frame);
  CHECK(off_diagonal_barcode_key(direct_diagram) == off_diagonal_barcode_key(standard));
  CHECK(essential_barcode_key(direct_diagram) == essential_barcode_key(standard));

  const auto compact_frame = morseframes::build_morse_reference_frame(compact, strategy);
  validate_sequence_shape(compact, compact_frame.sequence);
  validate_reference_recurrence(compact, compact_frame.sequence, compact_frame.references);

  CHECK(sequence_signature(view, direct_frame.sequence) ==
         sequence_signature(compact, compact_frame.sequence));
}

void validate_public_strategies_on_simplex_tree(const GudhiSimplexTree& simplex_tree) {
  for (const auto& strategy : public_maintainer_strategies()) {
    validate_strategy_on_simplex_tree(simplex_tree, strategy);
  }
}

void validate_public_api_on_simplex_tree(const GudhiSimplexTree& simplex_tree) {
  const auto compact = morseframes::filtered_complex_from_simplex_tree(simplex_tree);
  const auto standard = morseframes::compute_standard_z2_persistence(compact);

  for (const auto& strategy : public_gudhi_strategies()) {
    const auto frame = mp::compute_morse_sequence_and_reference_map(simplex_tree, strategy);
    validate_sequence_shape(frame.view, frame.sequence());

    const auto diagram_from_frame = mp::compute_morse_persistence(simplex_tree, frame);
    CHECK(off_diagonal_barcode_key(diagram_from_frame) == off_diagonal_barcode_key(standard));
    CHECK(essential_barcode_key(diagram_from_frame) == essential_barcode_key(standard));

    const auto result = mp::compute_morse_persistence(simplex_tree, strategy);
    validate_sequence_shape(result.view, result.sequence);
    assert_handles_are_mappable(result);
    CHECK(off_diagonal_barcode_key(result.diagram) == off_diagonal_barcode_key(standard));
    CHECK(essential_barcode_key(result.diagram) == essential_barcode_key(standard));
    CHECK(off_diagonal_barcode_key(diagram_from_frame) ==
          off_diagonal_barcode_key(result.diagram));
    CHECK(essential_barcode_key(diagram_from_frame) == essential_barcode_key(result.diagram));
  }
}

void test_maintainer_single_vertex() {
  const auto simplex_tree = make_single_vertex();
  validate_public_api_on_simplex_tree(simplex_tree);
  validate_public_strategies_on_simplex_tree(simplex_tree);

  const auto result = morseframes::compute_simplex_tree_morse_reference_persistence(
      simplex_tree, morseframes::MorseSequenceStrategy::FMax);
  assert_handles_are_mappable(result);

  const auto finite = off_diagonal_pairs_from_diagram(result.diagram);
  CHECK(finite.empty());
  CHECK(result.diagram.essential.size() == 1);
  CHECK(count_essential_dim(result.diagram, 0) == 1);
  CHECK(close(result.diagram.essential.front().birth_value, 0.0));
}

void test_maintainer_increasing_edge() {
  const auto simplex_tree = make_increasing_edge();
  validate_public_api_on_simplex_tree(simplex_tree);
  validate_public_strategies_on_simplex_tree(simplex_tree);

  CHECK(mp::strategy_from_name("f-max") == mp::Morse_sequence_strategy::f_max);
  CHECK(mp::Morse_sequence_strategy::F_MAX == mp::Morse_sequence_strategy::f_max);
  CHECK(std::string(mp::strategy_name(mp::Morse_sequence_strategy::F_MAX)) == "f-max");
  bool rejected_internal_strategy = false;
  try {
    (void)mp::strategy_from_name("saturated");
  } catch (const std::invalid_argument&) {
    rejected_internal_strategy = true;
  }
  CHECK(rejected_internal_strategy);

  const auto frame = mp::compute_morse_sequence_and_reference_map(
      simplex_tree,
      mp::Morse_sequence_strategy::F_MAX);
  const auto diagram_from_frame = mp::compute_morse_persistence(simplex_tree, frame);

  const auto result = mp::compute_morse_persistence(
      simplex_tree,
      mp::Morse_sequence_strategy::F_MAX);
  assert_handles_are_mappable(result);
  CHECK(off_diagonal_barcode_key(diagram_from_frame) ==
         off_diagonal_barcode_key(result.diagram));
  CHECK(essential_barcode_key(diagram_from_frame) == essential_barcode_key(result.diagram));

  const auto finite = off_diagonal_pairs_from_diagram(result.diagram);
  CHECK(finite.size() == 1);
  assert_off_diagonal_interval(finite.front(), 0, 0.0, 1.0);
  CHECK(result.diagram.essential.size() == 1);
  CHECK(count_essential_dim(result.diagram, 0) == 1);
}

void test_maintainer_plateau_filled_triangle() {
  const auto simplex_tree = make_plateau_filled_triangle();
  validate_public_api_on_simplex_tree(simplex_tree);
  validate_public_strategies_on_simplex_tree(simplex_tree);

  const auto result = morseframes::compute_simplex_tree_morse_reference_persistence(
      simplex_tree, morseframes::MorseSequenceStrategy::FMax);
  assert_handles_are_mappable(result);

  const auto finite = off_diagonal_pairs_from_diagram(result.diagram);
  CHECK(finite.empty());
  CHECK(result.diagram.essential.size() == 1);
  CHECK(count_essential_dim(result.diagram, 0) == 1);
}

void test_maintainer_triangle_with_tail() {
  const auto simplex_tree = make_filtered_triangle_with_tail();
  validate_public_api_on_simplex_tree(simplex_tree);
  validate_public_strategies_on_simplex_tree(simplex_tree);

  const auto result = morseframes::compute_simplex_tree_morse_reference_persistence(
      simplex_tree, morseframes::MorseSequenceStrategy::FMax);
  assert_handles_are_mappable(result);

  const auto finite = off_diagonal_pairs_from_diagram(result.diagram);
  CHECK(finite.size() == 4);
  CHECK(count_finite_dim(finite, 0) == 3);
  CHECK(count_finite_dim(finite, 1) == 1);
  CHECK(count_essential_dim(result.diagram, 0) == 1);
}

void test_maintainer_two_components_joined_late() {
  const auto simplex_tree = make_two_components_joined_late();
  validate_public_api_on_simplex_tree(simplex_tree);
  validate_public_strategies_on_simplex_tree(simplex_tree);

  const auto result = morseframes::compute_simplex_tree_morse_reference_persistence(
      simplex_tree, morseframes::MorseSequenceStrategy::FMax);
  assert_handles_are_mappable(result);

  const auto finite = off_diagonal_pairs_from_diagram(result.diagram);
  CHECK(finite.size() == 3);
  CHECK(count_finite_dim(finite, 0) == 3);
  CHECK(count_finite_dim(finite, 1) == 0);
  CHECK(count_essential_dim(result.diagram, 0) == 1);
}

void test_maintainer_cycle_killed_by_later_triangle() {
  const auto simplex_tree = make_cycle_killed_by_later_triangle();
  validate_public_api_on_simplex_tree(simplex_tree);
  validate_public_strategies_on_simplex_tree(simplex_tree);

  const auto result = morseframes::compute_simplex_tree_morse_reference_persistence(
      simplex_tree, morseframes::MorseSequenceStrategy::FMax);
  assert_handles_are_mappable(result);

  const auto finite = off_diagonal_pairs_from_diagram(result.diagram);
  CHECK(finite.size() == 3);
  CHECK(count_finite_dim(finite, 0) == 2);
  CHECK(count_finite_dim(finite, 1) == 1);
  CHECK(count_essential_dim(result.diagram, 0) == 1);
  CHECK(count_essential_dim(result.diagram, 1) == 0);
}

void test_real_gudhi_simplex_tree_view() {
  auto simplex_tree = make_filtered_triangle_with_tail();
  CHECK(simplex_tree.num_simplices() == 9);

  morseframes::SimplexTreeComplexView<GudhiSimplexTree> view(simplex_tree);
  CHECK(view.size() == simplex_tree.num_simplices());
  CHECK(view.num_levels() == 5);

  const auto triangle = view.find_simplex({0, 1, 2});
  CHECK(triangle != morseframes::kInvalidSimplex);
  CHECK(close(view.filtration(triangle), 2.0));
  CHECK(view.handle(triangle) != GudhiSimplexTree::null_simplex());
  const auto& triangle_boundary = view.boundary(triangle);
  CHECK(triangle_boundary.size() == 3);
  CHECK(view.vertices(triangle_boundary[0]) == std::vector<morseframes::VertexId>({1, 2}));
  CHECK(view.vertices(triangle_boundary[1]) == std::vector<morseframes::VertexId>({0, 2}));
  CHECK(view.vertices(triangle_boundary[2]) == std::vector<morseframes::VertexId>({0, 1}));

  auto copied_view = view;
  CHECK(copied_view.vertices(triangle) == std::vector<morseframes::VertexId>({0, 1, 2}));
  auto moved_view = std::move(copied_view);
  CHECK(moved_view.vertices(triangle) == std::vector<morseframes::VertexId>({0, 1, 2}));

  const auto sequence = morseframes::FSequenceBuilder(view).build_f_max();
  validate_sequence_shape(view, sequence);

  const auto references =
      morseframes::MorseReferenceComputer(view, sequence).compute_full_references();
  validate_reference_recurrence(view, sequence, references);

  const auto frame = morseframes::MorseReferenceFrameBuilder(view).build_f_max();
  CHECK(frame.sequence.critical_simplices() == sequence.critical_simplices());
  CHECK(frame.references == references);

  const auto direct_diagram =
      morseframes::MorseReferencePersistenceReducer(view, frame.sequence, frame.references).compute();
  const auto compact_complex = morseframes::filtered_complex_from_simplex_tree(simplex_tree);
  const auto standard = morseframes::compute_standard_z2_persistence(compact_complex);
  CHECK(finite_barcode_key(direct_diagram) == finite_barcode_key(standard));
  CHECK(essential_barcode_key(direct_diagram) == essential_barcode_key(standard));

  const auto api_result = morseframes::compute_simplex_tree_morse_reference_persistence(
      simplex_tree, morseframes::MorseSequenceStrategy::FMax);
  CHECK(finite_barcode_key(api_result.diagram) == finite_barcode_key(standard));
  CHECK(essential_barcode_key(api_result.diagram) == essential_barcode_key(standard));
  CHECK(!api_result.sequence.critical_simplices().empty());
  CHECK(api_result.view.handle(api_result.sequence.critical_simplices().front()) !=
         GudhiSimplexTree::null_simplex());
}

void test_real_gudhi_import_to_compact_complex() {
  auto simplex_tree = make_filtered_triangle_with_tail();
  auto complex = morseframes::filtered_complex_from_simplex_tree(simplex_tree);

  CHECK(complex.size() == simplex_tree.num_simplices());
  const auto triangle = complex.find_simplex({0, 1, 2});
  CHECK(triangle != morseframes::kInvalidSimplex);
  CHECK(close(complex.filtration(triangle), 2.0));

  const auto frame = morseframes::MorseReferenceFrameBuilder(complex).build_f_max();
  const auto diagram =
      morseframes::MorseReferencePersistenceReducer(complex, frame.sequence, frame.references)
          .compute();
  const auto standard = morseframes::compute_standard_z2_persistence(complex);

  CHECK(morseframes::off_diagonal_pairs(diagram).size() ==
         morseframes::off_diagonal_pairs(standard).size());
  CHECK(diagram.essential.size() == standard.essential.size());
}

}  // namespace

int main() {
  test_maintainer_single_vertex();
  test_maintainer_increasing_edge();
  test_maintainer_plateau_filled_triangle();
  test_maintainer_triangle_with_tail();
  test_maintainer_two_components_joined_late();
  test_maintainer_cycle_killed_by_later_triangle();
  test_real_gudhi_simplex_tree_view();
  test_real_gudhi_import_to_compact_complex();

  std::cout << "GUDHI Simplex_tree view tests passed.\n";
  return 0;
}
