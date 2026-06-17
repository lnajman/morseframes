#include <algorithm>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <gudhi/Morse_persistence.h>
#include <gudhi/Simplex_tree.h>

namespace {

using Simplex_tree = Gudhi::Simplex_tree<>;
namespace mp = Gudhi::morse_persistence;

constexpr double kEps = 1e-12;

using Finite_bar = std::tuple<std::uint16_t, double, double>;
using Essential_bar = std::pair<std::uint16_t, double>;

void check(bool condition, const char* expression, const char* file, int line) {
  if (!condition) {
    throw std::logic_error(std::string("Check failed: ") + expression + " at " +
                           file + ":" + std::to_string(line));
  }
}

#define CHECK(expression) check(static_cast<bool>(expression), #expression, __FILE__, __LINE__)

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <= kEps;
}

void insert(Simplex_tree& simplex_tree,
            std::vector<Simplex_tree::Vertex_handle> vertices,
            double filtration) {
  std::sort(vertices.begin(), vertices.end());
  simplex_tree.insert_simplex_and_subfaces(vertices, filtration);
  simplex_tree.clear_filtration();
}

struct Simplex_spec {
  std::vector<Simplex_tree::Vertex_handle> vertices;
  double filtration = 0.0;
};

Simplex_tree make_simplex_tree(std::initializer_list<Simplex_spec> simplices) {
  Simplex_tree simplex_tree;
  for (const auto& simplex : simplices) {
    insert(simplex_tree, simplex.vertices, simplex.filtration);
  }
  simplex_tree.initialize_filtration();
  return simplex_tree;
}

Simplex_tree make_single_vertex() {
  return make_simplex_tree({
      {{0}, 0.0},
  });
}

Simplex_tree make_increasing_edge() {
  return make_simplex_tree({
      {{0}, 0.0},
      {{1}, 0.0},
      {{0, 1}, 1.0},
  });
}

Simplex_tree make_plateau_filled_triangle() {
  return make_simplex_tree({
      {{0, 1, 2}, 0.0},
  });
}

Simplex_tree make_filtered_triangle_with_tail() {
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

Simplex_tree make_two_components_joined_late() {
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

Simplex_tree make_cycle_killed_by_later_triangle() {
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

const std::vector<mp::Morse_sequence_strategy>& public_strategies() {
  static const std::vector<mp::Morse_sequence_strategy> strategies = {
      mp::Morse_sequence_strategy::SAME_LEVEL_REDUCTION,
      mp::Morse_sequence_strategy::F_MAX,
      mp::Morse_sequence_strategy::F_MIN,
      mp::Morse_sequence_strategy::PLATEAU_GREEDY,
  };
  return strategies;
}

template <class ComplexView, class MorseSequence>
void validate_sequence_shape(const ComplexView& complex, const MorseSequence& sequence) {
  std::vector<bool> seen(complex.size(), false);
  std::size_t regular_pairs = 0;

  for (const auto& step : sequence.steps()) {
    using Step_type = std::decay_t<decltype(step.type)>;
    CHECK(step.sigma < complex.size());
    CHECK(!seen[step.sigma]);

    if (step.type == Step_type::Critical) {
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

template <class Diagram>
std::vector<Finite_bar> off_diagonal_barcode_key(const Diagram& diagram) {
  std::vector<Finite_bar> result;
  for (const auto& interval : mp::off_diagonal_pairs(diagram)) {
    result.emplace_back(interval.dimension, interval.birth_value, interval.death_value);
  }
  std::sort(result.begin(), result.end());
  return result;
}

template <class Diagram>
std::vector<Essential_bar> essential_barcode_key(const Diagram& diagram) {
  std::vector<Essential_bar> result;
  result.reserve(diagram.essential.size());
  for (const auto& interval : diagram.essential) {
    result.emplace_back(interval.dimension, interval.birth_value);
  }
  std::sort(result.begin(), result.end());
  return result;
}

void check_finite_barcode(const std::vector<Finite_bar>& actual,
                          std::vector<Finite_bar> expected) {
  std::sort(expected.begin(), expected.end());
  CHECK(actual.size() == expected.size());
  for (std::size_t index = 0; index < expected.size(); ++index) {
    CHECK(std::get<0>(actual[index]) == std::get<0>(expected[index]));
    CHECK(close(std::get<1>(actual[index]), std::get<1>(expected[index])));
    CHECK(close(std::get<2>(actual[index]), std::get<2>(expected[index])));
  }
}

void check_essential_barcode(const std::vector<Essential_bar>& actual,
                             std::vector<Essential_bar> expected) {
  std::sort(expected.begin(), expected.end());
  CHECK(actual.size() == expected.size());
  for (std::size_t index = 0; index < expected.size(); ++index) {
    CHECK(actual[index].first == expected[index].first);
    CHECK(close(actual[index].second, expected[index].second));
  }
}

template <class Result>
void assert_handles_are_mappable(const Result& result) {
  for (const auto& interval : result.diagram.finite_pairs) {
    CHECK(result.simplex_tree_handle(interval.birth) != Simplex_tree::null_simplex());
    CHECK(result.simplex_tree_handle(interval.death) != Simplex_tree::null_simplex());
  }
  for (const auto& interval : result.diagram.essential) {
    CHECK(result.simplex_tree_handle(interval.birth) != Simplex_tree::null_simplex());
  }
}

void validate_case(const Simplex_tree& simplex_tree,
                   const std::vector<Finite_bar>& expected_finite,
                   const std::vector<Essential_bar>& expected_essential) {
  for (const auto& strategy : public_strategies()) {
    const auto frame = mp::compute_morse_sequence_and_reference_map(simplex_tree, strategy);
    validate_sequence_shape(frame.view, frame.sequence());

    const auto diagram_from_frame = mp::compute_morse_persistence(simplex_tree, frame);
    check_finite_barcode(off_diagonal_barcode_key(diagram_from_frame), expected_finite);
    check_essential_barcode(essential_barcode_key(diagram_from_frame), expected_essential);

    const auto result = mp::compute_morse_persistence(simplex_tree, strategy);
    validate_sequence_shape(result.view, result.sequence);
    assert_handles_are_mappable(result);
    check_finite_barcode(off_diagonal_barcode_key(result.diagram), expected_finite);
    check_essential_barcode(essential_barcode_key(result.diagram), expected_essential);
    CHECK(off_diagonal_barcode_key(result.diagram) ==
          off_diagonal_barcode_key(diagram_from_frame));
    CHECK(essential_barcode_key(result.diagram) == essential_barcode_key(diagram_from_frame));
  }
}

void test_strategy_names() {
  CHECK(mp::strategy_from_name("same-level-reduction") ==
        mp::Morse_sequence_strategy::SAME_LEVEL_REDUCTION);
  CHECK(mp::strategy_from_name("f-max") == mp::Morse_sequence_strategy::F_MAX);
  CHECK(mp::strategy_from_name("f-min") == mp::Morse_sequence_strategy::F_MIN);
  CHECK(mp::strategy_from_name("plateau-greedy") ==
        mp::Morse_sequence_strategy::PLATEAU_GREEDY);
  CHECK(std::string(mp::strategy_name(mp::Morse_sequence_strategy::F_MAX)) == "f-max");

  bool rejected_internal_strategy = false;
  try {
    (void)mp::strategy_from_name("saturated");
  } catch (const std::invalid_argument&) {
    rejected_internal_strategy = true;
  }
  CHECK(rejected_internal_strategy);
}

void test_single_vertex() {
  validate_case(make_single_vertex(), {}, {{0, 0.0}});
}

void test_increasing_edge() {
  validate_case(make_increasing_edge(), {{0, 0.0, 1.0}}, {{0, 0.0}});
}

void test_plateau_filled_triangle() {
  validate_case(make_plateau_filled_triangle(), {}, {{0, 0.0}});
}

void test_filtered_triangle_with_tail() {
  validate_case(make_filtered_triangle_with_tail(),
                {{0, 0.0, 1.0},
                 {0, 0.0, 1.0},
                 {0, 0.25, 2.5},
                 {1, 1.0, 2.0}},
                {{0, 0.0}});
}

void test_two_components_joined_late() {
  validate_case(make_two_components_joined_late(),
                {{0, 0.0, 1.0},
                 {0, 0.25, 0.5},
                 {0, 0.25, 2.0}},
                {{0, 0.0}});
}

void test_cycle_killed_by_later_triangle() {
  validate_case(make_cycle_killed_by_later_triangle(),
                {{0, 0.0, 1.0},
                 {0, 0.0, 1.0},
                 {1, 1.0, 2.0}},
                {{0, 0.0}});
}

}  // namespace

int main() {
  test_strategy_names();
  test_single_vertex();
  test_increasing_edge();
  test_plateau_filled_triangle();
  test_filtered_triangle_with_tail();
  test_two_components_joined_late();
  test_cycle_killed_by_later_triangle();

  std::cout << "GUDHI Morse persistence public API tests passed.\n";
  return 0;
}
