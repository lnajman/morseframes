#include <algorithm>
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "morseframes/coreference_persistence.hpp"
#include "morseframes/debug_checks.hpp"
#include "morseframes/field_annotation_store.hpp"
#include "morseframes/filtered_complex.hpp"
#include "morseframes/instrumentation.hpp"
#include "morseframes/inverse_annotation_store.hpp"
#include "morseframes/morse_sequence.hpp"
#include "morseframes/reference_persistence.hpp"
#include "morseframes/simplex_tree_builder.hpp"
#include "morseframes/standard_persistence.hpp"

namespace {

using morseframes::FilteredSimplicialComplex;
using morseframes::FSequenceBuilder;
using morseframes::PersistenceDiagram;

constexpr double kEps = 1e-12;

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <= kEps;
}

bool field_annotation_equals(
    const morseframes::FieldAnnotation& annotation,
    std::initializer_list<std::pair<morseframes::CriticalId, std::uint32_t>> expected) {
  if (annotation.size() != expected.size()) {
    return false;
  }
  std::size_t index = 0;
  for (const auto& entry : expected) {
    if (annotation[index].label != entry.first ||
        annotation[index].coefficient != entry.second) {
      return false;
    }
    ++index;
  }
  return true;
}

using FiniteKey = std::tuple<std::uint16_t, double, double>;
using EssentialKey = std::tuple<std::uint16_t, double>;

std::vector<FiniteKey> finite_barcode(const PersistenceDiagram& diagram) {
  std::vector<FiniteKey> result;
  for (const auto& pair : morseframes::off_diagonal_pairs(diagram)) {
    result.emplace_back(pair.dimension, pair.birth_value, pair.death_value);
  }
  std::sort(result.begin(), result.end());
  return result;
}

std::vector<EssentialKey> essential_barcode(const PersistenceDiagram& diagram) {
  std::vector<EssentialKey> result;
  for (const auto& interval : diagram.essential) {
    result.emplace_back(interval.dimension, interval.birth_value);
  }
  std::sort(result.begin(), result.end());
  return result;
}

void assert_same_barcode(const PersistenceDiagram& lhs, const PersistenceDiagram& rhs) {
  assert(finite_barcode(lhs) == finite_barcode(rhs));
  assert(essential_barcode(lhs) == essential_barcode(rhs));
}

void add_simplex(FilteredSimplicialComplex& complex,
                 std::initializer_list<morseframes::VertexId> vertices,
                 double filtration) {
  complex.add_simplex(std::vector<morseframes::VertexId>(vertices), filtration);
}

void add_weighted_closure(FilteredSimplicialComplex& complex,
                          const std::vector<morseframes::VertexId>& facet,
                          const std::vector<double>& vertex_values,
                          double dimension_offset = 0.0) {
  const std::uint32_t mask_limit = 1u << facet.size();
  for (std::uint32_t mask = 1; mask < mask_limit; ++mask) {
    std::vector<morseframes::VertexId> simplex;
    double filtration = -1.0;
    for (std::size_t bit = 0; bit < facet.size(); ++bit) {
      if ((mask & (1u << bit)) == 0) {
        continue;
      }
      const auto vertex = facet[bit];
      simplex.push_back(vertex);
      filtration = std::max(filtration, vertex_values.at(vertex));
    }
    filtration += dimension_offset * static_cast<double>(simplex.size() - 1);
    complex.add_simplex(simplex, filtration);
  }
}

PersistenceDiagram run_reference(FilteredSimplicialComplex& complex) {
  complex.finalize();
  auto sequence = FSequenceBuilder(complex).build_saturated();
  morseframes::validate_morse_sequence(complex, sequence);

  auto references = morseframes::MorseReferenceComputer(complex, sequence).compute_full_references();
  morseframes::validate_reference_invariants(complex, sequence, references);
  auto morse_diagram =
      morseframes::MorseReferencePersistenceReducer(complex, sequence, references).compute();

  auto frame = morseframes::MorseReferenceFrameBuilder(complex).build_saturated();
  morseframes::validate_morse_sequence(complex, frame.sequence);
  morseframes::validate_reference_invariants(complex, frame.sequence, frame.references);
  assert(frame.sequence.critical_simplices() == sequence.critical_simplices());
  assert(frame.references == references);
  auto fused_diagram =
      morseframes::MorseReferencePersistenceReducer(complex, frame.sequence, frame.references).compute();

  auto compact_input =
      morseframes::MorseReferenceFrameBuilder(complex).build_saturated_reduction_input();
  assert(compact_input.sequence.critical_simplices() == sequence.critical_simplices());
  const auto expected_plan =
      morseframes::build_reference_reduction_plan(complex, compact_input.sequence, references);
  auto actual_working_set = compact_input.reduction_plan.working_set;
  auto expected_working_set = expected_plan.working_set;
  std::sort(actual_working_set.begin(), actual_working_set.end());
  std::sort(expected_working_set.begin(), expected_working_set.end());
  assert(actual_working_set == expected_working_set);
  assert(compact_input.reduction_plan.boundary_candidates.size() ==
         expected_plan.boundary_candidates.size());
  for (std::size_t index = 0; index < expected_plan.boundary_candidates.size(); ++index) {
    assert(compact_input.reduction_plan.boundary_candidates[index].critical_id ==
           expected_plan.boundary_candidates[index].critical_id);
    assert(compact_input.reduction_plan.boundary_candidates[index].simplex ==
           expected_plan.boundary_candidates[index].simplex);
  }
  assert(compact_input.reduction_plan.zero_boundary_critical_ids ==
         expected_plan.zero_boundary_critical_ids);
  assert(compact_input.reduction_plan.boundary_face_scans ==
         expected_plan.boundary_face_scans);
  assert(compact_input.reduction_plan.zero_boundary_skipped_faces ==
         expected_plan.zero_boundary_skipped_faces);
  assert(compact_input.annotations.size() == compact_input.reduction_plan.working_set.size());
  for (std::size_t index = 0; index < compact_input.annotations.size(); ++index) {
    assert(compact_input.annotations[index] ==
           references.at(compact_input.reduction_plan.working_set[index]));
  }
  auto compact_diagram =
      morseframes::MorseReferencePersistenceReducer(complex,
                                              compact_input.sequence,
                                              std::move(compact_input.reduction_plan),
                                              std::move(compact_input.annotations))
          .compute();

  auto deferred_input =
      morseframes::MorseReferenceFrameBuilder(
          complex, false, morseframes::ReferenceFrameReleasePolicy::Deferred)
          .build_saturated_reduction_input();
  assert(deferred_input.sequence.critical_simplices() == sequence.critical_simplices());
  assert(deferred_input.reduction_plan.working_set == expected_plan.working_set);
  assert(deferred_input.reduction_plan.boundary_candidates.size() ==
         expected_plan.boundary_candidates.size());
  assert(deferred_input.reduction_plan.zero_boundary_critical_ids ==
         expected_plan.zero_boundary_critical_ids);
  assert(deferred_input.frame_metrics.released_annotations == 0);
  assert(deferred_input.frame_metrics.released_total_annotation_size == 0);
  assert(deferred_input.annotations.size() == deferred_input.reduction_plan.working_set.size());
  for (std::size_t index = 0; index < deferred_input.annotations.size(); ++index) {
    assert(deferred_input.annotations[index] ==
           references.at(deferred_input.reduction_plan.working_set[index]));
  }
  auto deferred_diagram =
      morseframes::MorseReferencePersistenceReducer(complex,
                                              deferred_input.sequence,
                                              std::move(deferred_input.reduction_plan),
                                              std::move(deferred_input.annotations))
          .compute();
  assert_same_barcode(compact_diagram, deferred_diagram);

  auto coreduction_sequence = FSequenceBuilder(complex).build_same_level_reduction();
  morseframes::validate_morse_sequence(complex, coreduction_sequence);
  auto coreduction_references =
      morseframes::MorseReferenceComputer(complex, coreduction_sequence).compute_full_references();
  morseframes::validate_reference_invariants(complex, coreduction_sequence, coreduction_references);
  auto coreduction_diagram =
      morseframes::MorseReferencePersistenceReducer(
          complex, coreduction_sequence, coreduction_references)
          .compute();

  std::vector<PersistenceDiagram> flooding_diagrams;
  auto check_flooding_sequence = [&](morseframes::MorseSequence sequence) {
    morseframes::validate_morse_sequence(complex, sequence);
    auto references =
        morseframes::MorseReferenceComputer(complex, sequence).compute_full_references();
    morseframes::validate_reference_invariants(complex, sequence, references);
    flooding_diagrams.push_back(
        morseframes::MorseReferencePersistenceReducer(complex, sequence, references).compute());
  };
  check_flooding_sequence(FSequenceBuilder(complex).build_flooding_max());
  check_flooding_sequence(FSequenceBuilder(complex).build_flooding_min());
  check_flooding_sequence(FSequenceBuilder(complex).build_flooding_minmax());
  check_flooding_sequence(FSequenceBuilder(complex).build_flooding_maxmin());
  check_flooding_sequence(FSequenceBuilder(complex).build_f_max());
  check_flooding_sequence(FSequenceBuilder(complex).build_f_min());

  auto coreduction_coreference_frame =
      morseframes::MorseCoreferenceFrameBuilder(complex).build_same_level_reduction();
  morseframes::validate_morse_sequence(complex, coreduction_coreference_frame.sequence);
  assert(coreduction_coreference_frame.sequence.critical_simplices() ==
         coreduction_sequence.critical_simplices());
  auto expected_coreduction_coreferences =
      morseframes::MorseCoreferenceComputer(
          complex, coreduction_coreference_frame.sequence)
          .compute_full_coreferences();
  assert(coreduction_coreference_frame.coreferences ==
         expected_coreduction_coreferences);
  morseframes::validate_coreference_invariants(
      complex,
      coreduction_coreference_frame.sequence,
      coreduction_coreference_frame.coreferences);
  auto coreduction_coreference_diagram =
      morseframes::MorseCoreferencePersistenceReducer(
          complex,
          coreduction_coreference_frame.sequence,
          coreduction_coreference_frame.coreferences)
          .compute();

  auto coreferences =
      morseframes::MorseCoreferenceComputer(complex, sequence).compute_full_coreferences();
  morseframes::validate_coreference_invariants(complex, sequence, coreferences);
  auto coreference_diagram =
      morseframes::MorseCoreferencePersistenceReducer(complex, sequence, coreferences).compute();

  auto standard_diagram = morseframes::compute_standard_z2_persistence(complex);
  morseframes::validate_persistence_diagram(morse_diagram);
  morseframes::validate_persistence_diagram(fused_diagram);
  morseframes::validate_persistence_diagram(compact_diagram);
  morseframes::validate_persistence_diagram(coreduction_diagram);
  for (const auto& diagram : flooding_diagrams) {
    morseframes::validate_persistence_diagram(diagram);
  }
  morseframes::validate_persistence_diagram(coreduction_coreference_diagram);
  morseframes::validate_persistence_diagram(coreference_diagram);
  morseframes::validate_persistence_diagram(standard_diagram);

  assert_same_barcode(morse_diagram, standard_diagram);
  assert_same_barcode(fused_diagram, standard_diagram);
  assert_same_barcode(compact_diagram, standard_diagram);
  assert_same_barcode(coreduction_diagram, standard_diagram);
  for (const auto& diagram : flooding_diagrams) {
    assert_same_barcode(diagram, standard_diagram);
  }
  assert_same_barcode(coreduction_coreference_diagram, standard_diagram);
  assert_same_barcode(coreference_diagram, standard_diagram);
  return morse_diagram;
}

void assert_field_reference_matches_standard(
    const FilteredSimplicialComplex& complex,
    const morseframes::MorseSequence& sequence,
    std::uint32_t modulus) {
  morseframes::validate_morse_sequence(complex, sequence);
  auto references =
      morseframes::MorseFieldReferenceComputer(complex, sequence, modulus).compute_full_references();
  auto full_morse_diagram =
      morseframes::MorseFieldReferencePersistenceReducer(
          complex, sequence, references, modulus)
          .compute();
  auto reduction_plan = morseframes::build_reference_reduction_plan(complex, sequence, references);
  std::vector<morseframes::FieldAnnotation> compact_annotations;
  compact_annotations.reserve(reduction_plan.working_set.size());
  for (morseframes::SimplexId simplex : reduction_plan.working_set) {
    compact_annotations.push_back(references.at(simplex));
  }
  auto compact_morse_diagram =
      morseframes::MorseCompactFieldReferencePersistenceReducer(
          complex, sequence, std::move(reduction_plan), std::move(compact_annotations), modulus)
          .compute();
  auto wrapper_morse_diagram =
      morseframes::compute_morse_reference_prime_field_persistence(complex, sequence, modulus);
  auto standard_diagram =
      morseframes::compute_standard_prime_field_persistence(complex, modulus);
  assert_same_barcode(full_morse_diagram, standard_diagram);
  assert_same_barcode(compact_morse_diagram, standard_diagram);
  assert_same_barcode(wrapper_morse_diagram, standard_diagram);
}

void assert_field_coreference_matches_standard(
    const FilteredSimplicialComplex& complex,
    const morseframes::MorseSequence& sequence,
    std::uint32_t modulus) {
  morseframes::validate_morse_sequence(complex, sequence);
  auto coreferences =
      morseframes::MorseFieldCoreferenceComputer(complex, sequence, modulus)
          .compute_full_coreferences();
  auto full_morse_diagram =
      morseframes::MorseFieldCoreferencePersistenceReducer(
          complex, sequence, coreferences, modulus)
          .compute();
  auto working_set = morseframes::coreference_working_set(complex, sequence);
  std::vector<morseframes::FieldAnnotation> compact_annotations;
  compact_annotations.reserve(working_set.size());
  for (morseframes::SimplexId simplex : working_set) {
    compact_annotations.push_back(coreferences.at(simplex));
  }
  auto compact_morse_diagram =
      morseframes::MorseCompactFieldCoreferencePersistenceReducer(
          complex, sequence, std::move(working_set), std::move(compact_annotations), modulus)
          .compute();
  auto wrapper_morse_diagram =
      morseframes::compute_morse_coreference_prime_field_persistence(complex, sequence, modulus);
  auto standard_diagram =
      morseframes::compute_standard_prime_field_persistence(complex, modulus);
  assert_same_barcode(full_morse_diagram, standard_diagram);
  assert_same_barcode(compact_morse_diagram, standard_diagram);
  assert_same_barcode(wrapper_morse_diagram, standard_diagram);
}

std::size_t count_essential_dim(const PersistenceDiagram& diagram, std::uint16_t dim) {
  std::size_t count = 0;
  for (const auto& interval : diagram.essential) {
    if (interval.dimension == dim) {
      ++count;
    }
  }
  return count;
}

std::size_t count_finite_dim(const std::vector<morseframes::PersistencePair>& pairs, std::uint16_t dim) {
  std::size_t count = 0;
  for (const auto& pair : pairs) {
    if (pair.dimension == dim) {
      ++count;
    }
  }
  return count;
}

void test_boundary_and_coboundary() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 0.0);
  add_simplex(complex, {1}, 0.0);
  add_simplex(complex, {0, 1}, 1.0);
  complex.finalize();

  const auto edge = complex.find_simplex({0, 1});
  const auto v0 = complex.find_simplex({0});
  const auto v1 = complex.find_simplex({1});

  assert(complex.boundary(edge).size() == 2);
  assert(complex.coboundary(v0).size() == 1);
  assert(complex.coboundary(v1).size() == 1);
  assert(complex.coboundary(v0).front() == edge);
  assert(complex.coboundary(v1).front() == edge);
}

void test_inverse_annotation_store() {
  std::vector<morseframes::Annotation> annotations = {
      morseframes::Annotation{0, 2},
      morseframes::Annotation{1},
      morseframes::Annotation{0, 1},
  };
  morseframes::InverseAnnotationStore store(std::move(annotations), 3);

  store.remove_label_from_all(0);
  assert((store.annotation(0) == morseframes::Annotation{2}));
  assert((store.annotation(1) == morseframes::Annotation{1}));
  assert((store.annotation(2) == morseframes::Annotation{1}));

  store.xor_into_all_containing(1, morseframes::Annotation{1, 2});
  assert((store.annotation(0) == morseframes::Annotation{2}));
  assert((store.annotation(1) == morseframes::Annotation{2}));
  assert((store.annotation(2) == morseframes::Annotation{2}));

  store.xor_into_all_containing(2, morseframes::Annotation{2});
  assert(store.annotation(0).empty());
  assert(store.annotation(1).empty());
  assert(store.annotation(2).empty());
  assert(store.metrics().initial_nonempty_annotations == 3);
  assert(store.metrics().initial_total_annotation_size == 5);
  assert(store.metrics().initial_max_annotation_size == 2);
  assert(store.metrics().initial_inverse_list_entries == 5);
  assert(store.metrics().remove_candidate_scans == 2);
  assert(store.metrics().remove_applied == 2);
  assert(store.metrics().remove_total_annotation_size == 4);
  assert(store.metrics().remove_max_annotation_size == 2);
  assert(store.metrics().xor_candidate_scans == 5);
  assert(store.metrics().xor_applied == 5);
  assert(store.metrics().xor_changed_labels == 7);
  assert(store.metrics().xor_total_input_size == 12);
  assert(store.metrics().xor_total_output_size == 2);
  assert(store.metrics().xor_max_input_size == 3);
  assert(store.metrics().xor_max_output_size == 1);
  assert(store.metrics().xor_inserted_labels == 2);
  assert(store.metrics().xor_removed_labels == 5);
  assert(store.metrics().inverse_list_appends == 2);

  std::vector<morseframes::Annotation> restricted_annotations = {
      morseframes::Annotation{0},
      morseframes::Annotation{1},
      morseframes::Annotation{0, 1},
  };
  morseframes::InverseAnnotationStore restricted_store(restricted_annotations, {0, 2}, 2);
  assert((restricted_store.annotation(0) == morseframes::Annotation{0}));
  assert((restricted_store.annotation(2) == morseframes::Annotation{0, 1}));

  bool rejected = false;
  try {
    (void)restricted_store.annotation(1);
  } catch (const std::out_of_range&) {
    rejected = true;
  }
  assert(rejected);
}

void test_field_annotation_store() {
  std::vector<morseframes::FieldAnnotation> annotations = {
      morseframes::FieldAnnotation{
          morseframes::FieldAnnotationEntry{0, 1}, morseframes::FieldAnnotationEntry{2, 2}},
      morseframes::FieldAnnotation{morseframes::FieldAnnotationEntry{1, 1}},
      morseframes::FieldAnnotation{
          morseframes::FieldAnnotationEntry{0, 2}, morseframes::FieldAnnotationEntry{1, 1}},
  };
  morseframes::FieldAnnotationStore store(std::move(annotations),
                                    std::vector<morseframes::SimplexId>{0, 1, 2},
                                    3,
                                    3,
                                    3);

  store.remove_label_from_all(0);
  assert(field_annotation_equals(store.annotation(0), {{2, 2}}));
  assert(field_annotation_equals(store.annotation(1), {{1, 1}}));
  assert(field_annotation_equals(store.annotation(2), {{1, 1}}));

  store.eliminate_pivot(
      1,
      morseframes::FieldAnnotation{
          morseframes::FieldAnnotationEntry{1, 1}, morseframes::FieldAnnotationEntry{2, 1}},
      1);
  assert(field_annotation_equals(store.annotation(0), {{2, 2}}));
  assert(field_annotation_equals(store.annotation(1), {{2, 2}}));
  assert(field_annotation_equals(store.annotation(2), {{2, 2}}));

  store.eliminate_pivot(
      2, morseframes::FieldAnnotation{morseframes::FieldAnnotationEntry{2, 2}}, 2);
  assert(store.annotation(0).empty());
  assert(store.annotation(1).empty());
  assert(store.annotation(2).empty());
  assert(store.metrics().initial_nonempty_annotations == 3);
  assert(store.metrics().initial_total_annotation_size == 5);
  assert(store.metrics().initial_max_annotation_size == 2);
  assert(store.metrics().initial_inverse_list_entries == 5);
  assert(store.metrics().remove_candidate_scans == 2);
  assert(store.metrics().remove_applied == 2);
  assert(store.metrics().remove_total_annotation_size == 4);
  assert(store.metrics().remove_max_annotation_size == 2);
  assert(store.metrics().xor_candidate_scans == 5);
  assert(store.metrics().xor_applied == 5);
  assert(store.metrics().xor_changed_labels == 7);
  assert(store.metrics().xor_total_input_size == 12);
  assert(store.metrics().xor_total_output_size == 2);
  assert(store.metrics().xor_max_input_size == 3);
  assert(store.metrics().xor_max_output_size == 1);
  assert(store.metrics().xor_inserted_labels == 2);
  assert(store.metrics().xor_removed_labels == 5);
  assert(store.metrics().inverse_list_appends == 2);

  std::vector<morseframes::FieldAnnotation> restricted_annotations = {
      morseframes::FieldAnnotation{morseframes::FieldAnnotationEntry{0, 1}},
      morseframes::FieldAnnotation{morseframes::FieldAnnotationEntry{1, 2}},
      morseframes::FieldAnnotation{
          morseframes::FieldAnnotationEntry{0, 1}, morseframes::FieldAnnotationEntry{1, 1}},
  };
  morseframes::FieldAnnotationStore restricted_store(
      restricted_annotations, std::vector<morseframes::SimplexId>{0, 2}, 2, 3);
  assert(field_annotation_equals(restricted_store.annotation(0), {{0, 1}}));
  assert(field_annotation_equals(restricted_store.annotation(2), {{0, 1}, {1, 1}}));

  bool rejected = false;
  try {
    (void)restricted_store.annotation(1);
  } catch (const std::out_of_range&) {
    rejected = true;
  }
  assert(rejected);
}

void test_monotonicity_rejection() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 2.0);
  add_simplex(complex, {1}, 0.0);
  add_simplex(complex, {0, 1}, 1.0);

  bool rejected = false;
  try {
    complex.finalize();
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  assert(rejected);
}

void test_simplex_tree_builder_gudhi_style_insert() {
  morseframes::SimplexTreeBuilder builder;

  const bool inserted_triangle = builder.insert(std::vector<morseframes::VertexId>{0, 1, 2}, 2.0);
  assert(inserted_triangle);
  assert(builder.size() == 7);
  assert(builder.num_simplices() == 7);
  assert(builder.num_vertices() == 3);
  assert(builder.max_dimension() == 2);
  assert(builder.contains({0}));
  assert(builder.contains({0, 1}));
  assert(builder.find_simplex({0, 1, 2}));
  assert(close(builder.simplex_filtration({0, 1, 2}), 2.0));

  const bool lowered_vertex = builder.insert({0}, 0.5);
  assert(lowered_vertex);
  assert(close(builder.filtration({0}), 0.5));
  const bool repeated_vertex_changed = builder.insert({0}, 0.5);
  assert(!repeated_vertex_changed);

  const auto filtration = builder.get_filtration();
  assert(filtration.size() == 7);
  assert((filtration.front().first == std::vector<morseframes::VertexId>{0}));
  assert(close(filtration.front().second, 0.5));

  auto complex = builder.finalize(false);
  assert(builder.size() == 7);
  assert(complex.size() == 7);
  assert(complex.find_simplex({0, 1, 2}) != morseframes::kInvalidSimplex);

  (void)builder.finalize(true);
  assert(builder.size() == 0);
  assert(builder.max_dimension() == -1);
}

void test_simplex_tree_builder_strict_duplicate_rejection() {
  morseframes::SimplexTreeBuilder builder("strict");
  const bool inserted = builder.insert({0}, 1.0);
  assert(inserted);

  bool rejected = false;
  try {
    (void)builder.insert({0}, 2.0);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  assert(rejected);
}

void test_simplex_tree_builder_explicit_insert_can_be_nonclosed() {
  morseframes::SimplexTreeBuilder builder;
  const bool inserted = builder.insert_simplex_only({0, 1}, 1.0);
  assert(inserted);
  assert(builder.size() == 1);
  assert(!builder.contains({0}));

  bool rejected = false;
  try {
    (void)builder.finalize(false);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  assert(rejected);
}

struct FakeGudhiLikeSimplexTree {
  using Handle = std::size_t;

  std::vector<Handle> filtration_simplex_range() const {
    return {0, 1, 2};
  }

  std::vector<morseframes::VertexId> simplex_vertex_range(Handle handle) const {
    return simplices.at(handle).first;
  }

  double filtration(Handle handle) const {
    return simplices.at(handle).second;
  }

  std::vector<std::pair<std::vector<morseframes::VertexId>, double>> simplices = {
      {{1}, 0.0},
      {{0}, 0.0},
      {{1, 0}, 1.0},
  };
};

void test_filtered_complex_from_simplex_tree_adapter() {
  FakeGudhiLikeSimplexTree simplex_tree;
  auto complex = morseframes::filtered_complex_from_simplex_tree(simplex_tree);

  assert(complex.size() == 3);
  assert(complex.find_simplex({0}) != morseframes::kInvalidSimplex);
  assert(complex.find_simplex({1}) != morseframes::kInvalidSimplex);
  const auto edge = complex.find_simplex({0, 1});
  assert(edge != morseframes::kInvalidSimplex);
  assert(close(complex.filtration(edge), 1.0));
}

void test_f_sequence_builder_accepts_simplex_tree_view() {
  FakeGudhiLikeSimplexTree simplex_tree;
  morseframes::SimplexTreeComplexView<FakeGudhiLikeSimplexTree> view(simplex_tree);

  assert(view.size() == 3);
  assert(view.num_levels() == 2);
  assert(view.find_simplex({0, 1}) != morseframes::kInvalidSimplex);

  const auto sequence = FSequenceBuilder(view).build_saturated();
  const auto references =
      morseframes::MorseReferenceComputer(view, sequence).compute_full_references();
  const auto frame = morseframes::MorseReferenceFrameBuilder(view).build_saturated();

  std::size_t regular_pairs = 0;
  for (const auto& step : sequence.steps()) {
    if (step.type == morseframes::MorseStepType::RegularPair) {
      ++regular_pairs;
      assert(view.level(step.sigma) == view.level(step.tau));
    }
  }

  assert(sequence.critical_simplices().size() + 2 * regular_pairs == view.size());
  assert(frame.sequence.critical_simplices() == sequence.critical_simplices());
  assert(frame.references == references);
}

void test_one_vertex() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 0.0);
  auto diagram = run_reference(complex);
  auto finite = morseframes::off_diagonal_pairs(diagram);

  assert(finite.empty());
  assert(diagram.essential.size() == 1);
  assert(diagram.essential.front().dimension == 0);
  assert(close(diagram.essential.front().birth_value, 0.0));
}

void test_reducer_skips_initially_zero_boundaries() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 0.0);
  complex.finalize();

  const auto sequence = FSequenceBuilder(complex).build_saturated();
  const auto references = morseframes::MorseReferenceComputer(complex, sequence).compute_full_references();
  const auto result =
      morseframes::MorseReferencePersistenceReducer(complex, sequence, references).compute_with_metrics();

  assert(result.metrics.boundary_plan_face_scans == 0);
  assert(result.metrics.boundary_annotation_candidate_criticals == 0);
  assert(result.metrics.boundary_annotation_zero_skipped_criticals == 1);
  assert(result.metrics.boundary_annotation_zero_skipped_faces == 0);
  assert(result.metrics.boundary_annotation_xors == 0);
  assert(result.diagram.essential.size() == 1);
}

void test_two_vertices_joined_by_later_edge() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 0.0);
  add_simplex(complex, {1}, 0.0);
  add_simplex(complex, {0, 1}, 1.0);
  auto diagram = run_reference(complex);
  auto finite = morseframes::off_diagonal_pairs(diagram);

  assert(finite.size() == 1);
  assert(count_finite_dim(finite, 0) == 1);
  assert(close(finite.front().birth_value, 0.0));
  assert(close(finite.front().death_value, 1.0));
  assert(count_essential_dim(diagram, 0) == 1);
}

void test_same_level_edge_has_no_off_diagonal_pair() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 0.0);
  add_simplex(complex, {1}, 0.0);
  add_simplex(complex, {0, 1}, 0.0);
  auto diagram = run_reference(complex);
  auto finite = morseframes::off_diagonal_pairs(diagram);

  assert(finite.empty());
  assert(count_essential_dim(diagram, 0) == 1);
}

void test_triangle_boundary() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 0.0);
  add_simplex(complex, {1}, 0.0);
  add_simplex(complex, {2}, 0.0);
  add_simplex(complex, {0, 1}, 1.0);
  add_simplex(complex, {0, 2}, 1.0);
  add_simplex(complex, {1, 2}, 1.0);
  auto diagram = run_reference(complex);
  auto finite = morseframes::off_diagonal_pairs(diagram);

  assert(finite.size() == 2);
  assert(count_finite_dim(finite, 0) == 2);
  assert(count_essential_dim(diagram, 0) == 1);
  assert(count_essential_dim(diagram, 1) == 1);
}

void test_filled_triangle() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 0.0);
  add_simplex(complex, {1}, 0.0);
  add_simplex(complex, {2}, 0.0);
  add_simplex(complex, {0, 1}, 1.0);
  add_simplex(complex, {0, 2}, 1.0);
  add_simplex(complex, {1, 2}, 1.0);
  add_simplex(complex, {0, 1, 2}, 2.0);
  auto diagram = run_reference(complex);
  auto finite = morseframes::off_diagonal_pairs(diagram);

  assert(finite.size() == 3);
  assert(count_finite_dim(finite, 0) == 2);
  assert(count_finite_dim(finite, 1) == 1);
  assert(count_essential_dim(diagram, 0) == 1);
  assert(count_essential_dim(diagram, 1) == 0);
}

void test_standard_prime_field_persistence() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 0.0);
  add_simplex(complex, {1}, 0.0);
  add_simplex(complex, {2}, 0.0);
  add_simplex(complex, {0, 1}, 1.0);
  add_simplex(complex, {0, 2}, 1.0);
  add_simplex(complex, {1, 2}, 1.0);
  add_simplex(complex, {0, 1, 2}, 2.0);
  complex.finalize();

  const auto z2 = morseframes::compute_standard_z2_persistence(complex);
  assert_same_barcode(morseframes::compute_standard_prime_field_persistence(complex, 2), z2);
  assert_same_barcode(morseframes::compute_standard_prime_field_persistence(complex, 3), z2);
  assert_same_barcode(morseframes::compute_standard_prime_field_persistence(complex, 5), z2);
  assert_same_barcode(morseframes::compute_standard_prime_field_persistence(complex, 7), z2);

  bool rejected = false;
  try {
    (void)morseframes::compute_standard_prime_field_persistence(complex, 6);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  assert(rejected);
}

void test_morse_reference_prime_field_persistence() {
  FilteredSimplicialComplex complex;
  const std::vector<double> values = {1.0, 0.0, 1.0, 0.0};
  add_weighted_closure(complex, {0, 1, 2}, values);
  add_weighted_closure(complex, {0, 2, 3}, values);
  complex.finalize();

  auto check_sequence = [&](morseframes::MorseSequence sequence) {
    assert_field_reference_matches_standard(complex, sequence, 3);
    assert_field_reference_matches_standard(complex, sequence, 5);
  };

  check_sequence(FSequenceBuilder(complex).build_saturated());
  check_sequence(FSequenceBuilder(complex).build_plateau_greedy());
  check_sequence(FSequenceBuilder(complex).build_same_level_reduction());
  check_sequence(FSequenceBuilder(complex).build_f_max());
  check_sequence(FSequenceBuilder(complex).build_f_min());
  check_sequence(FSequenceBuilder(complex).build_flooding_max());
  check_sequence(FSequenceBuilder(complex).build_flooding_min());
  check_sequence(FSequenceBuilder(complex).build_flooding_minmax());
  check_sequence(FSequenceBuilder(complex).build_flooding_maxmin());

  const auto coreduction_sequence = FSequenceBuilder(complex).build_same_level_reduction();
  const auto references =
      morseframes::MorseFieldReferenceComputer(complex, coreduction_sequence, 3)
          .compute_full_references();
  bool saw_signed_coefficient = false;
  for (const auto& annotation : references) {
    for (const auto& entry : annotation) {
      if (entry.coefficient == 2) {
        saw_signed_coefficient = true;
      }
    }
  }
  assert(saw_signed_coefficient);

  bool rejected = false;
  try {
    (void)morseframes::compute_morse_reference_prime_field_persistence(
        complex, coreduction_sequence, 6);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  assert(rejected);
}

void test_morse_coreference_prime_field_persistence() {
  FilteredSimplicialComplex complex;
  const std::vector<double> values = {1.0, 0.0, 1.0, 0.0};
  add_weighted_closure(complex, {0, 1, 2}, values);
  add_weighted_closure(complex, {0, 2, 3}, values);
  complex.finalize();

  auto check_sequence = [&](morseframes::MorseSequence sequence) {
    assert_field_coreference_matches_standard(complex, sequence, 3);
    assert_field_coreference_matches_standard(complex, sequence, 5);
  };

  check_sequence(FSequenceBuilder(complex).build_saturated());
  check_sequence(FSequenceBuilder(complex).build_plateau_greedy());
  check_sequence(FSequenceBuilder(complex).build_same_level_reduction());
  check_sequence(FSequenceBuilder(complex).build_f_max());
  check_sequence(FSequenceBuilder(complex).build_f_min());
  check_sequence(FSequenceBuilder(complex).build_flooding_max());
  check_sequence(FSequenceBuilder(complex).build_flooding_min());
  check_sequence(FSequenceBuilder(complex).build_flooding_minmax());
  check_sequence(FSequenceBuilder(complex).build_flooding_maxmin());

  const auto saturated_sequence = FSequenceBuilder(complex).build_saturated();
  const auto coreferences =
      morseframes::MorseFieldCoreferenceComputer(complex, saturated_sequence, 3)
          .compute_full_coreferences();
  bool saw_signed_coefficient = false;
  for (const auto& annotation : coreferences) {
    for (const auto& entry : annotation) {
      if (entry.coefficient == 2) {
        saw_signed_coefficient = true;
      }
    }
  }
  assert(saw_signed_coefficient);

  bool rejected = false;
  try {
    (void)morseframes::compute_morse_coreference_prime_field_persistence(
        complex, saturated_sequence, 6);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  assert(rejected);
}

void test_tetrahedron_boundary() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 0.0);
  add_simplex(complex, {1}, 0.0);
  add_simplex(complex, {2}, 0.0);
  add_simplex(complex, {3}, 0.0);
  add_simplex(complex, {0, 1}, 1.0);
  add_simplex(complex, {0, 2}, 1.0);
  add_simplex(complex, {0, 3}, 1.0);
  add_simplex(complex, {1, 2}, 1.0);
  add_simplex(complex, {1, 3}, 1.0);
  add_simplex(complex, {2, 3}, 1.0);
  add_simplex(complex, {0, 1, 2}, 2.0);
  add_simplex(complex, {0, 1, 3}, 2.0);
  add_simplex(complex, {0, 2, 3}, 2.0);
  add_simplex(complex, {1, 2, 3}, 2.0);
  auto diagram = run_reference(complex);
  auto finite = morseframes::off_diagonal_pairs(diagram);

  assert(finite.size() == 6);
  assert(count_finite_dim(finite, 0) == 3);
  assert(count_finite_dim(finite, 1) == 3);
  assert(count_essential_dim(diagram, 0) == 1);
  assert(count_essential_dim(diagram, 2) == 1);
}

void test_filled_tetrahedron() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 0.0);
  add_simplex(complex, {1}, 0.0);
  add_simplex(complex, {2}, 0.0);
  add_simplex(complex, {3}, 0.0);
  add_simplex(complex, {0, 1}, 1.0);
  add_simplex(complex, {0, 2}, 1.0);
  add_simplex(complex, {0, 3}, 1.0);
  add_simplex(complex, {1, 2}, 1.0);
  add_simplex(complex, {1, 3}, 1.0);
  add_simplex(complex, {2, 3}, 1.0);
  add_simplex(complex, {0, 1, 2}, 2.0);
  add_simplex(complex, {0, 1, 3}, 2.0);
  add_simplex(complex, {0, 2, 3}, 2.0);
  add_simplex(complex, {1, 2, 3}, 2.0);
  add_simplex(complex, {0, 1, 2, 3}, 3.0);
  auto diagram = run_reference(complex);
  auto finite = morseframes::off_diagonal_pairs(diagram);

  assert(finite.size() == 7);
  assert(count_finite_dim(finite, 0) == 3);
  assert(count_finite_dim(finite, 1) == 3);
  assert(count_finite_dim(finite, 2) == 1);
  assert(count_essential_dim(diagram, 0) == 1);
  assert(count_essential_dim(diagram, 1) == 0);
  assert(count_essential_dim(diagram, 2) == 0);
}

void test_same_level_filled_tetrahedron_is_contractible() {
  FilteredSimplicialComplex complex;
  const std::vector<double> values = {0.0, 0.0, 0.0, 0.0};
  add_weighted_closure(complex, {0, 1, 2, 3}, values);
  auto diagram = run_reference(complex);
  auto finite = morseframes::off_diagonal_pairs(diagram);

  assert(finite.empty());
  assert(diagram.essential.size() == 1);
  assert(count_essential_dim(diagram, 0) == 1);
}

void test_lower_star_two_triangle_strip() {
  FilteredSimplicialComplex complex;
  const std::vector<double> values = {0.3, 0.1, 0.8, 0.2};
  add_weighted_closure(complex, {0, 1, 2}, values);
  add_weighted_closure(complex, {1, 2, 3}, values);
  auto diagram = run_reference(complex);

  assert(count_essential_dim(diagram, 0) == 1);
}

void test_lower_star_three_dimensional_pair() {
  FilteredSimplicialComplex complex;
  const std::vector<double> values = {0.25, 0.9, 0.4, 0.1, 0.7};
  add_weighted_closure(complex, {0, 1, 2, 3}, values, 0.05);
  add_weighted_closure(complex, {1, 2, 3, 4}, values, 0.05);
  auto diagram = run_reference(complex);

  assert(count_essential_dim(diagram, 0) == 1);
}

void test_instrumentation_metrics() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 0.0);
  add_simplex(complex, {1}, 0.0);
  add_simplex(complex, {2}, 0.0);
  add_simplex(complex, {0, 1}, 1.0);
  add_simplex(complex, {0, 2}, 1.0);
  add_simplex(complex, {1, 2}, 1.0);
  add_simplex(complex, {0, 1, 2}, 2.0);
  complex.finalize();

  const auto result = morseframes::run_instrumented_persistence(complex);
  const auto& metrics = result.metrics;

  assert(metrics.structural.num_simplices == 7);
  assert(metrics.structural.num_critical + 2 * metrics.structural.num_regular_pairs ==
         metrics.structural.num_simplices);
  assert(metrics.structural.w_boundary_plus_size <= metrics.structural.num_simplices);
  assert(metrics.structural.w_coboundary_plus_size <= metrics.structural.num_simplices);
  assert(metrics.reference_annotations.count == metrics.structural.w_boundary_plus_size);
  assert(metrics.coreference_annotations.count == metrics.structural.w_coboundary_plus_size);
  assert(metrics.reference_annotations.average_size >= 0.0);
  assert(metrics.coreference_annotations.average_size >= 0.0);
  assert(metrics.timings.f_sequence_ms >= 0.0);
  assert(metrics.timings.reference_compute_ms >= 0.0);
  assert(metrics.timings.reference_reduce_ms >= 0.0);
  assert(metrics.timings.coreference_compute_ms >= 0.0);
  assert(metrics.timings.coreference_reduce_ms >= 0.0);
  assert(metrics.timings.standard_reduce_ms >= 0.0);

  assert_same_barcode(result.reference_diagram, result.standard_diagram);
  assert_same_barcode(result.coreference_diagram, result.standard_diagram);
}

}  // namespace

int main() {
  test_boundary_and_coboundary();
  test_inverse_annotation_store();
  test_field_annotation_store();
  test_monotonicity_rejection();
  test_simplex_tree_builder_gudhi_style_insert();
  test_simplex_tree_builder_strict_duplicate_rejection();
  test_simplex_tree_builder_explicit_insert_can_be_nonclosed();
  test_filtered_complex_from_simplex_tree_adapter();
  test_f_sequence_builder_accepts_simplex_tree_view();
  test_one_vertex();
  test_reducer_skips_initially_zero_boundaries();
  test_two_vertices_joined_by_later_edge();
  test_same_level_edge_has_no_off_diagonal_pair();
  test_triangle_boundary();
  test_filled_triangle();
  test_standard_prime_field_persistence();
  test_morse_reference_prime_field_persistence();
  test_morse_coreference_prime_field_persistence();
  test_tetrahedron_boundary();
  test_filled_tetrahedron();
  test_same_level_filled_tetrahedron_is_contractible();
  test_lower_star_two_triangle_strip();
  test_lower_star_three_dimensional_pair();
  test_instrumentation_metrics();

  std::cout << "All Morse persistence prototype tests passed.\n";
  return 0;
}
