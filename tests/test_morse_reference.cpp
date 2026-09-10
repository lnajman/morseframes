#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <future>
#include <initializer_list>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "morseframes/coreference_persistence.hpp"
#include "morseframes/debug_checks.hpp"
#include "morseframes/field_annotation_store.hpp"
#include "morseframes/filtered_complex.hpp"
#include "morseframes/instrumentation.hpp"
#include "morseframes/inverse_annotation_store.hpp"
#include "morseframes/lower_star_complex.hpp"
#include "morseframes/morse_reference_api.hpp"
#include "morseframes/morse_sequence.hpp"
#include "morseframes/reduction_kernel_sequence.hpp"
#include "morseframes/reference_persistence.hpp"
#include "morseframes/simplex_tree_builder.hpp"
#include "morseframes/standard_persistence.hpp"

namespace {

using morseframes::FilteredSimplicialComplex;
using morseframes::FSequenceBuilder;
using morseframes::PersistenceDiagram;

constexpr double kEps = 1e-12;

void test_bounded_task_executor() {
  morseframes::BoundedTaskExecutor executor(3);
  assert(executor.worker_count() == 3);

  std::atomic<int> active{0};
  std::atomic<int> peak{0};
  std::vector<std::future<int>> futures;
  for (int value = 0; value < 18; ++value) {
    futures.push_back(executor.submit([value, &active, &peak]() {
      const int current = active.fetch_add(1) + 1;
      int observed = peak.load();
      while (observed < current &&
             !peak.compare_exchange_weak(observed, current)) {
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
      active.fetch_sub(1);
      return value * value;
    }));
  }
  for (int value = 0; value < 18; ++value) {
    assert(executor.get(futures[value]) == value * value);
  }
  assert(peak.load() > 1);
  assert(peak.load() <= 3);

  auto nested = executor.submit([&executor]() {
    std::vector<std::future<int>> inner;
    for (int value = 1; value <= 8; ++value) {
      inner.push_back(executor.submit([value]() { return value; }));
    }
    int total = 0;
    for (auto& future : inner) {
      total += executor.get(future);
    }
    return total;
  });
  assert(executor.get(nested) == 36);

  auto failure = executor.submit([]() -> int {
    throw std::runtime_error("executor failure");
  });
  bool propagated = false;
  try {
    (void)executor.get(failure);
  } catch (const std::runtime_error&) {
    propagated = true;
  }
  assert(propagated);

  morseframes::BoundedTaskExecutor single_worker_executor(1);
  auto single_worker_nested = single_worker_executor.submit(
      [&single_worker_executor]() {
        auto inner = single_worker_executor.submit([]() { return 17; });
        return single_worker_executor.get(inner);
      });
  assert(single_worker_executor.get(single_worker_nested) == 17);
}

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
  auto morse_metrics_result =
      morseframes::MorseReferencePersistenceReducer(complex, sequence, references)
          .compute_with_metrics();
  auto public_morse_diagram =
      morseframes::compute_morse_reference_persistence(
          complex, morseframes::MorseSequenceStrategy::Saturated);
  assert_same_barcode(morse_diagram, morse_metrics_result.diagram);
  assert_same_barcode(morse_diagram, public_morse_diagram);

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
  auto deferred_working_set = deferred_input.reduction_plan.working_set;
  std::sort(deferred_working_set.begin(), deferred_working_set.end());
  assert(deferred_working_set == expected_working_set);
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
  check_flooding_sequence(
      FSequenceBuilder(complex).build_flooding_reduction_kernel());
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

template <class T, class = void>
struct HasFMaxBuilderMethod : std::false_type {};

template <class T>
struct HasFMaxBuilderMethod<T, std::void_t<decltype(std::declval<T>().build_f_max())>>
    : std::true_type {};

void test_reduction_kernel_lightweight_initialization() {
  using Lean = morseframes::ReductionKernelSequenceBuilder<>;
  static_assert(!HasFMaxBuilderMethod<Lean>::value);
  static_assert(HasFMaxBuilderMethod<FSequenceBuilder<>>::value);
  static_assert(!std::is_convertible_v<Lean*, FSequenceBuilder<>*>);

  // The low-level view/builder can be empty, although the owning complex's
  // public finalize() API requires at least one simplex.
  const FilteredSimplicialComplex empty;
  const Lean empty_builder(empty);
  assert(empty_builder.build_flooding_reduction_kernel().steps().empty());
  assert(empty_builder.build_flooding_reduction_kernel_parallel(2).steps().empty());

  struct CountingView : FilteredSimplicialComplex {
    mutable std::size_t level_reads = 0, dimension_reads = 0;
    bool override_order = false;
    std::vector<morseframes::SimplexId> order;
    const std::vector<morseframes::SimplexId>& filtration_order() const {
      return override_order ? order : FilteredSimplicialComplex::filtration_order();
    }
    morseframes::LevelId level(morseframes::SimplexId id) const {
      ++level_reads;
      return FilteredSimplicialComplex::level(id);
    }
    std::uint16_t dimension(morseframes::SimplexId id) const {
      ++dimension_reads;
      return FilteredSimplicialComplex::dimension(id);
    }
  };
  CountingView view;
  add_simplex(view, {0}, 0.0);
  add_simplex(view, {1}, 0.0);
  add_simplex(view, {0, 1}, 0.0);
  view.finalize();
  view.level_reads = view.dimension_reads = 0;
  FSequenceBuilder eager(view);
  assert(view.level_reads == view.size() && view.dimension_reads == view.size());
  view.level_reads = view.dimension_reads = 0;
  morseframes::ReductionKernelSequenceBuilder lean(view);
  assert(view.level_reads == 0 && view.dimension_reads == 0);

  // Both constructors reject the same malformed permutations before kernels
  // run, including the sentinel ID and duplicate entries that omit a simplex.
  view.override_order = true;
  for (const auto& order : std::vector<std::vector<morseframes::SimplexId>>{
           {0, 1}, {0, 1, 3}, {0, 1, morseframes::kInvalidSimplex}, {0, 0, 2}}) {
    view.order = order;
    std::string eager_error, lean_error;
    try { (void)FSequenceBuilder(view); }
    catch (const std::logic_error& e) { eager_error = e.what(); }
    try { (void)morseframes::ReductionKernelSequenceBuilder(view); }
    catch (const std::logic_error& e) { lean_error = e.what(); }
    assert(!eager_error.empty() && eager_error == lean_error);
  }

  // The lightweight path also accepts non-owning/generic complex views.
  FakeGudhiLikeSimplexTree tree;
  morseframes::SimplexTreeComplexView<FakeGudhiLikeSimplexTree> tree_view(tree);
  const auto sequence = morseframes::ReductionKernelSequenceBuilder(tree_view)
                            .build_flooding_reduction_kernel();
  const auto expected = FSequenceBuilder(tree_view).build_flooding_reduction_kernel();
  assert(sequence.steps().size() == expected.steps().size());
  for (std::size_t i = 0; i < sequence.steps().size(); ++i) {
    const auto& a = sequence.steps()[i];
    const auto& b = expected.steps()[i];
    assert(a.type == b.type && a.sigma == b.sigma && a.tau == b.tau && a.level == b.level);
  }
}

void test_reduction_kernel_lightweight_persistence() {
  std::vector<FilteredSimplicialComplex> inputs;
  for (std::size_t vertices : {3, 7, 8}) {
    for (bool multilevel : {false, true}) {
      FilteredSimplicialComplex complex;
      std::vector<double> values(vertices, 0.0);
      std::vector<morseframes::VertexId> facet;
      for (std::size_t v = 0; v < vertices; ++v) {
        facet.push_back(static_cast<morseframes::VertexId>(v));
        if (multilevel) values[v] = static_cast<double>(v % 3);
      }
      add_weighted_closure(complex, facet, values);
      inputs.push_back(std::move(complex));
    }
  }
  FilteredSimplicialComplex graph;
  const std::vector<double> graph_values(131, 0.0);
  for (morseframes::VertexId v = 1; v < graph_values.size(); ++v) {
    add_weighted_closure(graph, {0, v}, graph_values);
  }
  inputs.push_back(std::move(graph));

  for (auto original : inputs) {
    original.finalize();
    for (bool cached : {false, true}) {
      auto complex = original;
      if (cached) complex.prepare_same_level_closure_cache();
      const FSequenceBuilder eager(complex);
      const auto expected = eager.build_flooding_reduction_kernel();
      const auto f_max = eager.build_f_max();
      const auto compare = [&](const auto& a, const auto& b) {
        assert(a.steps().size() == b.steps().size());
        for (std::size_t i = 0; i < a.steps().size(); ++i) {
          const auto& x = a.steps()[i];
          const auto& y = b.steps()[i];
          assert(x.type == y.type && x.sigma == y.sigma && x.tau == y.tau &&
                 x.level == y.level);
        }
        morseframes::validate_morse_sequence(complex, b);
      };
      const morseframes::ReductionKernelSequenceBuilder lean(complex);
      const auto copied = lean;
      compare(expected, copied.build_flooding_reduction_kernel());
      // No mutable lazy cache: repeated const calls can run independently.
      auto concurrent = std::async(std::launch::async, [&]() {
        return lean.build_flooding_reduction_kernel();
      });
      compare(expected, lean.build_flooding_reduction_kernel());
      compare(expected, concurrent.get());
      for (std::size_t workers : {1, 2, 4, 8}) {
        std::size_t callbacks = 0;
        compare(expected, lean.build_flooding_reduction_kernel_parallel_with_step_callback(
            [&](const auto& sequence, const auto& step) {
              assert(sequence.steps().size() == ++callbacks);
              const auto& e = expected.steps()[callbacks - 1];
              assert(e.type == step.type && e.sigma == step.sigma &&
                     e.tau == step.tau && e.level == step.level);
            }, workers));
        assert(callbacks == expected.steps().size());
        for (bool detailed : {false, true}) {
          morseframes::MorseSequenceBuildMetrics metrics;
          compare(expected, morseframes::ReductionKernelSequenceBuilder(complex, &metrics, detailed)
                                .build_flooding_reduction_kernel_parallel(workers));
        }
      }
      if (!expected.steps().empty()) {
        bool propagated = false;
        try {
          lean.build_flooding_reduction_kernel_with_step_callback(
              [](const auto&, const auto&) { throw std::runtime_error("callback failure"); });
        } catch (const std::runtime_error& e) {
          propagated = std::string(e.what()) == "callback failure";
        }
        assert(propagated);
        compare(expected, lean.build_flooding_reduction_kernel());
      }

      const auto references = morseframes::MorseReferenceComputer(complex, expected)
                                  .compute_full_references();
      const auto standard = morseframes::compute_standard_z2_persistence(complex);
      const morseframes::MorseReferenceFrameBuilder frame_builder(complex);
      for (bool parallel : {false, true}) {
        const auto frame = parallel ? frame_builder.build_flooding_reduction_kernel_parallel(4)
                                    : frame_builder.build_flooding_reduction_kernel();
        compare(expected, frame.sequence);
        assert(frame.references == references);
        auto compact = parallel ? frame_builder.build_flooding_reduction_kernel_parallel_reduction_input(4)
                                : frame_builder.build_flooding_reduction_kernel_reduction_input();
        compare(expected, compact.sequence);
        auto reducer = morseframes::MorseReferencePersistenceReducer(
            complex, compact.sequence, std::move(compact.reduction_plan),
            std::move(compact.annotations));
        assert_same_barcode(standard, reducer.compute());
        const auto strategy = parallel ? morseframes::MorseSequenceStrategy::FloodingReductionKernelParallel
                                       : morseframes::MorseSequenceStrategy::FloodingReductionKernel;
        assert_same_barcode(standard, morseframes::compute_morse_reference_persistence(complex, strategy));
      }
      assert_field_reference_matches_standard(complex, expected, 3);
      assert_field_coreference_matches_standard(complex, expected, 3);
      compare(f_max, eager.build_f_max());
    }
  }
}

void test_process_lower_stars_triangle_boundary() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 2.0);
  add_simplex(complex, {1}, 1.0);
  add_simplex(complex, {2}, 0.0);
  add_simplex(complex, {0, 1}, 2.0);
  add_simplex(complex, {0, 2}, 2.0);
  add_simplex(complex, {1, 2}, 1.0);
  complex.finalize();

  const auto sequence = FSequenceBuilder(complex).build_process_lower_stars();
  morseframes::validate_morse_sequence(complex, sequence);
  assert(sequence.steps().size() == 4);

  const auto v0 = complex.find_simplex({0});
  const auto v1 = complex.find_simplex({1});
  const auto v2 = complex.find_simplex({2});
  const auto e01 = complex.find_simplex({0, 1});
  const auto e02 = complex.find_simplex({0, 2});
  const auto e12 = complex.find_simplex({1, 2});
  const auto& steps = sequence.steps();
  assert(steps[0].type == morseframes::MorseStepType::Critical &&
         steps[0].sigma == v2);
  assert(steps[1].type == morseframes::MorseStepType::RegularPair &&
         steps[1].sigma == v1 && steps[1].tau == e12);
  assert(steps[2].type == morseframes::MorseStepType::RegularPair &&
         steps[2].sigma == v0 && steps[2].tau == e02);
  assert(steps[3].type == morseframes::MorseStepType::Critical &&
         steps[3].sigma == e01);

  morseframes::MorseSequenceBuildMetrics parallel_metrics;
  const auto parallel_sequence =
      FSequenceBuilder(complex, &parallel_metrics)
          .build_process_lower_stars_parallel(2);
  morseframes::validate_morse_sequence(complex, parallel_sequence);
  assert(sequence.steps().size() == parallel_sequence.steps().size());
  for (std::size_t index = 0; index < sequence.steps().size(); ++index) {
    const auto& expected = sequence.steps()[index];
    const auto& actual = parallel_sequence.steps()[index];
    assert(expected.type == actual.type);
    assert(expected.sigma == actual.sigma);
    assert(expected.tau == actual.tau);
    assert(expected.level == actual.level);
  }
  assert(parallel_metrics.process_lower_stars_executor_workers == 2);
  assert(parallel_metrics.process_lower_stars_parallel_tasks == 2);
  assert(parallel_metrics.process_lower_stars_count == 3);
  assert(parallel_metrics.process_lower_stars_max_star_size == 3);
  assert(parallel_metrics.process_lower_stars_min_task_load == 3);
  assert(parallel_metrics.process_lower_stars_max_task_load == 3);
  assert(parallel_metrics.process_lower_stars_setup_nanoseconds > 0);
  assert(parallel_metrics.process_lower_stars_local_wall_nanoseconds > 0);
  assert(parallel_metrics.process_lower_stars_replay_nanoseconds > 0);
  assert(parallel_metrics.process_lower_stars_cleanup_nanoseconds > 0);
  assert(parallel_metrics.process_lower_stars_cleanup_nanoseconds ==
         parallel_metrics.process_lower_stars_events_index_cleanup_nanoseconds +
         parallel_metrics.process_lower_stars_keys_cleanup_nanoseconds +
         parallel_metrics.process_lower_stars_membership_cleanup_nanoseconds +
         parallel_metrics.process_lower_stars_executor_cleanup_nanoseconds +
         parallel_metrics.process_lower_stars_vertices_cleanup_nanoseconds);
  assert(parallel_metrics.process_lower_stars_setup_nanoseconds >=
         parallel_metrics.process_lower_stars_output_init_nanoseconds +
         parallel_metrics.process_lower_stars_vertex_order_nanoseconds +
         parallel_metrics.process_lower_stars_executor_init_nanoseconds +
         parallel_metrics.process_lower_stars_storage_init_nanoseconds +
         parallel_metrics.process_lower_stars_owner_keys_nanoseconds +
         parallel_metrics.process_lower_stars_partition_nanoseconds);
  assert(parallel_metrics.process_lower_stars_local_wall_nanoseconds >=
         parallel_metrics.process_lower_stars_schedule_nanoseconds +
         parallel_metrics.process_lower_stars_execution_nanoseconds);
  assert(parallel_metrics.process_lower_stars_cumulative_task_nanoseconds > 0);
  assert(parallel_metrics.process_lower_stars_max_task_nanoseconds >=
         parallel_metrics.process_lower_stars_min_task_nanoseconds);
  const auto single_worker_sequence =
      FSequenceBuilder(complex).build_process_lower_stars_parallel(1);
  assert(single_worker_sequence.steps().size() == sequence.steps().size());
  for (std::size_t index = 0; index < sequence.steps().size(); ++index) {
    assert(single_worker_sequence.steps()[index].type == steps[index].type);
    assert(single_worker_sequence.steps()[index].sigma == steps[index].sigma);
    assert(single_worker_sequence.steps()[index].tau == steps[index].tau);
  }

  const auto diagram = morseframes::compute_morse_reference_persistence(
      complex, morseframes::MorseSequenceStrategy::ProcessLowerStars);
  assert_same_barcode(diagram,
                      morseframes::compute_standard_z2_persistence(complex));
  assert(morseframes::morse_sequence_strategy_from_name("process-lower-stars") ==
         morseframes::MorseSequenceStrategy::ProcessLowerStars);
  assert(morseframes::morse_sequence_strategy_from_name(
             "process-lower-stars-parallel") ==
         morseframes::MorseSequenceStrategy::ProcessLowerStarsParallel);

  FilteredSimplicialComplex tied_vertices;
  add_simplex(tied_vertices, {0}, 0.0);
  add_simplex(tied_vertices, {1}, 0.0);
  add_simplex(tied_vertices, {0, 1}, 0.0);
  tied_vertices.finalize();
  bool rejected_tie = false;
  try {
    (void)FSequenceBuilder(tied_vertices).build_process_lower_stars();
  } catch (const std::invalid_argument&) {
    rejected_tie = true;
  }
  assert(rejected_tie);

  FilteredSimplicialComplex delayed_edge;
  add_simplex(delayed_edge, {0}, 0.0);
  add_simplex(delayed_edge, {1}, 1.0);
  add_simplex(delayed_edge, {0, 1}, 2.0);
  delayed_edge.finalize();
  bool rejected_extension = false;
  try {
    (void)FSequenceBuilder(delayed_edge).build_process_lower_stars();
  } catch (const std::invalid_argument&) {
    rejected_extension = true;
  }
  assert(rejected_extension);
}

// Deliberately slow oracle: rescan the remaining local boundary on every step.
// It shares neither the production counters/XORs nor its heaps/direct-index map.
morseframes::MorseSequence scan_process_lower_stars(
    const FilteredSimplicialComplex& complex) {
  using morseframes::SimplexId;
  std::map<morseframes::VertexId, std::size_t> ranks;
  for (auto id : complex.filtration_order()) {
    if (complex.dimension(id) == 0) {
      const auto rank = ranks.size();
      ranks.emplace(complex.vertices(id)[0], rank);
    }
  }
  std::vector<std::vector<std::size_t>> keys(complex.size());
  std::vector<std::vector<SimplexId>> stars(ranks.size());
  for (SimplexId id = 0; id < complex.size(); ++id) {
    for (auto vertex : complex.vertices(id)) keys[id].push_back(ranks.at(vertex));
    std::sort(keys[id].begin(), keys[id].end(), std::greater<std::size_t>());
    stars[keys[id][0]].push_back(id);
  }
  const auto before = [&](SimplexId a, SimplexId b) {
    return keys[a] != keys[b] ? keys[a] < keys[b] : a < b;
  };
  morseframes::MorseSequence result(complex.size());
  std::vector<bool> classified(complex.size(), false);
  for (const auto& star : stars) {
    auto ordered = star;
    std::sort(ordered.begin(), ordered.end(), before);
    std::size_t remaining = star.size();
    while (remaining) {
      SimplexId pair = morseframes::kInvalidSimplex, face = pair, critical = pair;
      for (auto id : ordered) {
        if (classified[id]) continue;
        std::vector<SimplexId> boundary;
        for (auto f : complex.boundary(id)) {
          if (!classified[f] && keys[f][0] == keys[id][0]) boundary.push_back(f);
        }
        if (boundary.size() == 1) { pair = id; face = boundary[0]; break; }
        if (boundary.empty() && critical == morseframes::kInvalidSimplex) critical = id;
      }
      if (pair != morseframes::kInvalidSimplex) {
        result.add_regular_pair(face, pair, complex.level(pair));
        classified[face] = classified[pair] = true;
        remaining -= 2;
      } else {
        assert(critical != morseframes::kInvalidSimplex);
        result.add_critical(critical, complex.level(critical));
        classified[critical] = true;
        --remaining;
      }
    }
  }
  return result;
}

void test_process_lower_stars_workspace_and_dimensions() {
  const auto compare_step = [](const auto& a, const auto& b) {
    assert(a.type == b.type && a.sigma == b.sigma && a.tau == b.tau && a.level == b.level);
  };
  for (unsigned dimension = 0; dimension <= 7; ++dimension) {
    for (unsigned seed : {0u, 2u, 7u}) {
      // Sparse vertex identifiers and filtration ranks unrelated to their order.
      const unsigned vertices = dimension + 4;
      std::vector<double> values(vertices);
      std::iota(values.begin(), values.end(), -3.0);
      std::mt19937 rng(seed);
      std::shuffle(values.begin(), values.end(), rng);
      FilteredSimplicialComplex complex;
      const auto add_cell = [&](const std::vector<unsigned>& cell) {
        for (unsigned mask = 1; mask < (1u << cell.size()); ++mask) {
          std::vector<morseframes::VertexId> face;
          double value = -std::numeric_limits<double>::infinity();
          for (unsigned i = 0; i < cell.size(); ++i) {
            if (mask & (1u << i)) {
              face.push_back(10 + 17 * cell[i]);
              value = std::max(value, values[cell[i]]);
            }
          }
          complex.add_simplex(face, value);
        }
      };
      for (unsigned i = 0; i < vertices; ++i) add_cell({i});
      std::vector<unsigned> cell(dimension + 1);
      std::iota(cell.begin(), cell.end(), 0);
      add_cell(cell);
      for (auto& v : cell) ++v;
      add_cell(cell);
      if (dimension) add_cell({0, vertices - 1}); // Non-pure when dimension > 1.
      complex.finalize();
      const auto expected = scan_process_lower_stars(complex);
      morseframes::validate_morse_sequence(complex, expected);
      FSequenceBuilder builder(complex);
      for (unsigned workers : {1u, 2u, 4u, 8u}) {
        std::size_t callbacks = 0;
        const auto callback = [&](const auto& prefix, const auto& step) {
          compare_step(expected.steps()[callbacks], step);
          ++callbacks;
          assert(prefix.steps().size() == callbacks);
        };
        const auto actual = builder.build_process_lower_stars_parallel_with_step_callback(callback, workers);
        assert(callbacks == expected.steps().size());
        assert(actual.steps().size() == expected.steps().size());
        morseframes::validate_morse_sequence(complex, actual);
        const auto repeated = builder.build_process_lower_stars();
        assert(repeated.steps().size() == expected.steps().size());
        for (std::size_t i = 0; i < expected.steps().size(); ++i) {
          compare_step(expected.steps()[i], actual.steps()[i]);
          compare_step(expected.steps()[i], repeated.steps()[i]);
        }
      }
    }
  }
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

void test_non_flooding_f_sequence_is_reduced_in_flooding_order() {
  FilteredSimplicialComplex complex;
  add_simplex(complex, {0}, 0.0);
  add_simplex(complex, {1}, 1.0);
  add_simplex(complex, {0, 1}, 1.0);
  complex.finalize();

  const auto v0 = complex.find_simplex({0});
  const auto v1 = complex.find_simplex({1});
  const auto edge = complex.find_simplex({0, 1});

  morseframes::MorseSequence sequence(complex.size());
  sequence.add_critical(v1, complex.level(v1));
  sequence.add_critical(v0, complex.level(v0));
  sequence.add_critical(edge, complex.level(edge));
  morseframes::validate_morse_sequence(complex, sequence);

  const auto standard = morseframes::compute_standard_z2_persistence(complex);
  auto references =
      morseframes::MorseReferenceComputer(complex, sequence).compute_full_references();
  morseframes::validate_reference_invariants(complex, sequence, references);
  const auto reference =
      morseframes::MorseReferencePersistenceReducer(complex, sequence, references).compute();
  const auto reference_with_metrics =
      morseframes::MorseReferencePersistenceReducer(complex, sequence, references)
          .compute_with_metrics()
          .diagram;

  auto coreferences =
      morseframes::MorseCoreferenceComputer(complex, sequence).compute_full_coreferences();
  morseframes::validate_coreference_invariants(complex, sequence, coreferences);
  const auto coreference =
      morseframes::MorseCoreferencePersistenceReducer(complex, sequence, coreferences).compute();

  assert_same_barcode(reference, standard);
  assert_same_barcode(reference_with_metrics, standard);
  assert_same_barcode(coreference, standard);
  assert(morseframes::off_diagonal_pairs(reference).empty());
  assert(reference.essential.size() == 1);
  assert(reference.essential.front().dimension == 0);
  assert(close(reference.essential.front().birth_value, 0.0));

  assert_field_reference_matches_standard(complex, sequence, 3);
  assert_field_coreference_matches_standard(complex, sequence, 3);
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
  check_sequence(FSequenceBuilder(complex).build_flooding_reduction_kernel());
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
  check_sequence(FSequenceBuilder(complex).build_flooding_reduction_kernel());
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

void test_flooding_reduction_kernel_on_shared_facets() {
  const auto assert_same_sequence = [](const auto& expected_sequence,
                                       const auto& actual_sequence) {
    assert(expected_sequence.steps().size() == actual_sequence.steps().size());
    for (std::size_t index = 0; index < expected_sequence.steps().size();
         ++index) {
      const auto& expected = expected_sequence.steps()[index];
      const auto& actual = actual_sequence.steps()[index];
      assert(expected.type == actual.type);
      assert(expected.sigma == actual.sigma);
      assert(expected.tau == actual.tau);
      assert(expected.level == actual.level);
    }
  };

  FilteredSimplicialComplex complex;
  const std::vector<double> values = {0.0, 0.0, 0.0, 0.0};
  add_weighted_closure(complex, {0, 1, 2}, values);
  add_weighted_closure(complex, {1, 2, 3}, values);
  complex.finalize();

  const auto sequence =
      FSequenceBuilder(complex).build_flooding_reduction_kernel();
  morseframes::validate_morse_sequence(complex, sequence);
  assert(sequence.critical_simplices().size() == 1);

  morseframes::MorseSequenceBuildMetrics parallel_metrics;
  const auto parallel_sequence =
      FSequenceBuilder(complex, &parallel_metrics)
          .build_flooding_reduction_kernel_parallel_with_step_callback(
              [](const morseframes::MorseSequence&,
                 const morseframes::MorseStep&) {},
              2);
  morseframes::validate_morse_sequence(complex, parallel_sequence);
  assert_same_sequence(sequence, parallel_sequence);
  // A pair of tiny packed facets no longer warrants task dispatch.
  assert(parallel_metrics.reduction_kernel_parallel_batches == 0);
  assert(parallel_metrics.reduction_kernel_max_parallel_facets == 0);
  assert(parallel_metrics.reduction_kernel_facet_parallel_tasks == 0);
  assert(parallel_metrics.reduction_kernel_executor_workers == 2);
  assert(parallel_metrics.reduction_kernel_facet_discovery_parallel_tasks == 0);
  assert(parallel_metrics.reduction_kernel_facet_discovery_mask_tests > 0);
  assert(parallel_metrics.reduction_kernel_local_coboundary_mask_tests > 0);
  assert(parallel_metrics.reduction_kernel_local_coboundary_visits == 0);
  // Small packed levels compute protected cores with word operations.
  assert(parallel_metrics.reduction_kernel_essential_parallel_tasks == 0);
  assert(parallel_metrics.reduction_kernel_aggregation_rounds > 0);
  const auto single_worker_sequence =
      FSequenceBuilder(complex).build_flooding_reduction_kernel_parallel(1);
  morseframes::validate_morse_sequence(complex, single_worker_sequence);
  assert_same_sequence(sequence, single_worker_sequence);

  FilteredSimplicialComplex multilevel_complex;
  const std::vector<double> multilevel_values = {0.0, 1.0, 2.0, 3.0};
  add_weighted_closure(multilevel_complex, {0, 1, 2}, multilevel_values);
  add_weighted_closure(multilevel_complex, {1, 2, 3}, multilevel_values);
  multilevel_complex.finalize();
  const auto multilevel_sequence =
      FSequenceBuilder(multilevel_complex).build_flooding_reduction_kernel();
  morseframes::MorseSequenceBuildMetrics multilevel_parallel_metrics;
  const auto multilevel_parallel_sequence =
      FSequenceBuilder(multilevel_complex, &multilevel_parallel_metrics)
          .build_flooding_reduction_kernel_parallel_with_step_callback(
              [](const morseframes::MorseSequence&,
                 const morseframes::MorseStep&) {},
              4);
  morseframes::validate_morse_sequence(multilevel_complex,
                                       multilevel_parallel_sequence);
  assert_same_sequence(multilevel_sequence, multilevel_parallel_sequence);
  assert(multilevel_parallel_metrics
             .reduction_kernel_parallel_level_batches > 0);
  assert(multilevel_parallel_metrics
             .reduction_kernel_facet_discovery_parallel_tasks == 0);
  assert(multilevel_parallel_metrics
             .reduction_kernel_essential_parallel_tasks == 0);
  assert(multilevel_parallel_metrics.reduction_kernel_parallel_batches == 0);
  assert(multilevel_parallel_metrics.reduction_kernel_facet_parallel_tasks == 0);
  assert(multilevel_parallel_metrics
             .reduction_kernel_aggregation_parallel_tasks == 0);
  assert(multilevel_parallel_metrics.reduction_kernel_max_parallel_levels ==
         std::min<std::size_t>(multilevel_complex.num_levels(), 4));
  assert(multilevel_parallel_metrics.reduction_kernel_executor_workers == 4);
  const auto automatic_parallel_sequence =
      FSequenceBuilder(multilevel_complex)
          .build_flooding_reduction_kernel_parallel();
  morseframes::validate_morse_sequence(multilevel_complex,
                                       automatic_parallel_sequence);
  assert_same_sequence(multilevel_sequence, automatic_parallel_sequence);

  FilteredSimplicialComplex four_facet_complex;
  for (std::uint32_t edge = 0; edge < 4; ++edge) {
    const std::uint32_t first = 2 * edge;
    add_simplex(four_facet_complex, {first}, 0.0);
    add_simplex(four_facet_complex, {first + 1}, 0.0);
    add_simplex(four_facet_complex, {first, first + 1}, 0.0);
  }
  four_facet_complex.finalize();
  const auto four_facet_sequence =
      FSequenceBuilder(four_facet_complex).build_flooding_reduction_kernel();
  morseframes::MorseSequenceBuildMetrics four_facet_parallel_metrics;
  const auto four_facet_parallel_sequence =
      FSequenceBuilder(four_facet_complex, &four_facet_parallel_metrics)
          .build_flooding_reduction_kernel_parallel_with_step_callback(
              [](const morseframes::MorseSequence&,
                 const morseframes::MorseStep&) {},
              4);
  morseframes::validate_morse_sequence(four_facet_complex,
                                       four_facet_parallel_sequence);
  assert_same_sequence(four_facet_sequence, four_facet_parallel_sequence);
  assert(four_facet_parallel_metrics.reduction_kernel_aggregation_rounds > 0);
  assert(four_facet_parallel_metrics
             .reduction_kernel_aggregation_parallel_tasks == 0);

  const auto strategy =
      morseframes::morse_sequence_strategy_from_name("reduction-kernel");
  assert(strategy ==
         morseframes::MorseSequenceStrategy::FloodingReductionKernel);
  assert(std::string(morseframes::morse_sequence_strategy_name(strategy)) ==
         "flooding-reduction-kernel");
  const auto parallel_strategy = morseframes::morse_sequence_strategy_from_name(
      "reduction-kernel-parallel");
  assert(parallel_strategy ==
         morseframes::MorseSequenceStrategy::FloodingReductionKernelParallel);

  const auto result = morseframes::compute_morse_reference_persistence(
      complex, strategy);
  assert_same_barcode(result,
                      morseframes::compute_standard_z2_persistence(complex));

  auto input = morseframes::MorseReferenceFrameBuilder(complex, true)
                   .build_flooding_reduction_kernel_reduction_input();
  const auto& metrics = input.frame_metrics;
  assert(metrics.sequence_reduction_kernel_levels == 1);
  assert(metrics.sequence_reduction_kernel_rounds > 0);
  assert(metrics.sequence_reduction_kernel_facet_kernels > 0);
  assert(metrics.sequence_reduction_kernel_reductions ==
         metrics.sequence_regular_pairs);
  assert(metrics.sequence_reduction_kernel_perforations ==
         metrics.sequence_criticals);
}

void test_reduction_kernel_packed_core_matches_sparse_cache() {
  const auto check = [](FilteredSimplicialComplex complex) {
    complex.finalize();
    auto cached = complex;
    cached.prepare_same_level_closure_cache();
    // The cache retains the independent sparse incidence/cell implementation.
    const auto expected =
        FSequenceBuilder(cached).build_flooding_reduction_kernel();
    const auto compare = [&](const auto& actual) {
      morseframes::validate_morse_sequence(complex, actual);
      assert(expected.steps().size() == actual.steps().size());
      for (std::size_t i = 0; i < expected.steps().size(); ++i) {
        const auto& a = expected.steps()[i];
        const auto& b = actual.steps()[i];
        assert(a.type == b.type && a.sigma == b.sigma && a.tau == b.tau &&
               a.level == b.level);
      }
    };
    compare(FSequenceBuilder(complex).build_flooding_reduction_kernel());
    compare(FSequenceBuilder(complex).build_flooding_reduction_kernel_parallel(4));
    morseframes::MorseSequenceBuildMetrics metrics;
    compare(FSequenceBuilder(complex, &metrics)
                .build_flooding_reduction_kernel_parallel(4));
    for (std::size_t workers : {1, 4}) {
      morseframes::MorseSequenceBuildMetrics coarse;
      FSequenceBuilder builder(complex, &coarse, false);
      compare(workers == 1 ? builder.build_flooding_reduction_kernel()
                          : builder.build_flooding_reduction_kernel_parallel(workers));
      assert(coarse.reduction_kernel_setup_nanoseconds > 0);
      assert(coarse.reduction_kernel_level_wall_nanoseconds > 0);
      assert(coarse.reduction_kernel_rounds == 0);
      assert(coarse.reduction_kernel_local_candidate_visits == 0);
      assert(coarse.reduction_kernel_facet_execution_nanoseconds == 0);
      if (workers > 1 && complex.num_levels() > 1) {
        assert(coarse.reduction_kernel_max_parallel_levels > 1);
        assert(coarse.reduction_kernel_level_chunks > 0);
      }
    }
  };

  // Exercise both word boundaries and the packed/sparse cutoff. Combining
  // these components also reuses worker scratch across different bucket sizes.
  FilteredSimplicialComplex multilevel;
  std::uint32_t offset = 0;
  std::size_t level = 0;
  for (std::size_t count : {63, 64, 65, 127, 128, 129, 255}) {
    const std::size_t vertices = count < 127 ? 6 : (count < 255 ? 7 : 8);
    FilteredSimplicialComplex single_level;
    for (std::size_t mask = 1; mask < (std::size_t{1} << vertices); ++mask) {
      std::vector<morseframes::VertexId> simplex;
      for (std::size_t v = 0; v < vertices; ++v) {
        if ((mask & (std::size_t{1} << v)) != 0) {
          simplex.push_back(offset + static_cast<std::uint32_t>(v));
        }
      }
      single_level.add_simplex(simplex, 0.0);
      multilevel.add_simplex(simplex, static_cast<double>(level));
    }
    const auto extras = count - ((std::size_t{1} << vertices) - 1);
    for (std::size_t i = 0; i < extras; ++i) {
      const auto v = offset + static_cast<std::uint32_t>(vertices + i);
      single_level.add_simplex({v}, 0.0);
      multilevel.add_simplex({v}, static_cast<double>(level));
    }
    check(single_level);
    offset += static_cast<std::uint32_t>(vertices + extras);
    ++level;
  }
  check(multilevel);

  // Shared high-dimensional facets exercise protected faces, tied weights,
  // perforations and repeated kernel rounds against the sparse oracle.
  for (unsigned seed = 0; seed < 12; ++seed) {
    std::mt19937 rng(seed);
    FilteredSimplicialComplex complex;
    std::vector<double> weights(9, 0.0);
    if (seed % 2 != 0) {
      for (std::size_t v = 0; v < weights.size(); ++v) {
        weights[v] = static_cast<double>((v * 13 + seed) % 4);
      }
    }
    for (std::size_t facet = 0; facet < 5; ++facet) {
      std::vector<morseframes::VertexId> vertices{0, 1, 2, 3, 4, 5, 6, 7, 8};
      std::shuffle(vertices.begin(), vertices.end(), rng);
      vertices.resize(4 + seed % 2);
      add_weighted_closure(complex, vertices, weights, 0.0);
    }
    check(complex);
  }
}

void test_reduction_kernel_linear_sparse_incidence() {
  const auto check = [](FilteredSimplicialComplex complex,
                        std::size_t max_closure_size) {
    complex.finalize();
    assert(complex.size() > 128);
    morseframes::MorseSequenceBuildMetrics sequential_metrics;
    const auto expected = FSequenceBuilder(complex, &sequential_metrics)
                              .build_flooding_reduction_kernel();
    morseframes::validate_morse_sequence(complex, expected);
    for (std::size_t workers : {1, 2, 4, 8}) {
      morseframes::MorseSequenceBuildMetrics metrics;
      const auto actual = FSequenceBuilder(complex, &metrics)
                              .build_flooding_reduction_kernel_parallel(workers);
      morseframes::validate_morse_sequence(complex, actual);
      assert(expected.steps().size() == actual.steps().size());
      for (std::size_t i = 0; i < expected.steps().size(); ++i) {
        const auto& a = expected.steps()[i];
        const auto& b = actual.steps()[i];
        assert(a.type == b.type && a.sigma == b.sigma && a.tau == b.tau &&
               a.level == b.level);
      }
      // These assertions detect a return to all-pairs incidence, independently
      // of wall-clock noise. Sparse entries are visited at most once per facet.
      assert(metrics.reduction_kernel_incidence_cell_visits <=
             max_closure_size * metrics.reduction_kernel_facet_kernels);
      assert(metrics.reduction_kernel_incidence_cell_visits ==
             sequential_metrics.reduction_kernel_incidence_cell_visits);
      assert(metrics.reduction_kernel_essential_parallel_tasks == 0);
      if (workers > 1 && complex.num_levels() == 1) {
        // 136 triangle closures are below the local-execution threshold;
        // graph scans and the larger high-dimensional closures exceed it.
        assert((metrics.reduction_kernel_parallel_batches > 0) ==
               (max_closure_size != 7));
        assert(metrics.reduction_kernel_parallel_batches <=
               metrics.reduction_kernel_rounds);
        assert(metrics.reduction_kernel_facet_parallel_tasks <=
               workers * metrics.reduction_kernel_parallel_batches);
        assert(metrics.reduction_kernel_facet_parallel_tasks <
               metrics.reduction_kernel_facet_kernels);
      }
    }
  };

  // More than two facets share a face: incidence must saturate at two, and
  // remain correct as reductions expose lower-dimensional facets over rounds.
  for (std::size_t facet_vertices : {2, 3, 5}) {
    for (bool multiple_levels : {false, true}) {
      FilteredSimplicialComplex complex;
      std::vector<double> values(140, 0.0);
      for (std::uint32_t v = 4; v < 140; ++v) {
        if (multiple_levels) {
          values[v] = static_cast<double>(1 + v % 3);
        }
        std::vector<morseframes::VertexId> facet;
        for (std::size_t shared = 0; shared + 1 < facet_vertices; ++shared) {
          facet.push_back(static_cast<morseframes::VertexId>(shared));
        }
        facet.push_back(v);
        add_weighted_closure(complex, facet, values);
      }
      // An isolated vertex exercises incidence for a dimension-zero facet.
      complex.add_simplex({140}, 0.0);
      check(complex, (std::size_t{1} << facet_vertices) - 1);
    }
  }
}

void test_reduction_kernel_batched_facets() {
  // Unequal cells exercise chunk tails, result ordering, and both inline and
  // overflow result storage, including dimensions above three.
  for (std::size_t facet_count :
       {1, 2, 3, 7, 8, 9, 31, 32, 33, 65, 86, 87, 129, 257}) {
    FilteredSimplicialComplex complex;
    const std::vector<double> values(6 * facet_count + 1, 0.0);
    std::size_t closure_work = 1; // The additional isolated vertex.
    for (std::size_t i = 0; i < facet_count; ++i) {
      std::vector<morseframes::VertexId> facet;
      for (std::size_t j = 0; j < 2 + i % 5; ++j) {
        facet.push_back(static_cast<morseframes::VertexId>(6 * i + j));
      }
      closure_work += (std::size_t{1} << facet.size()) - 1;
      add_weighted_closure(complex, facet, values);
    }
    complex.add_simplex(
        {static_cast<morseframes::VertexId>(6 * facet_count)}, 0.0);
    complex.finalize();
    const auto expected =
        FSequenceBuilder(complex).build_flooding_reduction_kernel();
    for (std::size_t workers : {1, 2, 4, 8}) {
      const auto compare = [&](const auto& actual) {
        morseframes::validate_morse_sequence(complex, actual);
        assert(expected.steps().size() == actual.steps().size());
        for (std::size_t i = 0; i < expected.steps().size(); ++i) {
          const auto& a = expected.steps()[i];
          const auto& b = actual.steps()[i];
          assert(a.type == b.type && a.sigma == b.sigma && a.tau == b.tau &&
                 a.level == b.level);
        }
      };
      compare(FSequenceBuilder(complex)
                  .build_flooding_reduction_kernel_parallel(workers));
      for (bool detailed : {false, true}) {
        morseframes::MorseSequenceBuildMetrics metrics;
        compare(FSequenceBuilder(complex, &metrics, detailed)
                    .build_flooding_reduction_kernel_parallel(workers));
        assert(metrics.reduction_kernel_parallel_batches <=
               metrics.reduction_kernel_rounds);
        assert(metrics.reduction_kernel_facet_parallel_tasks <=
               workers * metrics.reduction_kernel_parallel_batches);
        assert(metrics.reduction_kernel_max_parallel_facets <= workers);
        // All disjoint facets collapse in the first round; remaining isolated
        // roots are cheap. Small graph-only fixtures are also below threshold.
        const auto tasks = std::min(workers, closure_work / 1024);
        assert(metrics.reduction_kernel_facet_parallel_tasks ==
               (detailed && tasks > 1 ? tasks : 0));
      }
    }
  }
}

void test_reduction_kernel_facet_work_scheduling() {
  // Check exact threshold boundaries, partial worker budgets, and shrinking
  // rounds without introducing a second policy in the metrics-free path.
  for (std::size_t work : {1024, 2047, 2048, 2049, 3072, 8192}) {
    FilteredSimplicialComplex original;
    const std::size_t edges = (work - 7) / 3;
    const std::size_t extras = (work - 7) % 3;
    const std::vector<double> values(edges + extras + 3, 0.0);
    add_weighted_closure(original, {0, 1, 2}, values);
    for (std::size_t i = 0; i < edges; ++i) {
      add_weighted_closure(
          original, {0, static_cast<morseframes::VertexId>(i + 3)}, values);
    }
    for (std::size_t i = 0; i < extras; ++i) {
      original.add_simplex({static_cast<morseframes::VertexId>(edges + 3 + i)}, 0.0);
    }
    original.finalize();
    assert(original.size() > 128);
    for (bool cached : {false, true}) {
      auto complex = original;
      if (cached) complex.prepare_same_level_closure_cache();
      morseframes::MorseSequenceBuildMetrics sequential_metrics;
      const auto expected = FSequenceBuilder(complex, &sequential_metrics)
                                .build_flooding_reduction_kernel();
      const auto compare = [&](const auto& actual) {
        morseframes::validate_morse_sequence(complex, actual);
        assert(expected.steps().size() == actual.steps().size());
        for (std::size_t i = 0; i < expected.steps().size(); ++i) {
          const auto& a = expected.steps()[i];
          const auto& b = actual.steps()[i];
          assert(a.type == b.type && a.sigma == b.sigma && a.tau == b.tau &&
                 a.level == b.level);
        }
      };
      for (std::size_t workers : {1, 2, 4, 8}) {
        compare(FSequenceBuilder(complex)
                    .build_flooding_reduction_kernel_parallel(workers));
        for (bool detailed : {false, true}) {
          morseframes::MorseSequenceBuildMetrics metrics;
          compare(FSequenceBuilder(complex, &metrics, detailed)
                      .build_flooding_reduction_kernel_parallel(workers));
          const auto tasks = std::min(workers, work / 1024);
          assert(metrics.reduction_kernel_facet_parallel_tasks ==
                 (detailed && tasks > 1 ? tasks : 0));
          assert(metrics.reduction_kernel_parallel_batches ==
                 (detailed && tasks > 1 ? 1 : 0));
          if (detailed) {
            assert(metrics.reduction_kernel_incidence_cell_visits ==
                   sequential_metrics.reduction_kernel_incidence_cell_visits);
            assert(metrics.reduction_kernel_facet_cell_visits ==
                   sequential_metrics.reduction_kernel_facet_cell_visits);
            assert(metrics.reduction_kernel_local_candidate_visits ==
                   sequential_metrics.reduction_kernel_local_candidate_visits);
          }
        }
      }
    }
  }

  // Two large closures can justify parallelism where hundreds of small ones
  // do not. Include dimensions above three and overflow cell/event storage.
  for (std::size_t extras : {0, 2}) {
    FilteredSimplicialComplex complex;
    const std::vector<double> values(20 + extras, 0.0);
    add_weighted_closure(complex, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, values);
    add_weighted_closure(complex, {10, 11, 12, 13, 14, 15, 16, 17, 18, 19}, values);
    for (std::size_t i = 0; i < extras; ++i) {
      complex.add_simplex({static_cast<morseframes::VertexId>(20 + i)}, 0.0);
    }
    complex.finalize();
    const auto expected = FSequenceBuilder(complex).build_flooding_reduction_kernel();
    for (bool detailed : {false, true}) {
      morseframes::MorseSequenceBuildMetrics metrics;
      const auto actual = FSequenceBuilder(complex, &metrics, detailed)
                              .build_flooding_reduction_kernel_parallel(8);
      morseframes::validate_morse_sequence(complex, actual);
      assert(expected.steps().size() == actual.steps().size());
      for (std::size_t i = 0; i < expected.steps().size(); ++i) {
        const auto& a = expected.steps()[i];
        const auto& b = actual.steps()[i];
        assert(a.type == b.type && a.sigma == b.sigma && a.tau == b.tau &&
               a.level == b.level);
      }
      // Initial work is 2 * 1023 plus the isolated vertices; later rounds
      // contain only the two roots and the additional isolated vertices.
      assert(metrics.reduction_kernel_facet_parallel_tasks ==
             (detailed && extras == 2 ? 2 : 0));
    }
  }
}

void test_reduction_kernel_discovery_granularity() {
  // Connected stars collapse to one vertex in a single reducing round. Thus
  // only the first discovery can dispatch tasks; later rounds must use the
  // remaining active count rather than the original bucket size.
  for (std::size_t size : {129, 8191, 8192, 8193, 12288, 32769}) {
    FilteredSimplicialComplex original;
    // One filled triangle enables closure storage, keeping this scheduling
    // test independent of the graph-only local cell's full-bucket scan.
    const std::size_t edges = (size - 7) / 2;
    const std::vector<double> values(edges + 4, 0.0);
    add_weighted_closure(original, {0, 1, 2}, values);
    for (morseframes::VertexId v = 3; v < edges + 3; ++v) {
      add_weighted_closure(original, {0, v}, values);
    }
    if (size % 2 == 0) {
      original.add_simplex({static_cast<morseframes::VertexId>(edges + 3)}, 0.0);
    }
    original.finalize();
    assert(original.size() == size);
    for (bool cached : {false, true}) {
      auto complex = original;
      if (cached) complex.prepare_same_level_closure_cache();
      morseframes::MorseSequenceBuildMetrics sequential_metrics;
      const auto expected = FSequenceBuilder(complex, &sequential_metrics)
                                .build_flooding_reduction_kernel();
      const auto compare = [&](const auto& actual) {
        morseframes::validate_morse_sequence(complex, actual);
        assert(expected.steps().size() == actual.steps().size());
        for (std::size_t i = 0; i < expected.steps().size(); ++i) {
          const auto& a = expected.steps()[i];
          const auto& b = actual.steps()[i];
          assert(a.type == b.type && a.sigma == b.sigma && a.tau == b.tau &&
                 a.level == b.level);
        }
      };
      for (std::size_t workers : {1, 2, 4, 8}) {
        morseframes::MorseSequenceBuildMetrics metrics;
        compare(FSequenceBuilder(complex, &metrics)
                    .build_flooding_reduction_kernel_parallel(workers));
        const auto tasks = std::min(workers, size / 4096);
        assert(metrics.reduction_kernel_facet_discovery_parallel_tasks ==
               (tasks > 1 ? tasks : 0));
        assert(metrics.reduction_kernel_facet_discovery_coboundary_visits ==
               sequential_metrics.reduction_kernel_facet_discovery_coboundary_visits);
        assert(metrics.reduction_kernel_incidence_cell_visits ==
               sequential_metrics.reduction_kernel_incidence_cell_visits);
        assert(metrics.reduction_kernel_rounds > 1);
      }
      compare(FSequenceBuilder(complex).build_flooding_reduction_kernel_parallel(8));
    }
  }
}

void test_reduction_kernel_discovery_failure_drains_tasks() {
  struct FailingDiscoveryView : FilteredSimplicialComplex {
    mutable std::atomic<std::size_t> failures{0};
    const std::vector<morseframes::SimplexId>& coboundary(
        morseframes::SimplexId) const {
      // A graph has no cached closure build, so the first coboundary access is
      // in discovery. Each static chunk throws, including after a peer fails.
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
      ++failures;
      throw std::runtime_error("discovery failure");
    }
  };
  FailingDiscoveryView complex;
  const std::vector<double> values(16385, 0.0);
  for (morseframes::VertexId v = 1; v < values.size(); ++v) {
    add_weighted_closure(complex, {0, v}, values);
  }
  complex.finalize();
  assert(complex.size() == 32769);
  for (std::size_t workers : {2, 4, 8}) {
    for (bool detailed : {false, true}) {
      complex.failures = 0;
      auto executor = std::make_shared<morseframes::BoundedTaskExecutor>(workers);
      morseframes::ReductionKernelExecutionOptions options;
      options.policy = morseframes::ReductionKernelExecutionPolicy::Parallel;
      options.collect_metrics = detailed;
      morseframes::ReductionKernelWorkspace<FailingDiscoveryView> workspace(
          complex, options, executor);
      bool propagated = false;
      try {
        (void)workspace.compute_level_isolated(0);
      } catch (const std::runtime_error& error) {
        propagated = std::string(error.what()) == "discovery failure";
      }
      assert(propagated);
      assert(complex.failures == workers); // Before implicit teardown joins.
      auto following = executor->submit([]() { return 17; });
      assert(executor->get(following) == 17);
    }
  }
}

void test_reduction_kernel_facet_failure_drains_tasks() {
  struct FailingFacetView : FilteredSimplicialComplex {
    mutable std::atomic<std::size_t> failures{0};
    const std::vector<morseframes::VertexId>& vertices(
        morseframes::SimplexId) const {
      // In a sparse graph, only the local facet kernel needs vertices().
      // Every submitted task fails on its first chunk, even after a peer fails.
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
      ++failures;
      throw std::runtime_error("facet failure");
    }
  };
  FailingFacetView complex;
  const std::vector<double> values(81, 0.0);
  for (morseframes::VertexId v = 1; v < values.size(); ++v) {
    add_weighted_closure(complex, {0, v}, values);
  }
  complex.finalize();
  assert(complex.size() > 128);
  for (std::size_t workers : {2, 4, 8}) {
    for (bool detailed : {false, true}) {
      complex.failures = 0;
      auto executor = std::make_shared<morseframes::BoundedTaskExecutor>(workers);
      morseframes::ReductionKernelExecutionOptions options;
      options.policy = morseframes::ReductionKernelExecutionPolicy::Parallel;
      options.collect_metrics = detailed;
      morseframes::ReductionKernelWorkspace<FailingFacetView> workspace(
          complex, options, executor);
      bool propagated = false;
      try {
        (void)workspace.compute_level_isolated(0);
      } catch (const std::runtime_error& error) {
        propagated = std::string(error.what()) == "facet failure";
      }
      assert(propagated);
      // Check BEFORE workspace/executor teardown can implicitly join workers.
      assert(complex.failures == workers);
      auto following = executor->submit([]() { return 17; });
      assert(executor->get(following) == 17);
    }
  }
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

void test_complex_construction_contract() {
  std::mt19937 rng(17);
  for (unsigned dimension = 1; dimension <= 7; ++dimension) {
    std::vector<std::vector<morseframes::VertexId>> faces;
    for (unsigned mask = 1; mask < (1u << (dimension + 1)); ++mask) {
      std::vector<morseframes::VertexId> face;
      for (unsigned i = 0; i <= dimension; ++i)
        if (mask & (1u << i)) face.push_back(i * 1000003u);
      faces.push_back(face);
    }
    std::sort(faces.begin(), faces.end());
    auto shuffled = faces;
    std::shuffle(shuffled.begin(), shuffled.end(), rng);
    FilteredSimplicialComplex a, b;
    for (auto face : faces) a.add_simplex(face, double(face.size() - 1));
    for (auto face : shuffled) {
      const double value = double(face.size() - 1);
      std::reverse(face.begin(), face.end());
      b.add_simplex(face, value);
      b.add_simplex(face, value + 0.5e-12); // Keep the original value within tolerance.
    }
    a.finalize();
    morseframes::ComplexConstructionMetrics metrics;
    metrics.boundaries_seconds = -1;
    b.finalize_with_metrics(metrics);
    assert(metrics.reset_seconds >= 0 && metrics.index_and_simplices_seconds >= 0);
    assert(metrics.levels_seconds >= 0 && metrics.boundaries_seconds >= 0);
    assert(metrics.coboundaries_seconds >= 0 && metrics.orders_and_buckets_seconds >= 0);
    assert(a.size() == faces.size() && b.size() == a.size());
    assert(a.filtration_order() == b.filtration_order());
    assert(a.level_values() == b.level_values());
    for (morseframes::SimplexId id = 0; id < a.size(); ++id) {
      assert(a.vertices(id) == faces[id] && a.vertices(id) == b.vertices(id));
      assert(a.filtration(id) == b.filtration(id));
      assert(a.level(id) == b.level(id) && a.dimension(id) == b.dimension(id));
      assert(a.boundary(id) == b.boundary(id) && a.coboundary(id) == b.coboundary(id));
      assert(b.find_simplex(faces[id]) == id);
      auto reversed = faces[id]; std::reverse(reversed.begin(), reversed.end());
      assert(b.find_simplex(reversed) == id);
      assert(a.simplices_of_level(a.level(id)) == b.simplices_of_level(b.level(id)));
    }
    auto copied = b;
    auto moved = std::move(copied);
    b.prepare_same_level_closure_cache();
    b.add_simplex({4000000000u}, 10);
    assert(!b.has_same_level_closure_cache());
    b.finalize();
    assert(b.size() == a.size() + 1);
    assert(moved.size() == a.size());
    auto check = [&moved, &faces] {
      for (morseframes::SimplexId id = 0; id < faces.size(); ++id)
        assert(moved.find_simplex(faces[id]) == id);
    };
    auto future = std::async(std::launch::async, check); check(); future.get();
    moved.finalize();
    check();
    assert(moved.find_simplex({4000000000u}) == morseframes::kInvalidSimplex);
    bool rejected = false;
    try { moved.add_simplex({0}, 2); } catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
  }
  FilteredSimplicialComplex empty;
  morseframes::ComplexConstructionMetrics metrics;
  bool rejected = false;
  try { empty.finalize_with_metrics(metrics); } catch (const std::invalid_argument&) { rejected = true; }
  assert(rejected);
}

void test_compact_simplex_lookup() {
  using Vertices = std::vector<morseframes::VertexId>;
  const auto missing = morseframes::kInvalidSimplex;
  const auto largest = std::numeric_limits<morseframes::VertexId>::max();
  const Vertices vertices{0, 4, 17, 1000, 4000000000u, largest};
  const std::vector<Vertices> cells{{0, 17, 4000000000u}, {17, 1000, largest},
                                   {4}, {4000000000u, largest}};
  std::map<Vertices, double> oracle;
  FilteredSimplicialComplex complex;
  assert(complex.find_simplex({}) == missing);
  assert(complex.find_simplex({0}) == missing);
  for (const auto& cell : cells) {
    for (unsigned mask = 1; mask < (1u << cell.size()); ++mask) {
      Vertices face;
      for (unsigned i = 0; i < cell.size(); ++i) if (mask & (1u << i)) face.push_back(cell[i]);
      std::sort(face.begin(), face.end());
      oracle[face] = 0;
      std::reverse(face.begin(), face.end());
      complex.add_simplex(face, 0);
    }
  }
  assert(complex.find_simplex({0}) == missing); // Pending insertion is not finalization.
  const auto check = [&](const auto& view) {
    for (unsigned mask = 0; mask < (1u << vertices.size()); ++mask) {
      Vertices face;
      for (unsigned i = 0; i < vertices.size(); ++i) if (mask & (1u << i)) face.push_back(vertices[i]);
      const auto it = oracle.find(face);
      const auto expected = it == oracle.end() ? missing
          : static_cast<morseframes::SimplexId>(std::distance(oracle.begin(), it));
      assert(view.find_simplex(face) == expected);
      std::reverse(face.begin(), face.end());
      assert(view.find_simplex(face) == expected);
    }
    for (auto vertex : {1u, 3u, 5u, 999u, 4000000001u})
      assert(view.find_simplex({vertex}) == missing);
    bool rejected = false;
    try { (void)view.find_simplex({17, 17}); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
  };
  complex.finalize(); check(complex);
  auto copied = complex;
  auto moved = std::move(copied);
  check(moved);
  auto concurrent = std::async(std::launch::async, [&] { check(moved); });
  check(moved); concurrent.get();
  complex.add_simplex({2}, 0);
  complex.add_simplex({2, 17}, 0);
  assert(complex.find_simplex({2}) == missing);
  oracle[{2}] = 0; oracle[{2, 17}] = 0;
  complex.finalize(); check(complex); // New prefix inserted between existing ranges.
  assert(complex.find_simplex({2}) != missing);
  assert(complex.find_simplex({17, 2}) != missing);
  complex.finalize(); check(complex);
  // Missing faces must still be detected, including a missing singleton at
  // the start of an otherwise present first-vertex range.
  for (const auto& bad : {std::vector<Vertices>{{2, 9}, {9}},
                          std::vector<Vertices>{{2}, {9}, {2, 9, 17}}}) {
    FilteredSimplicialComplex invalid;
    for (const auto& face : bad) invalid.add_simplex(face, 0);
    bool rejected = false;
    try { invalid.finalize(); } catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
  }
}

void test_boundary_and_filtration_order() {
  using Vertices = std::vector<morseframes::VertexId>;
  using Id = morseframes::SimplexId;
  const auto largest = std::numeric_limits<morseframes::VertexId>::max();
  std::mt19937 rng(509);
  for (unsigned dimension = 0; dimension <= 7; ++dimension) {
    for (unsigned trial = 0; trial < 4; ++trial) {
      Vertices vertices{0, 2, 19, 1000, 1000003, 2000007, 3000017,
                        4000000000u, largest - 2, largest};
      std::map<morseframes::VertexId, double> weights;
      for (auto vertex : vertices) weights[vertex] = trial % 2 ? double(rng() % 4) : 0;
      std::vector<Vertices> cells{
          Vertices(vertices.begin(), vertices.begin() + dimension + 1),
          Vertices(vertices.begin() + 1, vertices.begin() + dimension + 2),
          {largest - 2, largest}, {vertices[dimension + 1]}, {largest}};
      std::map<Vertices, double> faces;
      for (const auto& cell : cells) {
        for (unsigned mask = 1; mask < (1u << cell.size()); ++mask) {
          Vertices face;
          double value = 0;
          for (unsigned i = 0; i < cell.size(); ++i) if (mask & (1u << i)) {
            face.push_back(cell[i]); value = std::max(value, weights.at(cell[i]));
          }
          // Also exercise monotone filtrations which are not vertex lower stars.
          if (trial >= 2) value += double(face.size() - 1);
          faces[face] = value;
        }
      }
      std::vector<Vertices> ordered, shuffled;
      std::map<Vertices, Id> ids;
      for (const auto& entry : faces) {
        ids[entry.first] = static_cast<Id>(ordered.size());
        ordered.push_back(entry.first);
      }
      shuffled = ordered;
      std::shuffle(shuffled.begin(), shuffled.end(), rng);
      FilteredSimplicialComplex complex;
      for (auto face : shuffled) {
        const auto value = faces.at(face);
        std::reverse(face.begin(), face.end());
        complex.add_simplex(face, value);
      }
      for (unsigned pass = 0; pass < 2; ++pass) {
        if (pass) {
          morseframes::ComplexConstructionMetrics metrics;
          complex.finalize_with_metrics(metrics);
        } else complex.finalize();
        std::vector<std::vector<Id>> expected_coboundaries(ordered.size());
        std::vector<Id> expected_order;
        std::vector<double> levels;
        for (Id id = 0; id < ordered.size(); ++id) {
          assert(complex.vertices(id) == ordered[id]);
          std::vector<Id> boundary;
          if (ordered[id].size() > 1) {
            for (std::size_t removed = 0; removed < ordered[id].size(); ++removed) {
              auto face = ordered[id]; face.erase(face.begin() + removed);
              boundary.push_back(ids.at(face));
              expected_coboundaries[ids.at(face)].push_back(id);
            }
          }
          assert(complex.boundary(id) == boundary); // Deletion order, not sorted IDs.
          expected_order.push_back(id);
          levels.push_back(faces.at(ordered[id]));
        }
        for (Id id = 0; id < ordered.size(); ++id)
          assert(complex.coboundary(id) == expected_coboundaries[id]);
        std::sort(levels.begin(), levels.end());
        levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
        assert(complex.level_values() == levels);
        // Independent legacy comparator: filtration, dimension, vertex vector.
        std::sort(expected_order.begin(), expected_order.end(), [&](Id a, Id b) {
          const auto av = faces.at(ordered[a]), bv = faces.at(ordered[b]);
          if (av != bv) return av < bv;
          if (ordered[a].size() != ordered[b].size()) return ordered[a].size() < ordered[b].size();
          return ordered[a] < ordered[b];
        });
        assert(complex.filtration_order() == expected_order);
        for (std::size_t level = 0; level < levels.size(); ++level) {
          std::vector<Id> bucket;
          for (Id id : expected_order)
            if (faces.at(ordered[id]) == levels[level]) bucket.push_back(id);
          assert(complex.simplices_of_level(level) == bucket);
        }
      }
    }
  }
  // Detect absent and non-monotone facets in both the general lookup (removed
  // first vertex) and the reused-range lookup (every other removed vertex).
  const Vertices cell{2, 19, 4000000000u, largest};
  for (unsigned removed = 0; removed < cell.size(); ++removed) {
    auto bad_face = cell; bad_face.erase(bad_face.begin() + removed);
    for (bool missing : {false, true}) {
      FilteredSimplicialComplex invalid;
      for (unsigned mask = 1; mask < (1u << cell.size()); ++mask) {
        Vertices face;
        for (unsigned i = 0; i < cell.size(); ++i) if (mask & (1u << i)) face.push_back(cell[i]);
        if (missing && face == bad_face) continue;
        invalid.add_simplex(face, face == bad_face ? 1 : 0);
      }
      bool rejected = false;
      try { invalid.finalize(); } catch (const std::invalid_argument&) { rejected = true; }
      assert(rejected);
    }
  }
}

void test_bulk_lower_star_construction() {
  using Cells = std::vector<std::vector<morseframes::VertexId>>;
  const auto legacy = [](auto& complex, const auto& values, const auto& cells) {
    for (const auto& cell : cells) {
      for (std::size_t mask = 1; mask < (std::size_t{1} << cell.size()); ++mask) {
        std::vector<morseframes::VertexId> face;
        double value = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < cell.size(); ++i) if (mask & (std::size_t{1} << i)) {
          face.push_back(cell[i]);
          value = std::max(value, values[cell[i]]);
        }
        complex.add_simplex(std::move(face), value);
      }
    }
  };
  const auto same = [](const auto& a, const auto& b) {
    assert(a.size() == b.size());
    assert(a.filtration_order() == b.filtration_order());
    assert(a.level_values() == b.level_values());
    for (std::size_t i = 0; i < a.num_levels(); ++i) {
      assert(std::signbit(a.level_values()[i]) == std::signbit(b.level_values()[i]));
      assert(a.simplices_of_level(i) == b.simplices_of_level(i));
    }
    for (morseframes::SimplexId id = 0; id < a.size(); ++id) {
      assert(a.vertices(id) == b.vertices(id));
      assert(a.filtration(id) == b.filtration(id));
      assert(std::signbit(a.filtration(id)) == std::signbit(b.filtration(id)));
      assert(a.level(id) == b.level(id) && a.dimension(id) == b.dimension(id));
      assert(a.boundary(id) == b.boundary(id) && a.coboundary(id) == b.coboundary(id));
      assert(a.find_simplex(a.vertices(id)) == b.find_simplex(a.vertices(id)));
    }
  };
  std::mt19937 rng(290);
  for (unsigned dimension = 0; dimension <= 7; ++dimension) {
    for (unsigned trial = 0; trial < 8; ++trial) {
      // Shared faces, repeated/reversed cells, a lower-dimensional maximal cell
      // and an isolated vertex. Exercise the vector-backed path above 3D too.
      Cells cells(2);
      for (unsigned i = 0; i <= dimension; ++i) {
        cells[0].push_back(i);
        cells[1].push_back(i + 1);
      }
      cells.push_back({dimension + 2, dimension + 1});
      cells.push_back({dimension + 3});
      cells.push_back(cells[0]);
      for (auto& cell : cells) std::shuffle(cell.begin(), cell.end(), rng);
      std::shuffle(cells.begin(), cells.end(), rng);
      std::vector<double> values(dimension + 4);
      for (auto& value : values) value = trial % 2 ? double(rng() % 5) - 2 : -0.0;
      values[0] = +0.0;
      FilteredSimplicialComplex a, b, diagnostic;
      // Preexisting entries interleave with every dimension's sorted batch.
      a.add_simplex({dimension + 3}, values.back() + 0.5e-12);
      b.add_simplex({dimension + 3}, values.back() + 0.5e-12);
      diagnostic.add_simplex({dimension + 3}, values.back() + 0.5e-12);
      legacy(a, values, cells);
      morseframes::add_lower_star_cells(b, values, cells);
      morseframes::LowerStarConstructionMetrics metrics;
      metrics.generated_faces = 99999;
      morseframes::add_lower_star_cells_with_metrics(diagnostic, values, cells, metrics);
      std::size_t generated = 0;
      for (const auto& cell : cells) generated += (std::size_t{1} << cell.size()) - 1;
      assert(metrics.generated_faces == generated);
      assert(metrics.validation_seconds >= 0 && metrics.enumeration_seconds >= 0);
      assert(metrics.sort_and_dedup_seconds >= 0 && metrics.insertion_seconds >= 0);
      a.finalize(); b.finalize(); diagnostic.finalize();
      assert(metrics.unique_faces_submitted == b.size());
      same(a, b); same(a, diagnostic);
      const auto expected = FSequenceBuilder<FilteredSimplicialComplex>(a).build_f_max();
      const auto actual = FSequenceBuilder<FilteredSimplicialComplex>(b).build_f_max();
      assert(expected.steps().size() == actual.steps().size());
      for (std::size_t i = 0; i < expected.steps().size(); ++i) {
        const auto& x = expected.steps()[i]; const auto& y = actual.steps()[i];
        assert(x.type == y.type && x.sigma == y.sigma && x.tau == y.tau && x.level == y.level);
      }
      auto copied = b;
      auto moved = std::move(copied);
      same(a, moved);
      b.prepare_same_level_closure_cache();
      morseframes::add_lower_star_cells(b, values, cells);
      assert(!b.has_same_level_closure_cache());
      b.finalize(); same(a, b);
      b.prepare_same_level_closure_cache();
      morseframes::add_lower_star_cells_with_metrics(b, values, {}, metrics);
      assert(metrics.generated_faces == 0 && metrics.unique_faces_submitted == 0);
      assert(b.has_same_level_closure_cache());
      b.finalize(); same(a, b);
    }
  }
  // All validation precedes mutation, even if an earlier cell would be valid.
  FilteredSimplicialComplex invalid;
  invalid.add_simplex({0}, 0); invalid.finalize();
  const auto original = invalid;
  invalid.prepare_same_level_closure_cache();
  for (const auto& cells : {Cells{{1}, {}}, Cells{{1}, {0, 0}}, Cells{{1}, {2}}}) {
    bool rejected = false;
    try { morseframes::add_lower_star_cells(invalid, {0, 1}, cells); }
    catch (const std::exception&) { rejected = true; }
    assert(rejected && invalid.has_same_level_closure_cache());
    same(original, invalid);
  }
  bool rejected = false;
  try { morseframes::add_lower_star_cells(invalid, {0, std::nan("")}, {{1}}); }
  catch (const std::invalid_argument&) { rejected = true; }
  assert(rejected && invalid.has_same_level_closure_cache());
  invalid.finalize(); same(original, invalid); // No hidden pending insertions.
  rejected = false;
  try { morseframes::add_lower_star_cells(invalid, {2}, {{0}}); }
  catch (const std::invalid_argument&) { rejected = true; }
  assert(rejected);
  invalid.finalize(); same(original, invalid);
  const double infinity = std::numeric_limits<double>::infinity();
  FilteredSimplicialComplex a, b;
  legacy(a, std::vector<double>{-infinity, infinity}, Cells{{1, 0}, {0, 1}});
  morseframes::add_lower_star_cells(b, {-infinity, infinity}, {{1, 0}, {0, 1}});
  a.finalize(); b.finalize(); same(a, b);
  assert(morseframes::detail::LowerStarComplexBuilder::combinations(10, 3) == 120);
  rejected = false;
  try { (void)morseframes::detail::LowerStarComplexBuilder::combinations(1000, 500); }
  catch (const std::length_error&) { rejected = true; }
  assert(rejected);
}

}  // namespace

int main() {
  test_compact_simplex_lookup();
  test_boundary_and_filtration_order();
  test_bulk_lower_star_construction();
  test_complex_construction_contract();
  test_bounded_task_executor();
  test_boundary_and_coboundary();
  test_inverse_annotation_store();
  test_field_annotation_store();
  test_monotonicity_rejection();
  test_simplex_tree_builder_gudhi_style_insert();
  test_simplex_tree_builder_strict_duplicate_rejection();
  test_simplex_tree_builder_explicit_insert_can_be_nonclosed();
  test_filtered_complex_from_simplex_tree_adapter();
  test_f_sequence_builder_accepts_simplex_tree_view();
  test_process_lower_stars_triangle_boundary();
  test_process_lower_stars_workspace_and_dimensions();
  test_one_vertex();
  test_reducer_skips_initially_zero_boundaries();
  test_two_vertices_joined_by_later_edge();
  test_same_level_edge_has_no_off_diagonal_pair();
  test_non_flooding_f_sequence_is_reduced_in_flooding_order();
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
  test_flooding_reduction_kernel_on_shared_facets();
  test_reduction_kernel_packed_core_matches_sparse_cache();
  test_reduction_kernel_lightweight_initialization();
  test_reduction_kernel_lightweight_persistence();
  test_reduction_kernel_linear_sparse_incidence();
  test_reduction_kernel_batched_facets();
  test_reduction_kernel_facet_work_scheduling();
  test_reduction_kernel_discovery_granularity();
  test_reduction_kernel_discovery_failure_drains_tasks();
  test_reduction_kernel_facet_failure_drains_tasks();
  test_instrumentation_metrics();

  std::cout << "All Morse persistence prototype tests passed.\n";
  return 0;
}
