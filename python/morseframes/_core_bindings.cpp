#include <cctype>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "morseframes/coreference_persistence.hpp"
#include "morseframes/filtered_complex.hpp"
#include "morseframes/morse_sequence.hpp"
#include "morseframes/reference_persistence.hpp"
#include "morseframes/simplex_tree_builder.hpp"
#include "morseframes/standard_persistence.hpp"

namespace nb = nanobind;

namespace {

using SimplexInput = std::pair<std::vector<morseframes::VertexId>, double>;
using Clock = std::chrono::steady_clock;
using FiniteBarcodeKey = std::tuple<std::uint16_t, double, double>;
using EssentialBarcodeKey = std::tuple<std::uint16_t, double>;

std::uint64_t elapsed_nanoseconds(Clock::time_point start, Clock::time_point end) {
  return static_cast<std::uint64_t>(
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
}

struct PyFilteredComplex {
  morseframes::FilteredSimplicialComplex complex;
  bool finalized = false;

  void add_simplex(std::vector<morseframes::VertexId> vertices, double filtration) {
    complex.add_simplex(std::move(vertices), filtration);
    finalized = false;
  }

  void finalize() {
    complex.finalize();
    finalized = true;
  }
};

struct PySimplexTreeBuilder {
  explicit PySimplexTreeBuilder(const std::string& duplicate_policy)
      : builder(duplicate_policy) {}

  bool insert(std::vector<morseframes::VertexId> vertices,
              double filtration,
              bool include_faces) {
    return include_faces
        ? builder.insert(std::move(vertices), filtration)
        : builder.insert_simplex_only(std::move(vertices), filtration);
  }

  bool insert_simplex_only(std::vector<morseframes::VertexId> vertices, double filtration) {
    return builder.insert_simplex_only(std::move(vertices), filtration);
  }

  PyFilteredComplex to_filtered_complex(bool finalize) const {
    PyFilteredComplex result;
    result.complex = builder.to_filtered_complex(finalize);
    result.finalized = finalize;
    return result;
  }

  PyFilteredComplex finalize(bool clear_builder) {
    PyFilteredComplex result;
    result.complex = builder.finalize(clear_builder);
    result.finalized = true;
    return result;
  }

  morseframes::SimplexTreeBuilder builder;
};

struct PyMorseSequence {
  explicit PyMorseSequence(morseframes::MorseSequence sequence) : sequence(std::move(sequence)) {}

  morseframes::MorseSequence sequence;
};

struct PyReferenceMap {
  explicit PyReferenceMap(std::vector<morseframes::Annotation> references)
      : references(std::move(references)) {}

  std::vector<morseframes::Annotation> references;
};

struct PyFieldReferenceMap {
  explicit PyFieldReferenceMap(std::vector<morseframes::FieldAnnotation> references)
      : references(std::move(references)) {}

  std::vector<morseframes::FieldAnnotation> references;
};

struct PyMorseReferenceFrame {
  explicit PyMorseReferenceFrame(morseframes::MorseReferenceFrame frame)
      : sequence(std::move(frame.sequence)), references(std::move(frame.references)) {}

  morseframes::MorseSequence sequence;
  std::vector<morseframes::Annotation> references;
};

struct PyMorseCoreferenceFrame {
  explicit PyMorseCoreferenceFrame(morseframes::MorseCoreferenceFrame frame)
      : sequence(std::move(frame.sequence)), coreferences(std::move(frame.coreferences)) {}

  morseframes::MorseSequence sequence;
  std::vector<morseframes::Annotation> coreferences;
};

void require_finalized(const PyFilteredComplex& complex) {
  if (!complex.finalized) {
    throw std::logic_error("FilteredComplex.finalize() must be called before querying.");
  }
}

nb::list annotation_table_to_python(const std::vector<morseframes::Annotation>& annotations) {
  nb::list result;
  for (const auto& annotation : annotations) {
    nb::list py_annotation;
    for (morseframes::CriticalId label : annotation) {
      py_annotation.append(label);
    }
    result.append(py_annotation);
  }
  return result;
}

nb::list annotation_to_python(const morseframes::Annotation& annotation) {
  nb::list result;
  for (morseframes::CriticalId label : annotation) {
    result.append(label);
  }
  return result;
}

morseframes::Annotation annotation_from_python(nb::handle annotation) {
  nb::sequence labels = nb::cast<nb::sequence>(annotation);
  morseframes::Annotation result;
  result.reserve(nb::len(labels));
  for (nb::handle label : labels) {
    result.push_back(nb::cast<morseframes::CriticalId>(label));
  }
  return result;
}

std::vector<morseframes::Annotation> annotation_table_from_python(nb::sequence annotations) {
  std::vector<morseframes::Annotation> result;
  result.reserve(nb::len(annotations));
  for (nb::handle annotation : annotations) {
    result.push_back(annotation_from_python(annotation));
  }
  return result;
}

nb::list field_annotation_table_to_python(
    const std::vector<morseframes::FieldAnnotation>& annotations) {
  nb::list result;
  for (const auto& annotation : annotations) {
    nb::list py_annotation;
    for (const auto& entry : annotation) {
      py_annotation.append(nb::make_tuple(entry.label, entry.coefficient));
    }
    result.append(py_annotation);
  }
  return result;
}

nb::list sequence_steps_to_python(const PyMorseSequence& sequence) {
  nb::list result;
  for (const auto& step : sequence.sequence.steps()) {
    if (step.type == morseframes::MorseStepType::Critical) {
      result.append(nb::make_tuple("critical", step.sigma, nb::none(), step.level));
    } else {
      result.append(nb::make_tuple("regular_pair", step.sigma, step.tau, step.level));
    }
  }
  return result;
}

nb::list paired_with_to_python(const morseframes::FilteredSimplicialComplex& complex,
                               const morseframes::MorseSequence& sequence) {
  std::vector<morseframes::SimplexId> paired_with(
      complex.size(), static_cast<morseframes::SimplexId>(morseframes::kInvalidSimplex));
  for (const auto& step : sequence.steps()) {
    if (step.type != morseframes::MorseStepType::RegularPair) {
      continue;
    }
    paired_with.at(step.sigma) = step.tau;
    paired_with.at(step.tau) = step.sigma;
  }

  nb::list result;
  for (morseframes::SimplexId simplex = 0; simplex < complex.size(); ++simplex) {
    if (paired_with[simplex] == morseframes::kInvalidSimplex) {
      result.append(nb::none());
    } else {
      result.append(paired_with[simplex]);
    }
  }
  return result;
}

nb::list paired_with_to_python(const PyFilteredComplex& complex, const PyMorseSequence& sequence) {
  require_finalized(complex);
  return paired_with_to_python(complex.complex, sequence.sequence);
}

nb::list finite_pairs_to_python(const std::vector<morseframes::PersistencePair>& pairs) {
  nb::list result;
  for (const auto& pair : pairs) {
    result.append(nb::make_tuple(pair.birth,
                                 pair.death,
                                 pair.dimension,
                                 pair.birth_value,
                                 pair.death_value));
  }
  return result;
}

nb::list essential_to_python(const std::vector<morseframes::EssentialInterval>& intervals) {
  nb::list result;
  for (const auto& interval : intervals) {
    result.append(nb::make_tuple(interval.birth, interval.dimension, interval.birth_value));
  }
  return result;
}

nb::dict diagram_to_python(const morseframes::PersistenceDiagram& diagram) {
  nb::dict result;
  result["finite_pairs"] = finite_pairs_to_python(diagram.finite_pairs);
  result["essential"] = essential_to_python(diagram.essential);
  return result;
}

nb::dict reduction_metrics_to_python(const morseframes::MorseReferenceReductionMetrics& metrics) {
  nb::dict inverse;
  inverse["initial_nonempty_annotations"] =
      metrics.inverse_store.initial_nonempty_annotations;
  inverse["initial_total_annotation_size"] =
      metrics.inverse_store.initial_total_annotation_size;
  inverse["initial_max_annotation_size"] =
      metrics.inverse_store.initial_max_annotation_size;
  inverse["initial_inverse_list_entries"] =
      metrics.inverse_store.initial_inverse_list_entries;
  inverse["remove_candidate_scans"] = metrics.inverse_store.remove_candidate_scans;
  inverse["remove_applied"] = metrics.inverse_store.remove_applied;
  inverse["remove_total_annotation_size"] =
      metrics.inverse_store.remove_total_annotation_size;
  inverse["remove_max_annotation_size"] =
      metrics.inverse_store.remove_max_annotation_size;
  inverse["xor_candidate_scans"] = metrics.inverse_store.xor_candidate_scans;
  inverse["xor_applied"] = metrics.inverse_store.xor_applied;
  inverse["xor_changed_labels"] = metrics.inverse_store.xor_changed_labels;
  inverse["xor_total_input_size"] = metrics.inverse_store.xor_total_input_size;
  inverse["xor_total_output_size"] = metrics.inverse_store.xor_total_output_size;
  inverse["xor_max_input_size"] = metrics.inverse_store.xor_max_input_size;
  inverse["xor_max_output_size"] = metrics.inverse_store.xor_max_output_size;
  inverse["xor_inserted_labels"] = metrics.inverse_store.xor_inserted_labels;
  inverse["xor_removed_labels"] = metrics.inverse_store.xor_removed_labels;
  inverse["inverse_list_appends"] = metrics.inverse_store.inverse_list_appends;

  nb::dict result;
  result["working_set_size"] = metrics.working_set_size;
  result["critical_count"] = metrics.critical_count;
  result["reducer_setup_nanoseconds"] = metrics.reducer_setup_nanoseconds;
  result["reducer_compute_nanoseconds"] = metrics.reducer_compute_nanoseconds;
  result["boundary_plan_face_scans"] = metrics.boundary_plan_face_scans;
  result["boundary_annotation_candidate_criticals"] =
      metrics.boundary_annotation_candidate_criticals;
  result["boundary_annotation_zero_skipped_criticals"] =
      metrics.boundary_annotation_zero_skipped_criticals;
  result["boundary_annotation_zero_skipped_faces"] =
      metrics.boundary_annotation_zero_skipped_faces;
  result["boundary_annotation_xors"] = metrics.boundary_annotation_xors;
  result["boundary_annotation_total_input_size"] =
      metrics.boundary_annotation_total_input_size;
  result["boundary_annotation_total_output_size"] =
      metrics.boundary_annotation_total_output_size;
  result["boundary_annotation_max_size"] = metrics.boundary_annotation_max_size;
  result["boundary_annotation_max_output_size"] =
      metrics.boundary_annotation_max_output_size;
  result["pivot_eliminations"] = metrics.pivot_eliminations;
  result["finite_pairs"] = metrics.finite_pairs;
  result["essential_intervals"] = metrics.essential_intervals;
  result["inverse_store"] = inverse;
  return result;
}

nb::dict frame_metrics_to_python(const morseframes::MorseReferenceFrameMetrics& metrics) {
  nb::dict result;
  result["remaining_cofaces_nanoseconds"] = metrics.remaining_cofaces_nanoseconds;
  result["sequence_total_nanoseconds"] = metrics.sequence_total_nanoseconds;
  result["sequence_core_nanoseconds"] = metrics.sequence_core_nanoseconds;
  result["reference_update_nanoseconds"] = metrics.reference_update_nanoseconds;
  result["reduction_plan_nanoseconds"] = metrics.reduction_plan_nanoseconds;
  result["release_cleanup_nanoseconds"] = metrics.release_cleanup_nanoseconds;
  result["working_set_pack_nanoseconds"] = metrics.working_set_pack_nanoseconds;
  result["local_index_nanoseconds"] = metrics.local_index_nanoseconds;
  result["sequence_init_nanoseconds"] = metrics.sequence_init_nanoseconds;
  result["sequence_candidate_seed_nanoseconds"] =
      metrics.sequence_candidate_seed_nanoseconds;
  result["sequence_candidate_loop_nanoseconds"] =
      metrics.sequence_candidate_loop_nanoseconds;
  result["sequence_emit_nanoseconds"] = metrics.sequence_emit_nanoseconds;
  result["sequence_callback_nanoseconds"] = metrics.sequence_callback_nanoseconds;
  result["sequence_replay_nanoseconds"] = metrics.sequence_replay_nanoseconds;
  result["sequence_candidate_pushes"] = metrics.sequence_candidate_pushes;
  result["sequence_candidate_pops"] = metrics.sequence_candidate_pops;
  result["sequence_stale_candidate_skips"] = metrics.sequence_stale_candidate_skips;
  result["sequence_level_mismatch_skips"] = metrics.sequence_level_mismatch_skips;
  result["sequence_regular_pairs"] = metrics.sequence_regular_pairs;
  result["sequence_criticals"] = metrics.sequence_criticals;
  result["final_live_nonempty_annotations"] = metrics.final_live_nonempty_annotations;
  result["final_live_total_annotation_size"] = metrics.final_live_total_annotation_size;
  result["peak_live_nonempty_annotations"] = metrics.peak_live_nonempty_annotations;
  result["peak_live_total_annotation_size"] = metrics.peak_live_total_annotation_size;
  result["released_annotations"] = metrics.released_annotations;
  result["released_total_annotation_size"] = metrics.released_total_annotation_size;
  return result;
}

nb::dict reduction_result_to_python(const morseframes::MorseReferenceReductionResult& result) {
  nb::dict py_result;
  py_result["diagram"] = diagram_to_python(result.diagram);
  py_result["metrics"] = reduction_metrics_to_python(result.metrics);
  return py_result;
}

nb::dict reference_profile_to_python(const morseframes::MorseReferenceProfile& profile) {
  nb::dict result;
  result["num_simplices"] = profile.num_simplices;
  result["num_levels"] = profile.num_levels;
  result["num_critical_simplices"] = profile.num_critical_simplices;
  result["num_regular_pairs"] = profile.num_regular_pairs;
  result["frame_metrics"] = frame_metrics_to_python(profile.frame_metrics);
  result["metrics"] = reduction_metrics_to_python(profile.reduction_metrics);
  result["estimated_reducer_work"] = profile.estimated_reducer_work;
  return result;
}

std::vector<FiniteBarcodeKey> finite_barcode_signature(
    const morseframes::PersistenceDiagram& diagram) {
  std::vector<FiniteBarcodeKey> signature;
  signature.reserve(diagram.finite_pairs.size());
  for (const auto& pair : morseframes::off_diagonal_pairs(diagram)) {
    signature.emplace_back(pair.dimension, pair.birth_value, pair.death_value);
  }
  std::sort(signature.begin(), signature.end());
  return signature;
}

std::vector<EssentialBarcodeKey> essential_barcode_signature(
    const morseframes::PersistenceDiagram& diagram) {
  std::vector<EssentialBarcodeKey> signature;
  signature.reserve(diagram.essential.size());
  for (const auto& interval : diagram.essential) {
    signature.emplace_back(interval.dimension, interval.birth_value);
  }
  std::sort(signature.begin(), signature.end());
  return signature;
}

nb::dict core_benchmark_result_to_python(
    const morseframes::MorseReferenceReductionResult& morse_result,
    const morseframes::PersistenceDiagram& standard_diagram,
    std::size_t num_critical_simplices,
    std::uint64_t sequence_nanoseconds,
    std::uint64_t reference_nanoseconds,
    std::uint64_t morse_reduction_nanoseconds,
    std::uint64_t morse_nanoseconds,
    std::uint64_t standard_nanoseconds,
    const morseframes::MorseReferenceFrameMetrics& frame_metrics,
    const std::string& frame_mode) {
  const auto morse_finite = finite_barcode_signature(morse_result.diagram);
  const auto standard_finite = finite_barcode_signature(standard_diagram);
  const auto morse_essential = essential_barcode_signature(morse_result.diagram);
  const auto standard_essential = essential_barcode_signature(standard_diagram);

  nb::dict result;
  result["num_critical_simplices"] = num_critical_simplices;
  result["frame_mode"] = frame_mode;
  result["sequence_nanoseconds"] = sequence_nanoseconds;
  result["reference_nanoseconds"] = reference_nanoseconds;
  result["morse_reduction_nanoseconds"] = morse_reduction_nanoseconds;
  result["morse_nanoseconds"] = morse_nanoseconds;
  result["standard_nanoseconds"] = standard_nanoseconds;
  result["finite_signature_matches"] = morse_finite == standard_finite;
  result["essential_signature_matches"] = morse_essential == standard_essential;
  result["morse_finite_interval_count"] = morse_finite.size();
  result["morse_essential_interval_count"] = morse_essential.size();
  result["standard_finite_interval_count"] = standard_finite.size();
  result["standard_essential_interval_count"] = standard_essential.size();
  result["frame_metrics"] = frame_metrics_to_python(frame_metrics);
  result["metrics"] = reduction_metrics_to_python(morse_result.metrics);
  return result;
}

nb::dict direction_benchmark_result_to_python(
    const std::string& direction,
    const morseframes::MorseReferenceReductionResult& morse_result,
    const morseframes::PersistenceDiagram& standard_diagram,
    std::size_t num_simplices,
    std::size_t num_levels,
    std::size_t num_critical_simplices,
    std::uint64_t frame_nanoseconds,
    std::uint64_t morse_reduction_nanoseconds,
    std::uint64_t morse_nanoseconds,
    std::uint64_t standard_nanoseconds) {
  const auto morse_finite = finite_barcode_signature(morse_result.diagram);
  const auto morse_essential = essential_barcode_signature(morse_result.diagram);
  const auto standard_finite = finite_barcode_signature(standard_diagram);
  const auto standard_essential = essential_barcode_signature(standard_diagram);

  nb::dict result;
  result["direction"] = direction;
  result["num_simplices"] = num_simplices;
  result["num_levels"] = num_levels;
  result["num_critical_simplices"] = num_critical_simplices;
  result["frame_nanoseconds"] = frame_nanoseconds;
  result["morse_reduction_nanoseconds"] = morse_reduction_nanoseconds;
  result["morse_nanoseconds"] = morse_nanoseconds;
  result["standard_nanoseconds"] = standard_nanoseconds;
  result["finite_signature_matches"] = morse_finite == standard_finite;
  result["essential_signature_matches"] = morse_essential == standard_essential;
  result["morse_finite_interval_count"] = morse_finite.size();
  result["morse_essential_interval_count"] = morse_essential.size();
  result["standard_finite_interval_count"] = standard_finite.size();
  result["standard_essential_interval_count"] = standard_essential.size();
  result["metrics"] = reduction_metrics_to_python(morse_result.metrics);
  return result;
}

std::string normalize_sequence_algorithm(std::string algorithm) {
  for (char& character : algorithm) {
    if (character == '_') {
      character = '-';
    } else {
      character = static_cast<char>(
          std::tolower(static_cast<unsigned char>(character)));
    }
  }
  if (algorithm == "f-sequence" || algorithm == "saturated-f-sequence") {
    return "saturated";
  }
  if (algorithm == "f-max" ||
      algorithm == "paper-max" ||
      algorithm == "max-s-f" ||
      algorithm == "max-sf" ||
      algorithm == "maximal-f-sequence") {
    return "f-max";
  }
  if (algorithm == "f-min" ||
      algorithm == "paper-min" ||
      algorithm == "min-s-f" ||
      algorithm == "min-sf" ||
      algorithm == "minimal-f-sequence") {
    return "f-min";
  }
  if (algorithm == "plateau" || algorithm == "plateau-greedy-f-sequence") {
    return "plateau-greedy";
  }
  if (algorithm == "coreduction" ||
      algorithm == "same-level-coreduction" ||
      algorithm == "coreduction-f-sequence" ||
      algorithm == "reduction-reverse" ||
      algorithm == "same-level-reduction-reverse") {
    return "same-level-reduction";
  }
  if (algorithm == "flooding-maximal" || algorithm == "maximal-flooding") {
    return "flooding-max";
  }
  if (algorithm == "flooding-minimal" || algorithm == "minimal-flooding") {
    return "flooding-min";
  }
  if (algorithm == "flooding" ||
      algorithm == "flooding-minmax" || algorithm == "flooding-min-max" ||
      algorithm == "min-max" || algorithm == "minmax") {
    return "flooding-minmax";
  }
  if (algorithm == "flooding-maxmin" || algorithm == "flooding-max-min" ||
      algorithm == "max-min" || algorithm == "maxmin") {
    return "flooding-maxmin";
  }
  return algorithm;
}

bool is_implemented_sequence_algorithm(const std::string& algorithm) {
  return algorithm == "saturated" ||
         algorithm == "f-max" ||
         algorithm == "f-min" ||
         algorithm == "plateau-greedy" ||
         algorithm == "same-level-reduction" ||
         algorithm == "flooding-max" ||
         algorithm == "flooding-min" ||
         algorithm == "flooding-minmax" ||
         algorithm == "flooding-maxmin";
}

PyMorseSequence build_sequence(const PyFilteredComplex& complex,
                               const std::string& algorithm) {
  require_finalized(complex);
  const std::string normalized = normalize_sequence_algorithm(algorithm);
  if (normalized == "flooding" || normalized == "stack-flooding") {
    throw std::logic_error(
        "This Morse sequence algorithm is reserved for a future implementation.");
  }
  if (normalized == "plateau-greedy") {
    return PyMorseSequence{morseframes::FSequenceBuilder(complex.complex).build_plateau_greedy()};
  }
  if (normalized == "same-level-reduction") {
    return PyMorseSequence{
        morseframes::FSequenceBuilder(complex.complex).build_same_level_reduction()};
  }
  if (normalized == "f-max") {
    return PyMorseSequence{morseframes::FSequenceBuilder(complex.complex).build_f_max()};
  }
  if (normalized == "f-min") {
    return PyMorseSequence{morseframes::FSequenceBuilder(complex.complex).build_f_min()};
  }
  if (normalized == "flooding-max") {
    return PyMorseSequence{morseframes::FSequenceBuilder(complex.complex).build_flooding_max()};
  }
  if (normalized == "flooding-min") {
    return PyMorseSequence{morseframes::FSequenceBuilder(complex.complex).build_flooding_min()};
  }
  if (normalized == "flooding-minmax") {
    return PyMorseSequence{morseframes::FSequenceBuilder(complex.complex).build_flooding_minmax()};
  }
  if (normalized == "flooding-maxmin") {
    return PyMorseSequence{morseframes::FSequenceBuilder(complex.complex).build_flooding_maxmin()};
  }
  if (normalized != "saturated") {
    throw std::invalid_argument("Unknown Morse sequence algorithm: " + algorithm);
  }
  return PyMorseSequence{morseframes::FSequenceBuilder(complex.complex).build_saturated()};
}

PyMorseReferenceFrame build_sequence_and_reference_map(const PyFilteredComplex& complex,
                                                       const std::string& algorithm) {
  require_finalized(complex);
  const std::string normalized = normalize_sequence_algorithm(algorithm);
  if (normalized == "flooding" || normalized == "stack-flooding") {
    throw std::logic_error(
        "This Morse sequence algorithm is reserved for a future implementation.");
  }
  if (normalized == "plateau-greedy") {
    return PyMorseReferenceFrame{
        morseframes::MorseReferenceFrameBuilder(complex.complex).build_plateau_greedy()};
  }
  if (normalized == "same-level-reduction") {
    return PyMorseReferenceFrame{
        morseframes::MorseReferenceFrameBuilder(complex.complex).build_same_level_reduction()};
  }
  if (normalized == "f-max") {
    return PyMorseReferenceFrame{
        morseframes::MorseReferenceFrameBuilder(complex.complex).build_f_max()};
  }
  if (normalized == "f-min") {
    return PyMorseReferenceFrame{
        morseframes::MorseReferenceFrameBuilder(complex.complex).build_f_min()};
  }
  if (normalized == "flooding-max") {
    return PyMorseReferenceFrame{
        morseframes::MorseReferenceFrameBuilder(complex.complex).build_flooding_max()};
  }
  if (normalized == "flooding-min") {
    return PyMorseReferenceFrame{
        morseframes::MorseReferenceFrameBuilder(complex.complex).build_flooding_min()};
  }
  if (normalized == "flooding-minmax") {
    return PyMorseReferenceFrame{
        morseframes::MorseReferenceFrameBuilder(complex.complex).build_flooding_minmax()};
  }
  if (normalized == "flooding-maxmin") {
    return PyMorseReferenceFrame{
        morseframes::MorseReferenceFrameBuilder(complex.complex).build_flooding_maxmin()};
  }
  if (normalized != "saturated") {
    throw std::invalid_argument("Unknown Morse sequence algorithm: " + algorithm);
  }
  return PyMorseReferenceFrame{
      morseframes::MorseReferenceFrameBuilder(complex.complex).build_saturated()};
}

PyMorseCoreferenceFrame build_sequence_and_coreference_map(const PyFilteredComplex& complex,
                                                           const std::string& algorithm) {
  require_finalized(complex);
  const std::string normalized = normalize_sequence_algorithm(algorithm);
  if (normalized == "flooding" || normalized == "stack-flooding") {
    throw std::logic_error(
        "This Morse sequence algorithm is reserved for a future implementation.");
  }
  if (normalized == "same-level-reduction") {
    return PyMorseCoreferenceFrame{
        morseframes::MorseCoreferenceFrameBuilder(complex.complex).build_same_level_reduction()};
  }

  if (!is_implemented_sequence_algorithm(normalized)) {
    throw std::invalid_argument("Unknown Morse sequence algorithm: " + algorithm);
  }
  morseframes::MorseSequence sequence =
      normalized == "plateau-greedy"
          ? morseframes::FSequenceBuilder(complex.complex).build_plateau_greedy()
      : normalized == "same-level-reduction"
          ? morseframes::FSequenceBuilder(complex.complex).build_same_level_reduction()
      : normalized == "f-max"
          ? morseframes::FSequenceBuilder(complex.complex).build_f_max()
      : normalized == "f-min"
          ? morseframes::FSequenceBuilder(complex.complex).build_f_min()
      : normalized == "flooding-max"
          ? morseframes::FSequenceBuilder(complex.complex).build_flooding_max()
      : normalized == "flooding-min"
          ? morseframes::FSequenceBuilder(complex.complex).build_flooding_min()
      : normalized == "flooding-minmax"
          ? morseframes::FSequenceBuilder(complex.complex).build_flooding_minmax()
      : normalized == "flooding-maxmin"
          ? morseframes::FSequenceBuilder(complex.complex).build_flooding_maxmin()
          : morseframes::FSequenceBuilder(complex.complex).build_saturated();
  auto coreferences =
      morseframes::MorseCoreferenceComputer(complex.complex, sequence).compute_full_coreferences();
  return PyMorseCoreferenceFrame{
      morseframes::MorseCoreferenceFrame{std::move(sequence), std::move(coreferences)}};
}

PyReferenceMap compute_reference_map_object(const PyFilteredComplex& complex,
                                            const PyMorseSequence& sequence) {
  require_finalized(complex);
  auto references =
      morseframes::MorseReferenceComputer(complex.complex, sequence.sequence).compute_full_references();
  return PyReferenceMap{std::move(references)};
}

PyFieldReferenceMap compute_reference_map_modp_object(const PyFilteredComplex& complex,
                                                      const PyMorseSequence& sequence,
                                                      std::uint32_t modulus) {
  require_finalized(complex);
  auto references =
      morseframes::MorseFieldReferenceComputer(complex.complex, sequence.sequence, modulus)
          .compute_full_references();
  return PyFieldReferenceMap{std::move(references)};
}

PyReferenceMap compute_coreference_map_object(const PyFilteredComplex& complex,
                                              const PyMorseSequence& sequence) {
  require_finalized(complex);
  auto coreferences =
      morseframes::MorseCoreferenceComputer(complex.complex, sequence.sequence)
          .compute_full_coreferences();
  return PyReferenceMap{std::move(coreferences)};
}

PyFieldReferenceMap compute_coreference_map_modp_object(const PyFilteredComplex& complex,
                                                        const PyMorseSequence& sequence,
                                                        std::uint32_t modulus) {
  require_finalized(complex);
  auto coreferences =
      morseframes::MorseFieldCoreferenceComputer(complex.complex, sequence.sequence, modulus)
          .compute_full_coreferences();
  return PyFieldReferenceMap{std::move(coreferences)};
}

nb::list reference_map_to_python(const PyReferenceMap& references) {
  return annotation_table_to_python(references.references);
}

nb::list compute_reference_map(const PyFilteredComplex& complex, const PyMorseSequence& sequence) {
  return reference_map_to_python(compute_reference_map_object(complex, sequence));
}

nb::list field_reference_map_to_python(const PyFieldReferenceMap& references) {
  return field_annotation_table_to_python(references.references);
}

nb::list compute_reference_map_modp(const PyFilteredComplex& complex,
                                    const PyMorseSequence& sequence,
                                    std::uint32_t modulus) {
  return field_reference_map_to_python(
      compute_reference_map_modp_object(complex, sequence, modulus));
}

nb::list compute_coreference_map(const PyFilteredComplex& complex, const PyMorseSequence& sequence) {
  return reference_map_to_python(compute_coreference_map_object(complex, sequence));
}

nb::list compute_coreference_map_modp(const PyFilteredComplex& complex,
                                      const PyMorseSequence& sequence,
                                      std::uint32_t modulus) {
  return field_reference_map_to_python(
      compute_coreference_map_modp_object(complex, sequence, modulus));
}

nb::dict compute_morse_persistence(const PyFilteredComplex& complex,
                                   const PyMorseSequence& sequence) {
  require_finalized(complex);
  auto references =
      morseframes::MorseReferenceComputer(complex.complex, sequence.sequence).compute_full_references();
  auto diagram =
      morseframes::MorseReferencePersistenceReducer(complex.complex, sequence.sequence, references)
          .compute();
  return diagram_to_python(diagram);
}

nb::dict compute_morse_persistence_modp(const PyFilteredComplex& complex,
                                        const PyMorseSequence& sequence,
                                        std::uint32_t modulus) {
  require_finalized(complex);
  auto diagram =
      morseframes::compute_morse_reference_prime_field_persistence(
          complex.complex, sequence.sequence, modulus);
  return diagram_to_python(diagram);
}

nb::dict compute_morse_coreference_persistence(const PyFilteredComplex& complex,
                                               const PyMorseSequence& sequence) {
  require_finalized(complex);
  auto coreferences =
      morseframes::MorseCoreferenceComputer(complex.complex, sequence.sequence)
          .compute_full_coreferences();
  auto diagram =
      morseframes::MorseCoreferencePersistenceReducer(
          complex.complex, sequence.sequence, coreferences)
          .compute();
  return diagram_to_python(diagram);
}

nb::dict compute_morse_coreference_persistence_modp(const PyFilteredComplex& complex,
                                                    const PyMorseSequence& sequence,
                                                    std::uint32_t modulus) {
  require_finalized(complex);
  auto diagram =
      morseframes::compute_morse_coreference_prime_field_persistence(
          complex.complex, sequence.sequence, modulus);
  return diagram_to_python(diagram);
}

nb::dict reduce_morse_persistence(const PyFilteredComplex& complex,
                                  const PyMorseSequence& sequence,
                                  nb::sequence references) {
  require_finalized(complex);
  auto reference_table = annotation_table_from_python(references);
  auto diagram =
      morseframes::MorseReferencePersistenceReducer(complex.complex, sequence.sequence, reference_table)
          .compute();
  return diagram_to_python(diagram);
}

nb::dict reduce_morse_persistence_with_metrics(
    const PyFilteredComplex& complex,
    const PyMorseSequence& sequence,
    nb::sequence references) {
  require_finalized(complex);
  auto reference_table = annotation_table_from_python(references);
  const auto setup_started = Clock::now();
  morseframes::MorseReferencePersistenceReducer reducer(
      complex.complex, sequence.sequence, reference_table);
  const auto compute_started = Clock::now();
  auto result = reducer.compute_with_metrics();
  const auto finished = Clock::now();
  result.metrics.reducer_setup_nanoseconds =
      elapsed_nanoseconds(setup_started, compute_started);
  result.metrics.reducer_compute_nanoseconds =
      elapsed_nanoseconds(compute_started, finished);
  return reduction_result_to_python(result);
}

nb::dict reduce_morse_persistence_object(const PyFilteredComplex& complex,
                                         const PyMorseSequence& sequence,
                                         const PyReferenceMap& references) {
  require_finalized(complex);
  auto diagram =
      morseframes::MorseReferencePersistenceReducer(
          complex.complex, sequence.sequence, references.references)
          .compute();
  return diagram_to_python(diagram);
}

nb::dict reduce_morse_coreference_persistence_object(const PyFilteredComplex& complex,
                                                     const PyMorseSequence& sequence,
                                                     const PyReferenceMap& coreferences) {
  require_finalized(complex);
  auto diagram =
      morseframes::MorseCoreferencePersistenceReducer(
          complex.complex, sequence.sequence, coreferences.references)
          .compute();
  return diagram_to_python(diagram);
}

nb::dict reduce_morse_persistence_object_with_metrics(const PyFilteredComplex& complex,
                                                      const PyMorseSequence& sequence,
                                                      const PyReferenceMap& references) {
  require_finalized(complex);
  const auto setup_started = Clock::now();
  morseframes::MorseReferencePersistenceReducer reducer(
      complex.complex, sequence.sequence, references.references);
  const auto compute_started = Clock::now();
  auto result = reducer.compute_with_metrics();
  const auto finished = Clock::now();
  result.metrics.reducer_setup_nanoseconds =
      elapsed_nanoseconds(setup_started, compute_started);
  result.metrics.reducer_compute_nanoseconds =
      elapsed_nanoseconds(compute_started, finished);
  return reduction_result_to_python(result);
}

nb::dict reduce_morse_reference_frame_object(const PyFilteredComplex& complex,
                                             const PyMorseReferenceFrame& frame) {
  require_finalized(complex);
  auto diagram =
      morseframes::MorseReferencePersistenceReducer(complex.complex, frame.sequence, frame.references)
          .compute();
  return diagram_to_python(diagram);
}

nb::dict reduce_morse_reference_frame_object_with_metrics(const PyFilteredComplex& complex,
                                                          const PyMorseReferenceFrame& frame) {
  require_finalized(complex);
  const auto setup_started = Clock::now();
  morseframes::MorseReferencePersistenceReducer reducer(complex.complex, frame.sequence, frame.references);
  const auto compute_started = Clock::now();
  auto result = reducer.compute_with_metrics();
  const auto finished = Clock::now();
  result.metrics.reducer_setup_nanoseconds =
      elapsed_nanoseconds(setup_started, compute_started);
  result.metrics.reducer_compute_nanoseconds =
      elapsed_nanoseconds(compute_started, finished);
  return reduction_result_to_python(result);
}

nb::dict reduce_morse_coreference_frame_object(const PyFilteredComplex& complex,
                                               const PyMorseCoreferenceFrame& frame) {
  require_finalized(complex);
  auto diagram =
      morseframes::MorseCoreferencePersistenceReducer(
          complex.complex, frame.sequence, frame.coreferences)
          .compute();
  return diagram_to_python(diagram);
}

nb::dict reduce_morse_coreference_frame_object_with_metrics(
    const PyFilteredComplex& complex,
    const PyMorseCoreferenceFrame& frame) {
  require_finalized(complex);
  const auto setup_started = Clock::now();
  morseframes::MorseCoreferencePersistenceReducer reducer(
      complex.complex, frame.sequence, frame.coreferences);
  const auto compute_started = Clock::now();
  auto result = reducer.compute_with_metrics();
  const auto finished = Clock::now();
  result.metrics.reducer_setup_nanoseconds =
      elapsed_nanoseconds(setup_started, compute_started);
  result.metrics.reducer_compute_nanoseconds =
      elapsed_nanoseconds(compute_started, finished);
  return reduction_result_to_python(result);
}

nb::dict compute_standard_persistence(const PyFilteredComplex& complex) {
  require_finalized(complex);
  return diagram_to_python(morseframes::compute_standard_z2_persistence(complex.complex));
}

nb::dict compute_standard_persistence_modp(const PyFilteredComplex& complex,
                                           std::uint32_t modulus) {
  require_finalized(complex);
  return diagram_to_python(
      morseframes::compute_standard_prime_field_persistence(complex.complex, modulus));
}

nb::dict benchmark_morse_reference_core(const PyFilteredComplex& complex,
                                        const std::string& algorithm,
                                        const std::string& frame_mode) {
  require_finalized(complex);
  const std::string normalized_algorithm = normalize_sequence_algorithm(algorithm);
  if (normalized_algorithm == "flooding" || normalized_algorithm == "stack-flooding") {
    throw std::logic_error(
        "This Morse sequence algorithm is reserved for a future implementation.");
  }
  if (!is_implemented_sequence_algorithm(normalized_algorithm)) {
    throw std::invalid_argument("Unknown Morse sequence algorithm: " + algorithm);
  }

  std::string normalized_frame_mode = frame_mode;
  for (char& character : normalized_frame_mode) {
    character = static_cast<char>(
        std::tolower(static_cast<unsigned char>(character)));
  }

  const auto morse_started = Clock::now();
  morseframes::MorseReferenceReductionResult morse_result;
  std::uint64_t sequence_nanoseconds = 0;
  std::uint64_t reference_nanoseconds = 0;
  std::uint64_t morse_reduction_nanoseconds = 0;
  std::size_t num_critical_simplices = 0;
  morseframes::MorseReferenceFrameMetrics frame_metrics;

  if (normalized_frame_mode == "fused") {
    morseframes::MorseReferenceReductionInput input =
        normalized_algorithm == "plateau-greedy"
            ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                  .build_plateau_greedy_reduction_input()
        : normalized_algorithm == "same-level-reduction"
            ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                  .build_same_level_reduction_reduction_input()
        : normalized_algorithm == "f-max"
            ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                  .build_f_max_reduction_input()
        : normalized_algorithm == "f-min"
            ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                  .build_f_min_reduction_input()
        : normalized_algorithm == "flooding-max"
            ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                  .build_flooding_max_reduction_input()
        : normalized_algorithm == "flooding-min"
            ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                  .build_flooding_min_reduction_input()
        : normalized_algorithm == "flooding-minmax"
            ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                  .build_flooding_minmax_reduction_input()
        : normalized_algorithm == "flooding-maxmin"
            ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                  .build_flooding_maxmin_reduction_input()
            : morseframes::MorseReferenceFrameBuilder(complex.complex)
                  .build_saturated_reduction_input();
    const auto references_finished = Clock::now();
    sequence_nanoseconds = elapsed_nanoseconds(morse_started, references_finished);
    reference_nanoseconds = 0;
    num_critical_simplices = input.sequence.critical_simplices().size();
    frame_metrics = input.frame_metrics;

    const auto setup_started = Clock::now();
    morseframes::MorseReferencePersistenceReducer reducer(complex.complex,
                                                    input.sequence,
                                                    std::move(input.reduction_plan),
                                                    std::move(input.annotations));
    const auto compute_started = Clock::now();
    morse_result = reducer.compute_with_metrics();
    const auto reduction_finished = Clock::now();
    morse_result.metrics.reducer_setup_nanoseconds =
        elapsed_nanoseconds(setup_started, compute_started);
    morse_result.metrics.reducer_compute_nanoseconds =
        elapsed_nanoseconds(compute_started, reduction_finished);
    morse_reduction_nanoseconds =
        elapsed_nanoseconds(setup_started, reduction_finished);
  } else if (normalized_frame_mode == "separate") {
    auto sequence = normalized_algorithm == "plateau-greedy"
        ? morseframes::FSequenceBuilder(complex.complex).build_plateau_greedy()
        : normalized_algorithm == "same-level-reduction"
            ? morseframes::FSequenceBuilder(complex.complex).build_same_level_reduction()
        : normalized_algorithm == "f-max"
            ? morseframes::FSequenceBuilder(complex.complex).build_f_max()
        : normalized_algorithm == "f-min"
            ? morseframes::FSequenceBuilder(complex.complex).build_f_min()
        : normalized_algorithm == "flooding-max"
            ? morseframes::FSequenceBuilder(complex.complex).build_flooding_max()
        : normalized_algorithm == "flooding-min"
            ? morseframes::FSequenceBuilder(complex.complex).build_flooding_min()
        : normalized_algorithm == "flooding-minmax"
            ? morseframes::FSequenceBuilder(complex.complex).build_flooding_minmax()
        : normalized_algorithm == "flooding-maxmin"
            ? morseframes::FSequenceBuilder(complex.complex).build_flooding_maxmin()
            : morseframes::FSequenceBuilder(complex.complex).build_saturated();
    const auto sequence_finished = Clock::now();
    auto references =
        morseframes::MorseReferenceComputer(complex.complex, sequence).compute_full_references();
    const auto references_finished = Clock::now();
    sequence_nanoseconds = elapsed_nanoseconds(morse_started, sequence_finished);
    reference_nanoseconds = elapsed_nanoseconds(sequence_finished, references_finished);
    num_critical_simplices = sequence.critical_simplices().size();

    const auto setup_started = Clock::now();
    morseframes::MorseReferencePersistenceReducer reducer(
        complex.complex, sequence, references);
    const auto compute_started = Clock::now();
    morse_result = reducer.compute_with_metrics();
    const auto reduction_finished = Clock::now();
    morse_result.metrics.reducer_setup_nanoseconds =
        elapsed_nanoseconds(setup_started, compute_started);
    morse_result.metrics.reducer_compute_nanoseconds =
        elapsed_nanoseconds(compute_started, reduction_finished);
    morse_reduction_nanoseconds =
        elapsed_nanoseconds(setup_started, reduction_finished);
  } else {
    throw std::invalid_argument("Unknown Morse frame mode: " + frame_mode);
  }

  const auto morse_finished = Clock::now();
  const auto standard_started = Clock::now();
  auto standard_diagram = morseframes::compute_standard_z2_persistence(complex.complex);
  const auto standard_finished = Clock::now();

  return core_benchmark_result_to_python(
      morse_result,
      standard_diagram,
      num_critical_simplices,
      sequence_nanoseconds,
      reference_nanoseconds,
      morse_reduction_nanoseconds,
      elapsed_nanoseconds(morse_started, morse_finished),
      elapsed_nanoseconds(standard_started, standard_finished),
      frame_metrics,
      normalized_frame_mode);
}

nb::dict profile_morse_reference_frame_core(const PyFilteredComplex& complex,
                                            const std::string& algorithm) {
  require_finalized(complex);
  const std::string normalized_algorithm = normalize_sequence_algorithm(algorithm);
  if (normalized_algorithm == "flooding" || normalized_algorithm == "stack-flooding") {
    throw std::logic_error(
        "This Morse sequence algorithm is reserved for a future implementation.");
  }
  if (!is_implemented_sequence_algorithm(normalized_algorithm)) {
    throw std::invalid_argument("Unknown Morse sequence algorithm: " + algorithm);
  }

  const auto started = Clock::now();
  morseframes::MorseReferenceReductionInput input =
      normalized_algorithm == "plateau-greedy"
          ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                .build_plateau_greedy_reduction_input()
      : normalized_algorithm == "same-level-reduction"
          ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                .build_same_level_reduction_reduction_input()
      : normalized_algorithm == "f-max"
          ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                .build_f_max_reduction_input()
      : normalized_algorithm == "f-min"
          ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                .build_f_min_reduction_input()
      : normalized_algorithm == "flooding-max"
          ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                .build_flooding_max_reduction_input()
      : normalized_algorithm == "flooding-min"
          ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                .build_flooding_min_reduction_input()
      : normalized_algorithm == "flooding-minmax"
          ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                .build_flooding_minmax_reduction_input()
      : normalized_algorithm == "flooding-maxmin"
          ? morseframes::MorseReferenceFrameBuilder(complex.complex)
                .build_flooding_maxmin_reduction_input()
          : morseframes::MorseReferenceFrameBuilder(complex.complex)
                .build_saturated_reduction_input();
  auto profile = morseframes::profile_morse_reference_reduction_input(complex.complex, input);
  const auto finished = Clock::now();

  nb::dict result = reference_profile_to_python(profile);
  result["sequence_algorithm"] = normalized_algorithm;
  result["frame_mode"] = "fused";
  result["profile_nanoseconds"] = elapsed_nanoseconds(started, finished);
  return result;
}

nb::list benchmark_coreduction_directions_core(const PyFilteredComplex& complex) {
  require_finalized(complex);

  const auto standard_started = Clock::now();
  auto standard_diagram = morseframes::compute_standard_z2_persistence(complex.complex);
  const auto standard_finished = Clock::now();
  const auto standard_nanoseconds =
      elapsed_nanoseconds(standard_started, standard_finished);

  nb::list rows;

  {
    const auto morse_started = Clock::now();
    auto input =
        morseframes::MorseReferenceFrameBuilder(complex.complex)
            .build_same_level_reduction_reduction_input();
    const auto frame_finished = Clock::now();
    const auto num_critical_simplices = input.sequence.critical_simplices().size();

    const auto setup_started = Clock::now();
    morseframes::MorseReferencePersistenceReducer reducer(complex.complex,
                                                    input.sequence,
                                                    std::move(input.reduction_plan),
                                                    std::move(input.annotations));
    const auto compute_started = Clock::now();
    auto morse_result = reducer.compute_with_metrics();
    const auto reduction_finished = Clock::now();
    morse_result.metrics.reducer_setup_nanoseconds =
        elapsed_nanoseconds(setup_started, compute_started);
    morse_result.metrics.reducer_compute_nanoseconds =
        elapsed_nanoseconds(compute_started, reduction_finished);
    rows.append(direction_benchmark_result_to_python(
        "reference",
        morse_result,
        standard_diagram,
        complex.complex.size(),
        complex.complex.num_levels(),
        num_critical_simplices,
        elapsed_nanoseconds(morse_started, frame_finished),
        elapsed_nanoseconds(setup_started, reduction_finished),
        elapsed_nanoseconds(morse_started, reduction_finished),
        standard_nanoseconds));
  }

  {
    const auto morse_started = Clock::now();
    auto frame =
        morseframes::MorseCoreferenceFrameBuilder(complex.complex).build_same_level_reduction();
    const auto frame_finished = Clock::now();
    const auto num_critical_simplices = frame.sequence.critical_simplices().size();

    const auto setup_started = Clock::now();
    morseframes::MorseCoreferencePersistenceReducer reducer(
        complex.complex, frame.sequence, frame.coreferences);
    const auto compute_started = Clock::now();
    auto morse_result = reducer.compute_with_metrics();
    const auto reduction_finished = Clock::now();
    morse_result.metrics.reducer_setup_nanoseconds =
        elapsed_nanoseconds(setup_started, compute_started);
    morse_result.metrics.reducer_compute_nanoseconds =
        elapsed_nanoseconds(compute_started, reduction_finished);
    rows.append(direction_benchmark_result_to_python(
        "coreference",
        morse_result,
        standard_diagram,
        complex.complex.size(),
        complex.complex.num_levels(),
        num_critical_simplices,
        elapsed_nanoseconds(morse_started, frame_finished),
        elapsed_nanoseconds(setup_started, reduction_finished),
        elapsed_nanoseconds(morse_started, reduction_finished),
        standard_nanoseconds));
  }

  return rows;
}

nb::dict analyze(const std::vector<SimplexInput>& simplices) {
  morseframes::FilteredSimplicialComplex complex;
  for (const auto& [vertices, filtration] : simplices) {
    complex.add_simplex(vertices, filtration);
  }
  complex.finalize();

  PyMorseReferenceFrame frame{
      morseframes::MorseReferenceFrameBuilder(complex).build_saturated()};
  auto morse_diagram =
      morseframes::MorseReferencePersistenceReducer(complex, frame.sequence, frame.references).compute();
  auto standard_diagram = morseframes::compute_standard_z2_persistence(complex);

  nb::dict result;
  result["steps"] = sequence_steps_to_python(PyMorseSequence{frame.sequence});
  result["critical_simplices"] = frame.sequence.critical_simplices();
  result["critical_index_of_simplex"] = frame.sequence.critical_index_of_simplex();
  result["paired_with"] = paired_with_to_python(complex, frame.sequence);
  result["references"] = annotation_table_to_python(frame.references);
  result["morse"] = diagram_to_python(morse_diagram);
  result["standard"] = diagram_to_python(standard_diagram);
  return result;
}

}  // namespace

NB_MODULE(_morse_core, m) {
  m.doc() = "Nanobind bridge to the Morse persistence C++ prototype.";

  nb::class_<PyFilteredComplex>(m, "FilteredComplex")
      .def(nb::init<>())
      .def("add_simplex", &PyFilteredComplex::add_simplex, nb::arg("vertices"), nb::arg("filtration"))
      .def("finalize", &PyFilteredComplex::finalize)
      .def_prop_ro("finalized", [](const PyFilteredComplex& self) { return self.finalized; })
      .def_prop_ro("size", [](const PyFilteredComplex& self) {
        require_finalized(self);
        return self.complex.size();
      })
      .def_prop_ro("num_levels", [](const PyFilteredComplex& self) {
        require_finalized(self);
        return self.complex.num_levels();
      })
      .def_prop_ro("level_values", [](const PyFilteredComplex& self) {
        require_finalized(self);
        return self.complex.level_values();
      })
      .def_prop_ro("filtration_order", [](const PyFilteredComplex& self) {
        require_finalized(self);
        return self.complex.filtration_order();
      })
      .def("vertices", [](const PyFilteredComplex& self, morseframes::SimplexId simplex) {
        require_finalized(self);
        return self.complex.vertices(simplex);
      })
      .def("dimension", [](const PyFilteredComplex& self, morseframes::SimplexId simplex) {
        require_finalized(self);
        return self.complex.dimension(simplex);
      })
      .def("level", [](const PyFilteredComplex& self, morseframes::SimplexId simplex) {
        require_finalized(self);
        return self.complex.level(simplex);
      })
      .def("filtration", [](const PyFilteredComplex& self, morseframes::SimplexId simplex) {
        require_finalized(self);
        return self.complex.filtration(simplex);
      })
      .def("boundary", [](const PyFilteredComplex& self, morseframes::SimplexId simplex) {
        require_finalized(self);
        return self.complex.boundary(simplex);
      })
      .def("coboundary", [](const PyFilteredComplex& self, morseframes::SimplexId simplex) {
        require_finalized(self);
        return self.complex.coboundary(simplex);
      })
      .def("simplices_of_level", [](const PyFilteredComplex& self, morseframes::LevelId level) {
        require_finalized(self);
        return self.complex.simplices_of_level(level);
      });

  nb::class_<PySimplexTreeBuilder>(m, "SimplexTreeBuilder")
      .def(nb::init<const std::string&>(), nb::arg("merge") = "min")
      .def("insert",
           &PySimplexTreeBuilder::insert,
           nb::arg("vertices"),
           nb::arg("filtration") = 0.0,
           nb::arg("include_faces") = true)
      .def("add_simplex",
           &PySimplexTreeBuilder::insert,
           nb::arg("vertices"),
           nb::arg("filtration") = 0.0,
           nb::arg("include_faces") = true)
      .def("insert_simplex_only",
           &PySimplexTreeBuilder::insert_simplex_only,
           nb::arg("vertices"),
           nb::arg("filtration") = 0.0)
      .def("contains",
           [](const PySimplexTreeBuilder& self, std::vector<morseframes::VertexId> vertices) {
             return self.builder.contains(std::move(vertices));
           },
           nb::arg("vertices"))
      .def("find_simplex",
           [](const PySimplexTreeBuilder& self, std::vector<morseframes::VertexId> vertices) {
             return self.builder.find_simplex(std::move(vertices));
           },
           nb::arg("vertices"))
      .def("filtration",
           [](const PySimplexTreeBuilder& self, std::vector<morseframes::VertexId> vertices) {
             return self.builder.filtration(std::move(vertices));
           },
           nb::arg("vertices"))
      .def("simplex_filtration",
           [](const PySimplexTreeBuilder& self, std::vector<morseframes::VertexId> vertices) {
             return self.builder.simplex_filtration(std::move(vertices));
           },
           nb::arg("vertices"))
      .def_prop_ro("size", [](const PySimplexTreeBuilder& self) {
        return self.builder.size();
      })
      .def("num_simplices", [](const PySimplexTreeBuilder& self) {
        return self.builder.num_simplices();
      })
      .def("num_vertices", [](const PySimplexTreeBuilder& self) {
        return self.builder.num_vertices();
      })
      .def_prop_ro("max_dimension", [](const PySimplexTreeBuilder& self) {
        return self.builder.max_dimension();
      })
      .def("simplices", [](const PySimplexTreeBuilder& self) {
        return self.builder.simplices();
      })
      .def("get_filtration", [](const PySimplexTreeBuilder& self) {
        return self.builder.get_filtration();
      })
      .def("to_filtered_complex",
           &PySimplexTreeBuilder::to_filtered_complex,
           nb::arg("finalize") = true)
      .def("finalize",
           &PySimplexTreeBuilder::finalize,
           nb::arg("clear") = true)
      .def("clear", [](PySimplexTreeBuilder& self) { self.builder.clear(); });

  nb::class_<PyMorseSequence>(m, "MorseSequence")
      .def("steps", &sequence_steps_to_python)
      .def_prop_ro("critical_simplices", [](const PyMorseSequence& self) {
        return self.sequence.critical_simplices();
      })
      .def_prop_ro("critical_index_of_simplex", [](const PyMorseSequence& self) {
        return self.sequence.critical_index_of_simplex();
      })
      .def("critical_index", [](const PyMorseSequence& self, morseframes::SimplexId simplex) {
        return self.sequence.critical_index(simplex);
      })
      .def("is_critical", [](const PyMorseSequence& self, morseframes::SimplexId simplex) {
        return self.sequence.is_critical(simplex);
      })
      .def("paired_with",
           [](const PyMorseSequence& self, const PyFilteredComplex& complex) {
             return paired_with_to_python(complex, self);
           },
           nb::arg("complex"));

  nb::class_<PyReferenceMap>(m, "ReferenceMap")
      .def_prop_ro("size", [](const PyReferenceMap& self) { return self.references.size(); })
      .def("annotation", [](const PyReferenceMap& self, morseframes::SimplexId simplex) {
        return annotation_to_python(self.references.at(simplex));
      })
      .def("annotations", &reference_map_to_python);

  nb::class_<PyMorseReferenceFrame>(m, "MorseReferenceFrame")
      .def_prop_ro("sequence", [](const PyMorseReferenceFrame& self) {
        return PyMorseSequence{self.sequence};
      })
      .def_prop_ro("references", [](const PyMorseReferenceFrame& self) {
        return PyReferenceMap{self.references};
      });

  nb::class_<PyMorseCoreferenceFrame>(m, "MorseCoreferenceFrame")
      .def_prop_ro("sequence", [](const PyMorseCoreferenceFrame& self) {
        return PyMorseSequence{self.sequence};
      })
      .def_prop_ro("coreferences", [](const PyMorseCoreferenceFrame& self) {
        return PyReferenceMap{self.coreferences};
      });

  m.def("compute_morse_sequence",
        &build_sequence,
        nb::arg("complex"),
        nb::arg("algorithm") = "saturated");
  m.def("compute_morse_sequence_and_reference_map_object",
        &build_sequence_and_reference_map,
        nb::arg("complex"),
        nb::arg("algorithm") = "saturated");
  m.def("compute_morse_sequence_and_coreference_map_object",
        &build_sequence_and_coreference_map,
        nb::arg("complex"),
        nb::arg("algorithm") = "same-level-reduction");
  m.def("compute_reference_map_object",
        &compute_reference_map_object,
        nb::arg("complex"),
        nb::arg("sequence"));
  m.def("compute_coreference_map_object",
        &compute_coreference_map_object,
        nb::arg("complex"),
        nb::arg("sequence"));
  m.def("compute_reference_map", &compute_reference_map, nb::arg("complex"), nb::arg("sequence"));
  m.def("compute_reference_map_modp",
        &compute_reference_map_modp,
        nb::arg("complex"),
        nb::arg("sequence"),
        nb::arg("modulus"));
  m.def("compute_coreference_map",
        &compute_coreference_map,
        nb::arg("complex"),
        nb::arg("sequence"));
  m.def("compute_coreference_map_modp",
        &compute_coreference_map_modp,
        nb::arg("complex"),
        nb::arg("sequence"),
        nb::arg("modulus"));
  m.def("compute_morse_persistence",
        &compute_morse_persistence,
        nb::arg("complex"),
        nb::arg("sequence"));
  m.def("compute_morse_persistence_modp",
        &compute_morse_persistence_modp,
        nb::arg("complex"),
        nb::arg("sequence"),
        nb::arg("modulus"));
  m.def("compute_morse_coreference_persistence",
        &compute_morse_coreference_persistence,
        nb::arg("complex"),
        nb::arg("sequence"));
  m.def("compute_morse_coreference_persistence_modp",
        &compute_morse_coreference_persistence_modp,
        nb::arg("complex"),
        nb::arg("sequence"),
        nb::arg("modulus"));
  m.def("reduce_morse_persistence",
        &reduce_morse_persistence,
        nb::arg("complex"),
        nb::arg("sequence"),
        nb::arg("references"));
  m.def("reduce_morse_persistence_with_metrics",
        &reduce_morse_persistence_with_metrics,
        nb::arg("complex"),
        nb::arg("sequence"),
        nb::arg("references"));
  m.def("reduce_morse_persistence_object",
        &reduce_morse_persistence_object,
        nb::arg("complex"),
        nb::arg("sequence"),
        nb::arg("references"));
  m.def("reduce_morse_coreference_persistence_object",
        &reduce_morse_coreference_persistence_object,
        nb::arg("complex"),
        nb::arg("sequence"),
        nb::arg("coreferences"));
  m.def("reduce_morse_persistence_object_with_metrics",
        &reduce_morse_persistence_object_with_metrics,
        nb::arg("complex"),
        nb::arg("sequence"),
        nb::arg("references"));
  m.def("reduce_morse_reference_frame_object",
        &reduce_morse_reference_frame_object,
        nb::arg("complex"),
        nb::arg("frame"));
  m.def("reduce_morse_reference_frame_object_with_metrics",
        &reduce_morse_reference_frame_object_with_metrics,
        nb::arg("complex"),
        nb::arg("frame"));
  m.def("reduce_morse_coreference_frame_object",
        &reduce_morse_coreference_frame_object,
        nb::arg("complex"),
        nb::arg("frame"));
  m.def("reduce_morse_coreference_frame_object_with_metrics",
        &reduce_morse_coreference_frame_object_with_metrics,
        nb::arg("complex"),
        nb::arg("frame"));
  m.def("compute_standard_persistence", &compute_standard_persistence, nb::arg("complex"));
  m.def("compute_standard_persistence_modp",
        &compute_standard_persistence_modp,
        nb::arg("complex"),
        nb::arg("modulus"));
  m.def("benchmark_morse_reference_core",
        &benchmark_morse_reference_core,
        nb::arg("complex"),
        nb::arg("algorithm") = "saturated",
        nb::arg("frame_mode") = "fused");
  m.def("profile_morse_reference_frame_core",
        &profile_morse_reference_frame_core,
        nb::arg("complex"),
        nb::arg("algorithm") = "saturated");
  m.def("benchmark_coreduction_directions_core",
        &benchmark_coreduction_directions_core,
        nb::arg("complex"));
  m.def("analyze", &analyze, nb::arg("simplices"));
}
