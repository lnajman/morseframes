#pragma once

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "morseframes/annotation.hpp"
#include "morseframes/coreference_persistence.hpp"
#include "morseframes/filtered_complex.hpp"
#include "morseframes/morse_sequence.hpp"
#include "morseframes/reference_persistence.hpp"
#include "morseframes/standard_persistence.hpp"
#include "morseframes/working_sets.hpp"

namespace morseframes {

struct AnnotationStats {
  std::size_t count = 0;
  std::size_t total_size = 0;
  std::size_t max_size = 0;
  double average_size = 0.0;
};

struct StructuralMetrics {
  std::size_t num_simplices = 0;
  std::size_t num_critical = 0;
  std::size_t num_regular_pairs = 0;
  std::size_t w_boundary_plus_size = 0;
  std::size_t w_coboundary_plus_size = 0;
  std::vector<std::size_t> critical_by_dimension;
};

struct TimingMetrics {
  double f_sequence_ms = 0.0;
  double reference_compute_ms = 0.0;
  double reference_reduce_ms = 0.0;
  double coreference_compute_ms = 0.0;
  double coreference_reduce_ms = 0.0;
  double standard_reduce_ms = 0.0;
};

struct PipelineMetrics {
  StructuralMetrics structural;
  AnnotationStats reference_annotations;
  AnnotationStats coreference_annotations;
  TimingMetrics timings;
};

struct InstrumentedPersistenceResult {
  PersistenceDiagram reference_diagram;
  PersistenceDiagram coreference_diagram;
  PersistenceDiagram standard_diagram;
  PipelineMetrics metrics;
};

inline double elapsed_ms(std::chrono::steady_clock::time_point start,
                         std::chrono::steady_clock::time_point end) {
  return std::chrono::duration<double, std::milli>(end - start).count();
}

inline AnnotationStats compute_annotation_stats(const std::vector<Annotation>& annotations) {
  AnnotationStats stats;
  stats.count = annotations.size();

  for (const auto& annotation : annotations) {
    stats.total_size += annotation.size();
    if (annotation.size() > stats.max_size) {
      stats.max_size = annotation.size();
    }
  }

  if (stats.count != 0) {
    stats.average_size = static_cast<double>(stats.total_size) / static_cast<double>(stats.count);
  }
  return stats;
}

inline AnnotationStats compute_annotation_stats_for_simplices(
    const std::vector<Annotation>& annotations,
    const std::vector<SimplexId>& simplices) {
  AnnotationStats stats;
  stats.count = simplices.size();

  for (SimplexId simplex : simplices) {
    const auto& annotation = annotations.at(simplex);
    stats.total_size += annotation.size();
    if (annotation.size() > stats.max_size) {
      stats.max_size = annotation.size();
    }
  }

  if (stats.count != 0) {
    stats.average_size = static_cast<double>(stats.total_size) / static_cast<double>(stats.count);
  }
  return stats;
}

inline StructuralMetrics compute_structural_metrics(const FilteredSimplicialComplex& complex,
                                                     const MorseSequence& sequence) {
  StructuralMetrics metrics;
  metrics.num_simplices = complex.size();
  metrics.num_critical = sequence.critical_simplices().size();

  std::uint16_t max_dimension = 0;
  for (SimplexId simplex = 0; simplex < complex.size(); ++simplex) {
    if (complex.dimension(simplex) > max_dimension) {
      max_dimension = complex.dimension(simplex);
    }
  }
  metrics.critical_by_dimension.assign(static_cast<std::size_t>(max_dimension) + 1, 0);

  for (const MorseStep& step : sequence.steps()) {
    if (step.type == MorseStepType::RegularPair) {
      ++metrics.num_regular_pairs;
    }
  }

  for (SimplexId critical : sequence.critical_simplices()) {
    ++metrics.critical_by_dimension[complex.dimension(critical)];
  }

  metrics.w_boundary_plus_size = reference_working_set(complex, sequence).size();
  metrics.w_coboundary_plus_size = coreference_working_set(complex, sequence).size();

  return metrics;
}

inline InstrumentedPersistenceResult run_instrumented_persistence(
    const FilteredSimplicialComplex& complex) {
  InstrumentedPersistenceResult result;

  auto start = std::chrono::steady_clock::now();
  auto sequence = FSequenceBuilder(complex).build_saturated();
  auto end = std::chrono::steady_clock::now();
  result.metrics.timings.f_sequence_ms = elapsed_ms(start, end);
  result.metrics.structural = compute_structural_metrics(complex, sequence);

  start = std::chrono::steady_clock::now();
  auto references = MorseReferenceComputer(complex, sequence).compute_full_references();
  end = std::chrono::steady_clock::now();
  result.metrics.timings.reference_compute_ms = elapsed_ms(start, end);
  result.metrics.reference_annotations =
      compute_annotation_stats_for_simplices(references, reference_working_set(complex, sequence));

  start = std::chrono::steady_clock::now();
  result.reference_diagram =
      MorseReferencePersistenceReducer(complex, sequence, references).compute();
  end = std::chrono::steady_clock::now();
  result.metrics.timings.reference_reduce_ms = elapsed_ms(start, end);

  start = std::chrono::steady_clock::now();
  auto coreferences = MorseCoreferenceComputer(complex, sequence).compute_full_coreferences();
  end = std::chrono::steady_clock::now();
  result.metrics.timings.coreference_compute_ms = elapsed_ms(start, end);
  result.metrics.coreference_annotations = compute_annotation_stats_for_simplices(
      coreferences, coreference_working_set(complex, sequence));

  start = std::chrono::steady_clock::now();
  result.coreference_diagram =
      MorseCoreferencePersistenceReducer(complex, sequence, coreferences).compute();
  end = std::chrono::steady_clock::now();
  result.metrics.timings.coreference_reduce_ms = elapsed_ms(start, end);

  start = std::chrono::steady_clock::now();
  result.standard_diagram = compute_standard_z2_persistence(complex);
  end = std::chrono::steady_clock::now();
  result.metrics.timings.standard_reduce_ms = elapsed_ms(start, end);

  return result;
}

}  // namespace morseframes
