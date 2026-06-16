#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "morseframes/coreference_persistence.hpp"
#include "morseframes/filtered_complex.hpp"
#include "morseframes/instrumentation.hpp"
#include "morseframes/morse_sequence.hpp"
#include "morseframes/reference_persistence.hpp"
#include "morseframes/standard_persistence.hpp"

namespace {

using morseframes::FilteredSimplicialComplex;

void add_simplex(FilteredSimplicialComplex& complex,
                 std::initializer_list<morseframes::VertexId> vertices,
                 double filtration) {
  complex.add_simplex(std::vector<morseframes::VertexId>(vertices), filtration);
}

FilteredSimplicialComplex make_example(const std::string& name) {
  FilteredSimplicialComplex complex;

  if (name == "one_vertex") {
    add_simplex(complex, {0}, 0.0);
  } else if (name == "edge_later") {
    add_simplex(complex, {0}, 0.0);
    add_simplex(complex, {1}, 0.0);
    add_simplex(complex, {0, 1}, 1.0);
  } else if (name == "edge_same") {
    add_simplex(complex, {0}, 0.0);
    add_simplex(complex, {1}, 0.0);
    add_simplex(complex, {0, 1}, 0.0);
  } else if (name == "triangle_boundary") {
    add_simplex(complex, {0}, 0.0);
    add_simplex(complex, {1}, 0.0);
    add_simplex(complex, {2}, 0.0);
    add_simplex(complex, {0, 1}, 1.0);
    add_simplex(complex, {0, 2}, 1.0);
    add_simplex(complex, {1, 2}, 1.0);
  } else if (name == "filled_triangle") {
    add_simplex(complex, {0}, 0.0);
    add_simplex(complex, {1}, 0.0);
    add_simplex(complex, {2}, 0.0);
    add_simplex(complex, {0, 1}, 1.0);
    add_simplex(complex, {0, 2}, 1.0);
    add_simplex(complex, {1, 2}, 1.0);
    add_simplex(complex, {0, 1, 2}, 2.0);
  } else {
    throw std::invalid_argument("Unknown example: " + name);
  }

  return complex;
}

void print_diagram(const morseframes::PersistenceDiagram& diagram) {
  for (const auto& pair : morseframes::off_diagonal_pairs(diagram)) {
    std::cout << "finite\t" << pair.dimension << '\t' << pair.birth_value << '\t'
              << pair.death_value << '\n';
  }
  for (const auto& interval : diagram.essential) {
    std::cout << "essential\t" << interval.dimension << '\t' << interval.birth_value << "\tinf\n";
  }
}

void print_metrics(const morseframes::PipelineMetrics& metrics) {
  const auto& structural = metrics.structural;
  std::cout << "metric\tnum_simplices\t" << structural.num_simplices << '\n';
  std::cout << "metric\tnum_critical\t" << structural.num_critical << '\n';
  std::cout << "metric\tnum_regular_pairs\t" << structural.num_regular_pairs << '\n';
  std::cout << "metric\tw_boundary_plus_size\t" << structural.w_boundary_plus_size << '\n';
  std::cout << "metric\tw_coboundary_plus_size\t" << structural.w_coboundary_plus_size << '\n';

  for (std::size_t dim = 0; dim < structural.critical_by_dimension.size(); ++dim) {
    std::cout << "metric\tcritical_dim_" << dim << '\t'
              << structural.critical_by_dimension[dim] << '\n';
  }

  std::cout << "metric\treference_annotation_total\t"
            << metrics.reference_annotations.total_size << '\n';
  std::cout << "metric\treference_annotation_max\t" << metrics.reference_annotations.max_size
            << '\n';
  std::cout << "metric\treference_annotation_average\t"
            << metrics.reference_annotations.average_size << '\n';
  std::cout << "metric\tcoreference_annotation_total\t"
            << metrics.coreference_annotations.total_size << '\n';
  std::cout << "metric\tcoreference_annotation_max\t" << metrics.coreference_annotations.max_size
            << '\n';
  std::cout << "metric\tcoreference_annotation_average\t"
            << metrics.coreference_annotations.average_size << '\n';

  std::cout << "timing_ms\tf_sequence\t" << metrics.timings.f_sequence_ms << '\n';
  std::cout << "timing_ms\treference_compute\t" << metrics.timings.reference_compute_ms << '\n';
  std::cout << "timing_ms\treference_reduce\t" << metrics.timings.reference_reduce_ms << '\n';
  std::cout << "timing_ms\tcoreference_compute\t" << metrics.timings.coreference_compute_ms
            << '\n';
  std::cout << "timing_ms\tcoreference_reduce\t" << metrics.timings.coreference_reduce_ms
            << '\n';
  std::cout << "timing_ms\tstandard_reduce\t" << metrics.timings.standard_reduce_ms << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2 && argc != 3) {
    std::cerr << "usage: morse_example_barcode <example-name> "
                 "[reference|coreference|standard|metrics]\n";
    return 2;
  }

  const std::string algorithm = argc == 3 ? argv[2] : "reference";
  auto complex = make_example(argv[1]);
  complex.finalize();

  morseframes::PersistenceDiagram diagram;
  if (algorithm == "reference") {
    auto sequence = morseframes::FSequenceBuilder(complex).build_saturated();
    diagram = morseframes::compute_morse_reference_persistence(complex, sequence);
  } else if (algorithm == "coreference") {
    auto sequence = morseframes::FSequenceBuilder(complex).build_saturated();
    diagram = morseframes::compute_morse_coreference_persistence(complex, sequence);
  } else if (algorithm == "standard") {
    diagram = morseframes::compute_standard_z2_persistence(complex);
  } else if (algorithm == "metrics") {
    const auto result = morseframes::run_instrumented_persistence(complex);
    print_metrics(result.metrics);
    return 0;
  } else {
    std::cerr << "unknown algorithm: " << algorithm << '\n';
    return 2;
  }

  print_diagram(diagram);
  return 0;
}
