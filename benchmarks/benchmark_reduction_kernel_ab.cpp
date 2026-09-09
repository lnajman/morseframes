// Interactive worker for tools/benchmark_reduction_kernel_ab.py.
// Topology, validation and protocol I/O are outside the timed interval.
#include "morseframes/debug_checks.hpp"
#include "morseframes/filtered_complex.hpp"
#include "morseframes/morse_sequence.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Complex = morseframes::FilteredSimplicialComplex;
using Sequence = morseframes::MorseSequence;

Complex read_complex(const std::string& path) {
  std::ifstream input(path);
  std::string magic;
  std::size_t dimension = 0, vertices = 0, facets = 0;
  input >> magic >> dimension >> vertices >> facets;
  if (!input || magic != "morseframes-ttk-v1" || dimension < 1 ||
      dimension > 3 || vertices == 0 || facets == 0) {
    throw std::runtime_error("Invalid simplicial benchmark input.");
  }
  std::vector<double> values(vertices);
  for (auto& value : values) {
    input >> value;
  }
  Complex complex;
  std::vector<morseframes::VertexId> facet(dimension + 1);
  for (std::size_t index = 0; index < facets; ++index) {
    for (auto& vertex : facet) {
      input >> vertex;
      if (vertex >= vertices) {
        throw std::runtime_error("Invalid vertex identifier.");
      }
    }
    for (std::size_t mask = 1; mask < (std::size_t{1} << facet.size());
         ++mask) {
      std::vector<morseframes::VertexId> simplex;
      double weight = -std::numeric_limits<double>::infinity();
      for (std::size_t local = 0; local < facet.size(); ++local) {
        if ((mask & (std::size_t{1} << local)) != 0) {
          simplex.push_back(facet[local]);
          weight = std::max(weight, values[facet[local]]);
        }
      }
      complex.add_simplex(std::move(simplex), weight);
    }
  }
  if (!input) {
    throw std::runtime_error("Truncated input complex.");
  }
  complex.finalize();
  return complex;
}

Sequence build(const Complex& complex, std::size_t workers) {
  // Include a fresh builder, workspace, task pool and event replay in timing.
  morseframes::FSequenceBuilder<Complex> builder(complex);
  return workers == 1 ? builder.build_flooding_reduction_kernel()
                      : builder.build_flooding_reduction_kernel_parallel(workers);
}

void require_same_sequence(const Sequence& expected, const Sequence& actual) {
  if (expected.steps().size() != actual.steps().size()) {
    throw std::runtime_error("Sequence length changed.");
  }
  for (std::size_t i = 0; i < expected.steps().size(); ++i) {
    const auto& a = expected.steps()[i];
    const auto& b = actual.steps()[i];
    if (a.type != b.type || a.sigma != b.sigma || a.tau != b.tau ||
        a.level != b.level) {
      throw std::runtime_error("Sequence differs at step " + std::to_string(i));
    }
  }
}
}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 3) {
      throw std::runtime_error("Usage: worker INPUT SEQUENCE_DUMP");
    }
    const auto complex = read_complex(argv[1]);
    const auto reference = build(complex, 1);
    morseframes::validate_morse_sequence(complex, reference);
    std::vector<std::size_t> critical(4, 0);
    std::ofstream dump(argv[2]);
    for (const auto& step : reference.steps()) {
      dump << static_cast<int>(step.type) << ' ' << step.sigma << ' '
           << step.tau << ' ' << step.level << '\n';
      if (step.type == morseframes::MorseStepType::Critical) {
        ++critical[complex.dimension(step.sigma)];
      }
    }
    dump.close();
    if (!dump) {
      throw std::runtime_error("Failed to export the exact sequence.");
    }
    std::cout << "{\"num_simplices\":" << complex.size()
              << ",\"num_levels\":" << complex.num_levels()
              << ",\"critical_counts\":[" << critical[0] << ',' << critical[1]
              << ',' << critical[2] << ',' << critical[3] << "]}" << std::endl;

    std::string command;
    while (std::cin >> command) {
      if (command == "quit") {
        return 0;
      }
      std::size_t workers = 0, repeats = 0;
      std::cin >> workers >> repeats;
      if (!std::cin || command != "run" || workers == 0 || repeats == 0) {
        throw std::runtime_error("Invalid worker command.");
      }
      std::vector<double> samples;
      samples.reserve(repeats);
      for (std::size_t repeat = 0; repeat < repeats; ++repeat) {
        const auto start = std::chrono::steady_clock::now();
        const auto sequence = build(complex, workers);
        const auto stop = std::chrono::steady_clock::now();
        samples.push_back(std::chrono::duration<double>(stop - start).count());
        // Check every timed result; do not rely on counts or hashes alone.
        require_same_sequence(reference, sequence);
      }
      std::cout << '[' << std::setprecision(17);
      for (std::size_t i = 0; i < samples.size(); ++i) {
        if (i != 0) {
          std::cout << ',';
        }
        std::cout << samples[i];
      }
      std::cout << ']' << std::endl;
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
