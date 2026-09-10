// Resident-input construction A/B worker. Diagnostics and validation are not
// performance samples. Output destruction is outside both reported timers.
#include "morseframes/filtered_complex.hpp"
#if __has_include("morseframes/lower_star_complex.hpp")
#include "morseframes/lower_star_complex.hpp"
#define MORSEFRAMES_BULK_LOWER_STAR 1
#endif
#include "morseframes/morse_sequence.hpp"
#include "morseframes/reduction_kernel_sequence.hpp"
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <sys/resource.h>

namespace {
using Complex = morseframes::FilteredSimplicialComplex;
using Clock = std::chrono::steady_clock;
using Sequence = morseframes::MorseSequence;
double elapsed(Clock::time_point start, Clock::time_point stop) {
  return std::chrono::duration<double>(stop - start).count();
}
std::uint64_t peak_bytes() {
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0) throw std::runtime_error("getrusage failed");
#ifdef __APPLE__
  return usage.ru_maxrss;
#else
  return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024;
#endif
}
struct Input {
  std::vector<double> values;
  std::vector<std::vector<morseframes::VertexId>> cells;
  std::size_t attempts = 0;
};
Input read_input(const char* path) {
  std::ifstream stream(path);
  Input input;
  std::string magic;
  std::size_t dimension = 0, vertices = 0, cells = 0;
  stream >> magic >> dimension >> vertices >> cells;
  if (!stream || magic != "morseframes-ttk-v1" || dimension < 1 ||
      dimension > 15 || vertices == 0 || cells == 0)
    throw std::runtime_error("Invalid construction input");
  input.values.resize(vertices);
  for (auto& value : input.values) {
    stream >> value;
    if (!std::isfinite(value)) throw std::runtime_error("Non-finite value");
  }
  input.cells.assign(cells, std::vector<morseframes::VertexId>(dimension + 1));
  for (auto& cell : input.cells) for (auto& vertex : cell) {
    stream >> vertex;
    if (vertex >= vertices) throw std::runtime_error("Invalid vertex");
  }
  if (!stream) throw std::runtime_error("Truncated input");
  input.attempts = cells * ((std::size_t{1} << (dimension + 1)) - 1);
  return input;
}
template <bool Diagnostic>
double populate(const Input& input, Complex& complex) {
#ifdef MORSEFRAMES_BULK_LOWER_STAR
  static_assert(!Diagnostic, "Bulk diagnostics use named phase metrics.");
  morseframes::add_lower_star_cells(complex, input.values, input.cells);
  return 0;
#else
  double insertion = 0;
  for (const auto& cell : input.cells) {
    for (std::size_t mask = 1; mask < (std::size_t{1} << cell.size()); ++mask) {
      std::vector<morseframes::VertexId> simplex;
      double value = -std::numeric_limits<double>::infinity();
      for (std::size_t i = 0; i < cell.size(); ++i) if (mask & (std::size_t{1} << i)) {
        simplex.push_back(cell[i]);
        value = std::max(value, input.values[cell[i]]);
      }
      Clock::time_point start;
      if constexpr (Diagnostic) start = Clock::now();
      complex.add_simplex(std::move(simplex), value);
      if constexpr (Diagnostic) insertion += elapsed(start, Clock::now());
    }
  }
  return insertion;
#endif
}
struct Fingerprint {
  std::uint64_t value = 14695981039346656037ull;
  std::ostream* dump = nullptr;
  void add(std::uint64_t x) {
    if (dump) *dump << x << ' ';
    for (unsigned i = 0; i < 8; ++i) { value ^= (x >> (i * 8)) & 255; value *= 1099511628211ull; }
  }
  template <class T> void list(const std::vector<T>& xs) {
    add(xs.size());
    for (auto x : xs) add(x);
  }
  void real(double x) { std::uint64_t bits; std::memcpy(&bits, &x, sizeof bits); add(bits); }
};
std::uint64_t fingerprint(const Complex& complex, std::ostream* dump = nullptr) {
  Fingerprint hash; hash.dump = dump;
  hash.add(complex.size());
  for (morseframes::SimplexId id = 0; id < complex.size(); ++id) {
    hash.list(complex.vertices(id)); hash.add(complex.dimension(id));
    hash.add(complex.level(id)); hash.real(complex.filtration(id));
    hash.list(complex.boundary(id)); hash.list(complex.coboundary(id));
    if (complex.find_simplex(complex.vertices(id)) != id)
      throw std::runtime_error("Lookup does not preserve simplex IDs");
  }
  hash.list(complex.filtration_order());
  hash.add(complex.num_levels());
  for (std::size_t level = 0; level < complex.num_levels(); ++level) {
    hash.real(complex.level_values()[level]); hash.list(complex.simplices_of_level(level));
  }
  return hash.value;
}
std::uint64_t fingerprint(const Sequence& sequence, std::ostream* dump = nullptr) {
  Fingerprint hash; hash.dump = dump; hash.add(sequence.steps().size());
  for (const auto& s : sequence.steps()) {
    hash.add(static_cast<unsigned>(s.type)); hash.add(s.sigma); hash.add(s.tau); hash.add(s.level);
  }
  return hash.value;
}
Sequence gradient(const Complex& complex, bool rk, std::size_t workers) {
  if (!rk) return morseframes::FSequenceBuilder<Complex>(complex).build_f_max();
  morseframes::ReductionKernelSequenceBuilder<Complex> builder(complex, nullptr, false);
  return workers == 1 ? builder.build_flooding_reduction_kernel()
                      : builder.build_flooding_reduction_kernel_parallel(workers);
}
void profile(const Input& input) {
  Complex complex;
  auto start = Clock::now();
#ifdef MORSEFRAMES_BULK_LOWER_STAR
  morseframes::LowerStarConstructionMetrics bulk;
  morseframes::add_lower_star_cells_with_metrics(complex, input.values, input.cells, bulk);
#else
  double insertion = populate<true>(input, complex);
#endif
  double adapter = elapsed(start, Clock::now());
  morseframes::ComplexConstructionMetrics m;
  complex.finalize_with_metrics(m);
  std::cout << "{\"adapter_inclusive_seconds\":" << adapter;
#ifdef MORSEFRAMES_BULK_LOWER_STAR
  std::cout << ",\"bulk_validation_seconds\":" << bulk.validation_seconds
            << ",\"bulk_enumeration_seconds\":" << bulk.enumeration_seconds
            << ",\"bulk_sort_and_dedup_seconds\":" << bulk.sort_and_dedup_seconds
            << ",\"bulk_insertion_seconds\":" << bulk.insertion_seconds
            << ",\"generated_faces\":" << bulk.generated_faces
            << ",\"unique_faces_submitted\":" << bulk.unique_faces_submitted;
#else
  std::cout
            << ",\"canonicalization_and_dedup_seconds\":" << insertion
            << ",\"enumeration_and_clock_overhead_seconds\":" << adapter - insertion;
#endif
  std::cout
            << ",\"reset_seconds\":" << m.reset_seconds
            << ",\"index_and_simplices_seconds\":" << m.index_and_simplices_seconds
            << ",\"levels_seconds\":" << m.levels_seconds
            << ",\"boundaries_seconds\":" << m.boundaries_seconds
            << ",\"coboundaries_seconds\":" << m.coboundaries_seconds
            << ",\"orders_and_buckets_seconds\":" << m.orders_and_buckets_seconds
            << ",\"simplices\":" << complex.size() << ",\"insertion_attempts\":" << input.attempts
            << ",\"complex_fingerprint\":\"" << fingerprint(complex) << "\"}" << std::endl;
}
} // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 3) throw std::runtime_error("Usage: worker INPUT DUMP|--memory|--profile");
    std::cout << std::setprecision(17);
    const Input input = read_input(argv[1]); // All input arrays resident before timing.
    if (std::string(argv[2]) == "--profile") { profile(input); return 0; }
    if (std::string(argv[2]) == "--memory") {
      const auto before = peak_bytes();
      Complex complex; populate<false>(input, complex); complex.finalize();
      const auto after = peak_bytes(); // Before gradients, dumps or validation allocate.
      std::cout << "{\"input_peak_bytes\":" << before << ",\"construction_peak_bytes\":" << after
                << ",\"simplices\":" << complex.size() << "}\n";
      return 0;
    }
    std::uint64_t complex_hash = 0, sequence_hash[2]{};
    std::size_t simplex_count = 0;
    std::vector<std::size_t> critical[2];
    {
      Complex complex; populate<false>(input, complex); complex.finalize();
      std::ofstream dump(argv[2]);
      if (!dump) throw std::runtime_error("Cannot open exact dump");
      complex_hash = fingerprint(complex, &dump); simplex_count = complex.size();
      for (int rk = 0; rk < 2; ++rk) {
        const auto sequence = gradient(complex, rk, 1);
        sequence_hash[rk] = fingerprint(sequence, &dump);
        critical[rk].assign(input.cells.front().size(), 0);
        for (const auto& step : sequence.steps()) if (step.type == morseframes::MorseStepType::Critical)
          ++critical[rk][complex.dimension(step.sigma)];
      }
      dump.close(); if (!dump) throw std::runtime_error("Cannot write exact dump");
    }
    std::cout << "{\"simplices\":" << simplex_count << ",\"insertion_attempts\":" << input.attempts
              << ",\"complex_fingerprint\":\"" << complex_hash << "\",\"fmax_fingerprint\":\""
              << sequence_hash[0] << "\",\"rk_fingerprint\":\"" << sequence_hash[1] << "\"";
    for (int rk = 0; rk < 2; ++rk) {
      std::cout << (rk ? ",\"rk_critical_counts\":[" : ",\"fmax_critical_counts\":[");
      for (std::size_t i = 0; i < critical[rk].size(); ++i) {
        if (i) std::cout << ',';
        std::cout << critical[rk][i];
      }
      std::cout << ']';
    }
    std::cout << '}' << std::endl;
    std::string command;
    std::size_t iteration = 0;
    while (std::cin >> command && command != "quit") {
      if (command == "profile") { profile(input); continue; }
      std::size_t workers = 0, repeats = 0;
      if (command != "run" || !(std::cin >> workers >> repeats) || !workers || !repeats)
        throw std::runtime_error("Invalid command");
      std::cout << '[';
      for (std::size_t i = 0; i < repeats; ++i, ++iteration) {
        auto start = Clock::now();
        auto complex = std::make_unique<Complex>(); populate<false>(input, *complex);
        auto inserted = Clock::now(); complex->finalize(); auto ready = Clock::now();
        const double construction = elapsed(start, ready);
        const double insertion = elapsed(start, inserted);
        double times[2]{};
        for (int a = 0; a < 2; ++a) {
          const int rk = (a + iteration) % 2;
          start = Clock::now(); auto sequence = gradient(*complex, rk, workers);
          times[rk] = elapsed(start, Clock::now());
          if (fingerprint(sequence) != sequence_hash[rk]) throw std::runtime_error("Gradient changed");
        }
        if (fingerprint(*complex) != complex_hash) throw std::runtime_error("Complex changed");
        if (i) std::cout << ',';
        std::cout << "{\"construction_seconds\":" << construction
                  << ",\"enumeration_and_insertion_seconds\":" << insertion
                  << ",\"finalize_seconds\":" << elapsed(inserted, ready)
                  << ",\"fmax_seconds\":" << times[0] << ",\"rk_seconds\":" << times[1] << '}';
      }
      std::cout << ']' << std::endl;
    }
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
