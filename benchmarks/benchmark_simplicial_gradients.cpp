// Prepared-complex gradients in arbitrary tested dimensions. All function-
// specific PLS/RK setup and fresh builders stay inside algorithm timing.
#include "morseframes/debug_checks.hpp"
#include "morseframes/lower_star_complex.hpp"
#include "morseframes/reduction_kernel_sequence.hpp"
#include <array>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sys/resource.h>

namespace {
using Complex = morseframes::FilteredSimplicialComplex;
using Sequence = morseframes::MorseSequence;
using FBuilder = morseframes::FSequenceBuilder<Complex>;
using RBuilder = morseframes::ReductionKernelSequenceBuilder<Complex>;
using Clock = std::chrono::steady_clock;
constexpr std::array<const char*, 3> names{{"f_max", "process_lower_stars", "reduction_kernel"}};
constexpr std::array<std::array<int, 3>, 6> orders{{
    {{0, 1, 2}}, {{2, 1, 0}}, {{1, 2, 0}}, {{0, 2, 1}}, {{2, 0, 1}}, {{1, 0, 2}}}};
double seconds(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double>(b - a).count();
}
std::uint64_t peak_bytes() {
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage)) throw std::runtime_error("getrusage failed");
#ifdef __APPLE__
  return usage.ru_maxrss;
#else
  return std::uint64_t(usage.ru_maxrss) * 1024;
#endif
}
struct Input {
  unsigned dimension = 0;
  std::vector<double> values;
  std::vector<std::vector<morseframes::VertexId>> cells;
};
Input read(const char* path) {
  std::ifstream stream(path);
  Input input;
  std::string magic;
  std::size_t vertices = 0, cells = 0;
  stream >> magic >> input.dimension >> vertices >> cells;
  if (!stream || magic != "morseframes-ttk-v1" || input.dimension < 1 ||
      input.dimension > 15 || !vertices || !cells || vertices >= morseframes::kInvalidSimplex)
    throw std::runtime_error("Invalid simplicial input header");
  input.values.resize(vertices);
  for (auto& value : input.values) {
    stream >> value;
    if (!stream || !std::isfinite(value)) throw std::runtime_error("Invalid scalar value");
  }
  auto sorted = input.values;
  std::sort(sorted.begin(), sorted.end());
  if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
    throw std::runtime_error("This comparison requires injective vertex values");
  std::vector<bool> used(vertices, false);
  input.cells.assign(cells, std::vector<morseframes::VertexId>(input.dimension + 1));
  for (auto& cell : input.cells) {
    for (auto& vertex : cell) {
      stream >> vertex;
      if (!stream || vertex >= vertices) throw std::runtime_error("Invalid cell vertex");
      used[vertex] = true;
    }
    auto canonical = cell;
    std::sort(canonical.begin(), canonical.end());
    if (std::adjacent_find(canonical.begin(), canonical.end()) != canonical.end())
      throw std::runtime_error("Repeated cell vertex");
  }
  std::string trailing;
  if (stream >> trailing) throw std::runtime_error("Unexpected trailing input");
  if (std::find(used.begin(), used.end(), false) != used.end())
    throw std::runtime_error("Pure-mesh input has an unused vertex");
  return input;
}
struct Run {
  morseframes::MorseSequenceBuildMetrics metrics;
  std::unique_ptr<FBuilder> f;
  std::unique_ptr<RBuilder> rk;
  std::optional<Sequence> sequence;
  double builder_seconds = 0, kernel_seconds = 0, algorithm_seconds = 0;
};
std::unique_ptr<Run> run(const Complex& complex, int algorithm, std::size_t workers,
                         bool diagnostic = false) {
  const auto start = Clock::now();
  auto result = std::make_unique<Run>();
  auto* metrics = diagnostic ? &result->metrics : nullptr;
  if (algorithm == 2) result->rk = std::make_unique<RBuilder>(complex, metrics, false);
  else result->f = std::make_unique<FBuilder>(complex, metrics, false);
  const auto ready = Clock::now();
  if (algorithm == 0) result->sequence.emplace(result->f->build_f_max());
  else if (algorithm == 1) result->sequence.emplace(workers == 1
      ? result->f->build_process_lower_stars() : result->f->build_process_lower_stars_parallel(workers));
  else result->sequence.emplace(workers == 1
      ? result->rk->build_flooding_reduction_kernel()
      : result->rk->build_flooding_reduction_kernel_parallel(workers));
  const auto stop = Clock::now();
  result->builder_seconds = seconds(start, ready);
  result->kernel_seconds = seconds(ready, stop);
  result->algorithm_seconds = seconds(start, stop);
  return result; // Builders and output remain alive beyond the timed interval.
}
struct Hash {
  std::uint64_t value = 14695981039346656037ull;
  std::ostream* dump = nullptr;
  void add(std::uint64_t x) {
    if (dump) *dump << x << ' ';
    for (unsigned i = 0; i < 8; ++i) { value ^= (x >> (i * 8)) & 255; value *= 1099511628211ull; }
  }
  template <class Range> void list(const Range& values) {
    add(values.size()); for (auto value : values) add(value);
  }
  void real(double x) { std::uint64_t bits; std::memcpy(&bits, &x, sizeof bits); add(bits); }
};
std::uint64_t fingerprint(const Sequence& sequence, std::ostream* dump = nullptr) {
  Hash hash; hash.dump = dump; hash.add(sequence.steps().size());
  for (const auto& s : sequence.steps()) {
    hash.add(static_cast<unsigned>(s.type)); hash.add(s.sigma); hash.add(s.tau); hash.add(s.level);
  }
  return hash.value;
}
template <class T> void array(const std::vector<T>& values) {
  std::cout << '[';
  for (std::size_t i = 0; i < values.size(); ++i) { if (i) std::cout << ','; std::cout << values[i]; }
  std::cout << ']';
}
void timing(const Run& r) {
  std::cout << "{\"builder_seconds\":" << r.builder_seconds
            << ",\"kernel_seconds\":" << r.kernel_seconds
            << ",\"algorithm_seconds\":" << r.algorithm_seconds << '}';
}
} // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 3 && argc != 5) throw std::runtime_error("Usage: worker INPUT DUMP|--memory ALGORITHM WORKERS");
    std::cout << std::setprecision(17);
    const auto load_start = Clock::now();
    const auto input = read(argv[1]);
    const auto load_stop = Clock::now();
    Complex complex;
    morseframes::add_lower_star_cells(complex, input.values, input.cells);
    complex.finalize();
    const auto construct_stop = Clock::now();
    if (std::string(argv[2]) == "--memory") {
      if (argc != 5) throw std::runtime_error("Memory mode needs algorithm and workers");
      const int a = std::stoi(argv[3]); const auto workers = std::stoul(argv[4]);
      if (a < 0 || a > 2 || !workers) throw std::runtime_error("Invalid memory mode");
      const auto before = peak_bytes();
      const auto result = run(complex, a, workers);
      const auto after = peak_bytes(); // Before reference/validation allocations.
      std::cout << "{\"complex_peak_bytes\":" << before << ",\"gradient_peak_bytes\":" << after
                << ",\"simplices\":" << complex.size() << ",\"steps\":" << result->sequence->steps().size() << "}\n";
      return 0;
    }
    if (argc != 3) throw std::runtime_error("Unexpected arguments");
    std::ofstream dump(argv[2]);
    if (!dump) throw std::runtime_error("Cannot write reference dump");
    Hash topology; topology.dump = &dump; topology.add(complex.size());
    std::vector<std::size_t> stars(input.values.size(), 0);
    std::size_t incidences = 0;
    std::int64_t euler = 0;
    for (morseframes::SimplexId id = 0; id < complex.size(); ++id) {
      topology.list(complex.vertices(id)); topology.add(complex.dimension(id));
      topology.real(complex.filtration(id)); topology.add(complex.level(id));
      topology.list(complex.boundary(id)); topology.list(complex.coboundary(id));
      incidences += complex.boundary(id).size();
      euler += complex.dimension(id) % 2 ? -1 : 1;
      const auto& v = complex.vertices(id);
      ++stars[*std::max_element(v.begin(), v.end(), [&](auto a, auto b) { return input.values[a] < input.values[b]; })];
    }
    topology.list(complex.filtration_order());
    for (morseframes::LevelId l = 0; l < complex.num_levels(); ++l) {
      topology.real(complex.level_values()[l]); topology.list(complex.simplices_of_level(l));
    }
    std::array<std::uint64_t, 3> references{};
    std::cout << "{\"identity\":{\"dimension\":" << input.dimension << ",\"vertices\":" << input.values.size()
              << ",\"cells\":" << input.cells.size() << ",\"simplices\":" << complex.size()
              << ",\"incidences\":" << incidences << ",\"euler\":" << euler
              << ",\"complex_fingerprint\":\"" << topology.value << "\",\"lower_star_sizes\":";
    array(stars);
    std::cout << ",\"algorithms\":{";
    for (int a = 0; a < 3; ++a) {
      const auto result = run(complex, a, 1);
      morseframes::validate_morse_sequence(complex, *result->sequence);
      references[a] = fingerprint(*result->sequence, &dump);
      std::vector<std::size_t> counts(input.dimension + 1, 0);
      for (auto id : result->sequence->critical_simplices()) ++counts[complex.dimension(id)];
      if (a) std::cout << ',';
      std::cout << '"' << names[a] << "\":{\"fingerprint\":\"" << references[a] << "\",\"critical_counts\":";
      array(counts); std::cout << ",\"steps\":" << result->sequence->steps().size() << '}';
    }
    dump.close(); if (!dump) throw std::runtime_error("Reference dump failed");
    std::cout << "}},\"loading_seconds\":" << seconds(load_start, load_stop)
              << ",\"construction_seconds\":" << seconds(load_stop, construct_stop) << '}' << std::endl;
    std::string command;
    std::size_t round = 0;
    while (std::cin >> command && command != "quit") {
      std::size_t workers = 0, repeats = 0;
      if (command == "run") {
        std::cin >> workers >> repeats;
        if (!std::cin || !workers || !repeats) throw std::runtime_error("Invalid run command");
        std::cout << '[';
        for (std::size_t i = 0; i < repeats; ++i, ++round) {
          if (i) std::cout << ',';
          std::cout << '{'; bool first = true;
          for (int a : orders[round % orders.size()]) {
            const auto r = run(complex, a, workers);
            if (fingerprint(*r->sequence) != references[a]) throw std::runtime_error("Gradient differs from sequential reference");
            if (!first) std::cout << ','; first = false;
            std::cout << '"' << names[a] << "\":"; timing(*r);
          }
          std::cout << '}';
        }
        std::cout << ']' << std::endl;
      } else if (command == "profile") {
        std::cin >> workers;
        if (!std::cin || !workers) throw std::runtime_error("Invalid profile command");
        const auto r = run(complex, 1, workers, true);
        if (fingerprint(*r->sequence) != references[1]) throw std::runtime_error("Profile differs");
        const auto& m = r->metrics;
        std::cout << "{\"builder_seconds\":" << r->builder_seconds
                  << ",\"algorithm_seconds\":" << r->algorithm_seconds;
#define PLS_TIME(name) std::cout << ",\"" #name "_seconds\":" << 1e-9 * m.process_lower_stars_##name##_nanoseconds
        PLS_TIME(setup); PLS_TIME(local_wall); PLS_TIME(replay);
#undef PLS_TIME
        std::cout << ",\"stars\":" << m.process_lower_stars_count
                  << ",\"max_star_size\":" << m.process_lower_stars_max_star_size
                  << ",\"executor_workers\":" << m.process_lower_stars_executor_workers << '}' << std::endl;
      } else throw std::runtime_error("Unknown worker command");
    }
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n'; return 1;
  }
}
