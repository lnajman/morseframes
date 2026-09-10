// Fresh gradients from identical resident arrays, without prepared native data.
// Performance calls read clocks only at outer phase boundaries, with no local
// profiling. Separate diagnostic calls add internal gradient metrics.
#include <DiscreteGradient.h>
#include <OrderDisambiguation.h>
#include <Triangulation.h>

#include "morseframes/debug_checks.hpp"
#include "morseframes/lower_star_complex.hpp"
#include "morseframes/reduction_kernel_sequence.hpp"

#include <algorithm>
#include <array>
#include "../pls_profile.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
using Complex = morseframes::FilteredSimplicialComplex;
using Builder = morseframes::FSequenceBuilder<Complex>;
using RkBuilder = morseframes::ReductionKernelSequenceBuilder<Complex>;
using Sequence = morseframes::MorseSequence;
using Phases = std::map<std::string, double>;
constexpr const char* kSchema = "resident-gradient-v2";
constexpr std::array<const char*, 3> kAlgorithms{{"f_max", "reduction_kernel", "ttk"}};
constexpr std::array<std::array<int, 3>, 6> kOrders{{
    {{0, 1, 2}}, {{2, 1, 0}}, {{1, 2, 0}},
    {{0, 2, 1}}, {{2, 0, 1}}, {{1, 0, 2}}}};
constexpr const char* kPlsAlgorithm = "process_lower_stars";
// Williams design: every method occupies each position once; all 12 directed
// adjacent method pairs occur once within a four-round cycle.
constexpr std::array<std::array<int, 4>, 4> kPlsOrders{{
    {{0, 1, 3, 2}}, {{1, 2, 0, 3}}, {{2, 3, 1, 0}}, {{3, 0, 2, 1}}}};

double seconds(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double>(b - a).count();
}

// The sole common starting representation. No ordering, filtration extension,
// simplex enumeration, adjacency, native triangulation, or builder is cached.
struct Input {
  int dimension = 0;
  std::vector<double> values;
  std::vector<float> points;
  std::vector<std::vector<morseframes::VertexId>> cells;
};

Input read_input(const std::string& path) {
  std::ifstream stream(path);
  Input result;
  std::string magic;
  std::size_t vertices = 0, cells = 0;
  stream >> magic >> result.dimension >> vertices >> cells;
  if (!stream || magic != "morseframes-ttk-v1" ||
      result.dimension < 1 || result.dimension > 3 || vertices == 0 || cells == 0 ||
      vertices > std::numeric_limits<morseframes::VertexId>::max() ||
      vertices > static_cast<std::size_t>(std::numeric_limits<ttk::SimplexId>::max())) {
    throw std::runtime_error("Invalid resident benchmark input header.");
  }
  result.values.resize(vertices);
  result.points.assign(3 * vertices, 0.0F);
  for (std::size_t v = 0; v < vertices; ++v) {
    stream >> result.values[v];
    if (!stream || !std::isfinite(result.values[v])) {
      throw std::runtime_error("Invalid vertex value.");
    }
    result.points[3 * v] = static_cast<float>(v);
  }
  // This adapter represents pure 1D--3D meshes, just like the legacy driver.
  // Geometry is immaterial to the classic combinatorial gradient backend.
  std::vector<bool> used(vertices, false);
  for (std::size_t i = 0; i < cells; ++i) {
    std::vector<morseframes::VertexId> cell(result.dimension + 1);
    for (auto& v : cell) {
      stream >> v;
      if (!stream || v >= vertices) {
        throw std::runtime_error("Invalid cell vertex.");
      }
      used[v] = true;
    }
    auto sorted = cell;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
      throw std::runtime_error("Degenerate input cell.");
    }
    result.cells.push_back(std::move(cell));
  }
  if (std::find(used.begin(), used.end(), false) != used.end()) {
    throw std::runtime_error("This pure-mesh adapter does not support isolated vertices.");
  }
  return result;
}

struct Timing {
  double total = 0.0;
  Phases phases;
  Phases gradient_details;
  Phases pls_profile;
};

struct MorseRun {
  Complex complex;
  morseframes::MorseSequenceBuildMetrics metrics;
  std::unique_ptr<Builder> builder;
  std::unique_ptr<RkBuilder> rk_builder;
  std::optional<Sequence> sequence;
  Timing timing;
};

void populate_complex(const Input& input, Complex& complex) {
  morseframes::add_lower_star_cells(complex, input.values, input.cells);
  complex.finalize();
}

template <bool Diagnostic>
std::unique_ptr<MorseRun> run_morse(const Input& input, int algorithm, int workers) {
  const auto start = Clock::now();
  auto run = std::make_unique<MorseRun>();
  populate_complex(input, run->complex);
  const auto representation_stop = Clock::now();
  // Coarse RK profiling leaves local kernels uninstrumented. F-Max's existing
  // diagnostics are finer-grained and are used only in separate phase runs.
  if (algorithm == 0 || algorithm == 3) {
    run->builder = std::make_unique<Builder>(
        run->complex, Diagnostic ? &run->metrics : nullptr, false);
  } else {
    run->rk_builder = std::make_unique<RkBuilder>(
        run->complex, Diagnostic ? &run->metrics : nullptr, false);
  }
  const auto builder_stop = Clock::now();
  if (algorithm == 0) {
    run->sequence.emplace(run->builder->build_f_max());
  } else if (algorithm == 3) {
    run->sequence.emplace(workers == 1
        ? run->builder->build_process_lower_stars()
        : run->builder->build_process_lower_stars_parallel(workers));
  } else if (workers == 1) {
    run->sequence.emplace(run->rk_builder->build_flooding_reduction_kernel());
  } else {
    run->sequence.emplace(
        run->rk_builder->build_flooding_reduction_kernel_parallel(workers));
  }
  const auto stop = Clock::now(); // Gradient is now available; keep it alive.
  run->timing.total = seconds(start, stop);
  run->timing.phases = {
      {"representation_and_filtration", seconds(start, representation_stop)},
      {"builder_setup", seconds(representation_stop, builder_stop)},
      {"gradient", seconds(builder_stop, stop)}};
  if constexpr (Diagnostic) {
    const auto& m = run->metrics;
    if (algorithm == 0) {
      run->timing.gradient_details = {
          {"workspace_init", 1e-9 * m.init_nanoseconds},
          {"candidate_seeding", 1e-9 * m.candidate_seed_nanoseconds},
          {"candidate_selection", 1e-9 * m.candidate_loop_nanoseconds},
          {"emission_and_updates", 1e-9 * m.emit_nanoseconds},
          {"callbacks", 1e-9 * m.callback_nanoseconds}};
    } else if (algorithm == 3) {
      run->timing.pls_profile = pls_phase_profile(m);
      run->timing.gradient_details = {
          {"lower_star_setup", 1e-9 * m.process_lower_stars_setup_nanoseconds},
          {"local_processing", 1e-9 * m.process_lower_stars_local_wall_nanoseconds},
          {"replay", 1e-9 * m.process_lower_stars_replay_nanoseconds}};
    } else {
      run->timing.gradient_details = {
          {"workspace_and_pool", 1e-9 * m.reduction_kernel_setup_nanoseconds},
          {"level_processing", 1e-9 * m.reduction_kernel_level_wall_nanoseconds},
          {"replay", 1e-9 * m.reduction_kernel_replay_nanoseconds}};
    }
  }
  return run;
}

struct TtkRun {
  std::vector<ttk::SimplexId> offsets;
  std::vector<ttk::LongSimplexId> cells;
  ttk::Triangulation triangulation;
  ttk::dcg::DiscreteGradient gradient;
  Timing timing;
};

template <bool Diagnostic>
std::unique_ptr<TtkRun> run_ttk(const Input& input, int workers) {
  const auto start = Clock::now();
  auto run = std::make_unique<TtkRun>();
  const auto init_stop = Clock::now();
  // No precomputed ordering is supplied by the input parser. Ties use vertex
  // IDs consistently, without changing the original scalar values.
  run->offsets.resize(input.values.size());
  ttk::preconditionOrderArray(input.values.size(), input.values.data(),
                              run->offsets.data(), workers);
  const auto order_stop = Clock::now();
  run->cells.reserve(input.cells.size() * (input.dimension + 2));
  for (const auto& cell : input.cells) {
    run->cells.push_back(cell.size());
    run->cells.insert(run->cells.end(), cell.begin(), cell.end());
  }
  auto& triangulation = run->triangulation;
  triangulation.setDebugLevel(0);
  triangulation.setThreadNumber(workers);
  if (triangulation.setInputPoints(input.values.size(), input.points.data()) != 0 ||
      triangulation.setInputCells(input.cells.size(), run->cells.data()) != 0) {
    throw std::runtime_error("TTK rejected the resident mesh.");
  }
  auto& gradient = run->gradient;
  gradient.setDebugLevel(0);
  gradient.setThreadNumber(workers);
  gradient.setBackend(ttk::dcg::DiscreteGradient::BACKEND::CLASSIC_BACKEND);
  gradient.setInputScalarField(input.values.data(), 1);
  gradient.setInputOffsets(run->offsets.data());
  const auto setup_stop = Clock::now();
  gradient.preconditionTriangulation(&triangulation);
  const auto precondition_stop = Clock::now();
  if (gradient.buildGradient<double>(triangulation, true) != 0) {
    throw std::runtime_error("TTK gradient construction failed.");
  }
  const auto stop = Clock::now();
  run->timing.total = seconds(start, stop);
  run->timing.phases = {
      {"native_object_init", seconds(start, init_stop)},
      {"vertex_order", seconds(init_stop, order_stop)},
      {"representation_setup", seconds(order_stop, setup_stop)},
      {"connectivity_precondition", seconds(setup_stop, precondition_stop)},
      {"gradient", seconds(precondition_stop, stop)}};
  // Lower stars are built inside this call, not an excluded preparation.
  // Unmodified TTK does not expose separate lower-star/matching timers.
  return run;
}

std::vector<std::size_t> morse_counts(const MorseRun& run, int dimension) {
  std::vector<std::size_t> result(dimension + 1, 0);
  for (auto simplex : run.sequence->critical_simplices()) {
    ++result[run.complex.dimension(simplex)];
  }
  return result;
}

void compare_sequences(const Sequence& a, const Sequence& b) {
  if (a.steps().size() != b.steps().size()) {
    throw std::runtime_error("Gradient sequence length changed.");
  }
  for (std::size_t i = 0; i < a.steps().size(); ++i) {
    const auto& x = a.steps()[i];
    const auto& y = b.steps()[i];
    if (x.type != y.type || x.sigma != y.sigma || x.tau != y.tau || x.level != y.level) {
      throw std::runtime_error("Gradient differs from the one-worker reference.");
    }
  }
}

struct TtkSignature {
  std::vector<std::size_t> counts;
  std::vector<std::array<ttk::SimplexId, 3>> cells;
};

TtkSignature ttk_signature(const TtkRun& run, int dimension) {
  TtkSignature result;
  result.counts.assign(dimension + 1, 0);
  for (int dim = 0; dim <= dimension; ++dim) {
    const auto count = run.gradient.getNumberOfCells(dim, run.triangulation);
    for (ttk::SimplexId id = 0; id < count; ++id) {
      const bool critical = run.gradient.isCellCritical(dim, id);
      result.counts[dim] += critical;
      const ttk::dcg::Cell cell{dim, id};
      result.cells.push_back({critical ? 1 : 0,
          dim < dimension ? run.gradient.getPairedCell(cell, run.triangulation) : -1,
          dim > 0 ? run.gradient.getPairedCell(cell, run.triangulation, true) : -1});
    }
  }
  return result;
}

struct Options {
  std::string input;
  int workers = 1, repeats = 6, diagnostics = 3, warmups = 1;
  bool include_pls = false;
  int order_offset = 0;
};

Options parse_options(int argc, char** argv) {
  Options result;
  for (int i = 1; i < argc; i += 2) {
    if (i + 1 == argc) throw std::runtime_error("Missing option value.");
    const std::string key = argv[i], value = argv[i + 1];
    if (key == "--input") result.input = value;
    else if (key == "--workers") result.workers = std::stoi(value);
    else if (key == "--repeats") result.repeats = std::stoi(value);
    else if (key == "--diagnostics") result.diagnostics = std::stoi(value);
    else if (key == "--warmups") result.warmups = std::stoi(value);
    else if (key == "--include-pls" && (value == "0" || value == "1")) result.include_pls = value == "1";
    else if (key == "--order-offset") result.order_offset = std::stoi(value);
    else throw std::runtime_error("Unknown option: " + key);
  }
  if (result.input.empty() || result.workers < 1 || result.repeats < 1 ||
      result.diagnostics < 1 || result.warmups < 0 || result.order_offset < 0 ||
      result.order_offset >= (result.include_pls ? 4 : 6)) {
    throw std::runtime_error("Invalid resident benchmark options.");
  }
  return result;
}

void write_phases(const Phases& phases) {
  std::cout << '{';
  bool comma = false;
  for (const auto& [key, value] : phases) {
    if (comma) std::cout << ',';
    comma = true;
    std::cout << '"' << key << "\":" << value;
  }
  std::cout << '}';
}

void write_timing(const Timing& timing) {
  std::cout << "{\"total_seconds\":" << timing.total << ",\"phases_seconds\":";
  write_phases(timing.phases);
  std::cout << ",\"gradient_details_seconds\":";
  write_phases(timing.gradient_details);
  if (!timing.pls_profile.empty()) {
    std::cout << ",\"pls_profile_seconds\":";
    write_phases(timing.pls_profile);
  }
  std::cout << '}';
}

} // namespace

int main(int argc, char** argv) {
  try {
    const auto options = parse_options(argc, argv);
    std::vector<std::string> algorithms(kAlgorithms.begin(), kAlgorithms.end());
    std::vector<std::vector<int>> orders;
    if (options.include_pls) {
      algorithms.emplace_back(kPlsAlgorithm);
      for (const auto& order : kPlsOrders) orders.emplace_back(order.begin(), order.end());
    } else {
      for (const auto& order : kOrders) orders.emplace_back(order.begin(), order.end());
    }
    const auto order_at = [&](int round) -> const std::vector<int>& {
      return orders[(round + options.order_offset) % orders.size()];
    };
    const auto loading_start = Clock::now();
    const auto input = read_input(options.input);
    const auto loading_stop = Clock::now();
    // References and all comparisons are outside every recorded interval.
    auto f_reference = run_morse<false>(input, 0, 1);
    auto rk_reference = run_morse<false>(input, 1, 1);
    morseframes::validate_morse_sequence(f_reference->complex, *f_reference->sequence);
    morseframes::validate_morse_sequence(rk_reference->complex, *rk_reference->sequence);
    auto ttk_reference_run = run_ttk<false>(input, 1);
    const auto ttk_reference = ttk_signature(*ttk_reference_run, input.dimension);
    const auto simplex_count = f_reference->complex.size();
    std::int64_t euler = 0;
    for (morseframes::SimplexId id = 0; id < simplex_count; ++id) {
      euler += f_reference->complex.dimension(id) % 2 ? -1 : 1;
    }
    if (ttk_reference.cells.size() != simplex_count) {
      throw std::runtime_error("TTK and MorseFrames simplex counts differ.");
    }
    std::vector<std::vector<std::size_t>> counts{
        morse_counts(*f_reference, input.dimension),
        morse_counts(*rk_reference, input.dimension), ttk_reference.counts};
    // Only reference outputs need to remain resident, not prepared native data.
    const auto f_sequence = std::move(*f_reference->sequence);
    const auto rk_sequence = std::move(*rk_reference->sequence);
    std::optional<Sequence> pls_sequence;
    if (options.include_pls) {
      auto pls_reference = run_morse<false>(input, 3, 1);
      morseframes::validate_morse_sequence(pls_reference->complex, *pls_reference->sequence);
      counts.push_back(morse_counts(*pls_reference, input.dimension));
      pls_sequence.emplace(std::move(*pls_reference->sequence));
    }
    f_reference.reset();
    rk_reference.reset();
    ttk_reference_run.reset();
    for (const auto& critical_counts : counts) {
      std::int64_t critical_euler = 0;
      for (std::size_t dim = 0; dim < critical_counts.size(); ++dim) {
        critical_euler += dim % 2 ? -static_cast<std::int64_t>(critical_counts[dim])
                                 : static_cast<std::int64_t>(critical_counts[dim]);
      }
      if (critical_euler != euler) throw std::runtime_error("Critical-cell Euler characteristic differs.");
    }

    const auto measure = [&](int algorithm, bool diagnostic) {
      if (algorithm == 2) {
        auto run = diagnostic ? run_ttk<true>(input, options.workers)
                              : run_ttk<false>(input, options.workers);
        const auto actual = ttk_signature(*run, input.dimension);
        if (actual.cells != ttk_reference.cells || actual.counts != ttk_reference.counts) {
          throw std::runtime_error("TTK pairings differ from its one-worker reference.");
        }
        return run->timing;
      }
      auto run = diagnostic ? run_morse<true>(input, algorithm, options.workers)
                            : run_morse<false>(input, algorithm, options.workers);
      compare_sequences(algorithm == 0 ? f_sequence : algorithm == 1 ? rk_sequence : *pls_sequence,
                        *run->sequence);
      return run->timing;
    };
    for (int i = 0; i < options.warmups; ++i) {
      for (int algorithm : order_at(i)) (void)measure(algorithm, false);
    }
    std::vector<std::vector<Timing>> performance(algorithms.size()), diagnostics(algorithms.size());
    for (int i = 0; i < options.repeats; ++i) {
      for (int algorithm : order_at(i)) {
        performance[algorithm].push_back(measure(algorithm, false));
      }
    }
    for (int i = 0; i < options.warmups; ++i) {
      for (int algorithm : order_at(i)) (void)measure(algorithm, true);
    }
    for (int i = 0; i < options.diagnostics; ++i) {
      for (int algorithm : order_at(i)) {
        diagnostics[algorithm].push_back(measure(algorithm, true));
      }
    }

    std::cout << std::setprecision(17) << "{\"schema\":\""
              << (options.include_pls ? "resident-gradient-v3" : kSchema)
              << "\",\"ttk_revision\":\"" << MORSEFRAMES_TTK_REVISION
              << "\",\"dimension\":" << input.dimension
              << ",\"vertices\":" << input.values.size()
              << ",\"simplices\":" << simplex_count
              << ",\"workers\":" << options.workers
              << ",\"input_loading_seconds\":" << seconds(loading_start, loading_stop)
              << ",\"exact_reference_checks\":true,\"critical_counts_match\":"
              << (std::all_of(counts.begin(), counts.end(), [&](const auto& c) { return c == counts[0]; })
                      ? "true" : "false");
    if (options.include_pls) {
      std::cout << ",\"euler_characteristic\":" << euler << ",\"performance_orders\":[";
      for (int round = 0; round < options.repeats; ++round) {
        if (round) std::cout << ',';
        std::cout << '[';
        const auto& order = order_at(round);
        for (std::size_t position = 0; position < order.size(); ++position) {
          if (position) std::cout << ',';
          std::cout << '"' << algorithms[order[position]] << '"';
        }
        std::cout << ']';
      }
      std::cout << ']';
    }
    std::cout << ",\"algorithms\":{";
    for (std::size_t algorithm = 0; algorithm < algorithms.size(); ++algorithm) {
      if (algorithm != 0) std::cout << ',';
      std::cout << '"' << algorithms[algorithm] << "\":{\"critical_counts\":[";
      for (std::size_t i = 0; i < counts[algorithm].size(); ++i) {
        if (i != 0) std::cout << ',';
        std::cout << counts[algorithm][i];
      }
      std::cout << "],\"performance_seconds\":[";
      for (std::size_t i = 0; i < performance[algorithm].size(); ++i) {
        if (i != 0) std::cout << ',';
        std::cout << performance[algorithm][i].total;
      }
      std::cout << "],\"performance_phases_seconds\":[";
      for (std::size_t i = 0; i < performance[algorithm].size(); ++i) {
        if (i != 0) std::cout << ',';
        write_phases(performance[algorithm][i].phases);
      }
      std::cout << "],\"diagnostics\":[";
      for (std::size_t i = 0; i < diagnostics[algorithm].size(); ++i) {
        if (i != 0) std::cout << ',';
        write_timing(diagnostics[algorithm][i]);
      }
      std::cout << "]}";
    }
    std::cout << "}}\n";
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
