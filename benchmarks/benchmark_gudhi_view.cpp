#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include "morseframes/simplex_tree_morse.hpp"

namespace {

using GudhiSimplexTree = Gudhi::Simplex_tree<>;
using GudhiFieldZp = Gudhi::persistent_cohomology::Field_Zp;
using GudhiPersistentCohomology =
    Gudhi::persistent_cohomology::Persistent_cohomology<GudhiSimplexTree, GudhiFieldZp>;
using Clock = std::chrono::steady_clock;

template <class T>
struct TimedValue {
  T value;
  double seconds = 0.0;
};

struct GudhiPersistenceSummary {
  std::size_t finite_pairs = 0;
  std::size_t essential_intervals = 0;
};

template <class Callback>
auto time_value(Callback&& callback) {
  const auto start = Clock::now();
  auto value = callback();
  const auto stop = Clock::now();
  using Value = decltype(value);
  const std::chrono::duration<double> elapsed = stop - start;
  return TimedValue<Value>{std::move(value), elapsed.count()};
}

double milliseconds(double seconds) {
  return 1000.0 * seconds;
}

double milliseconds_from_nanoseconds(std::uint64_t nanoseconds) {
  return static_cast<double>(nanoseconds) / 1000000.0;
}

double ratio(double numerator, double denominator) {
  if (denominator == 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return numerator / denominator;
}

std::uint64_t frame_profile_nanoseconds(
    const morseframes::MorseReferenceFrameMetrics& metrics) {
  return metrics.remaining_cofaces_nanoseconds + metrics.sequence_total_nanoseconds +
         metrics.working_set_pack_nanoseconds + metrics.local_index_nanoseconds;
}

void insert_simplex(GudhiSimplexTree& simplex_tree,
                    std::vector<GudhiSimplexTree::Vertex_handle> simplex,
                    double filtration) {
  std::sort(simplex.begin(), simplex.end());
  simplex_tree.insert_simplex_and_subfaces(simplex, filtration);
}

std::size_t regular_pair_count(const morseframes::MorseSequence& sequence) {
  std::size_t count = 0;
  for (const auto& step : sequence.steps()) {
    if (step.type == morseframes::MorseStepType::RegularPair) {
      ++count;
    }
  }
  return count;
}

GudhiPersistenceSummary compute_gudhi_persistence_summary(GudhiSimplexTree& simplex_tree) {
  GudhiPersistentCohomology persistence(simplex_tree);
  persistence.init_coefficients(2);
  persistence.compute_persistent_cohomology(0.0);

  GudhiPersistenceSummary summary;
  for (const auto& interval : persistence.get_persistent_pairs()) {
    if (std::get<1>(interval) == simplex_tree.null_simplex()) {
      ++summary.essential_intervals;
    } else {
      ++summary.finite_pairs;
    }
  }
  return summary;
}

std::vector<std::tuple<std::uint16_t, double, double>> finite_barcode_key(
    const morseframes::PersistenceDiagram& diagram) {
  std::vector<std::tuple<std::uint16_t, double, double>> result;
  result.reserve(diagram.finite_pairs.size());
  for (const auto& interval : diagram.finite_pairs) {
    result.emplace_back(interval.dimension, interval.birth_value, interval.death_value);
  }
  std::sort(result.begin(), result.end());
  return result;
}

std::vector<std::pair<std::uint16_t, double>> essential_barcode_key(
    const morseframes::PersistenceDiagram& diagram) {
  std::vector<std::pair<std::uint16_t, double>> result;
  result.reserve(diagram.essential.size());
  for (const auto& interval : diagram.essential) {
    result.emplace_back(interval.dimension, interval.birth_value);
  }
  std::sort(result.begin(), result.end());
  return result;
}

template <class ComplexView>
std::vector<std::tuple<int,
                       morseframes::LevelId,
                       std::vector<morseframes::VertexId>,
                       std::vector<morseframes::VertexId>>>
sequence_signature(const ComplexView& complex, const morseframes::MorseSequence& sequence) {
  std::vector<std::tuple<int,
                         morseframes::LevelId,
                         std::vector<morseframes::VertexId>,
                         std::vector<morseframes::VertexId>>>
      result;
  result.reserve(sequence.steps().size());
  for (const auto& step : sequence.steps()) {
    result.emplace_back(step.type == morseframes::MorseStepType::Critical ? 0 : 1,
                        step.level,
                        complex.vertices(step.sigma),
                        step.type == morseframes::MorseStepType::Critical
                            ? std::vector<morseframes::VertexId>{}
                            : complex.vertices(step.tau));
  }
  return result;
}

template <class ComplexView>
morseframes::MorseReferenceReductionInput build_reduction_input(const ComplexView& complex,
                                                          const std::string& strategy,
                                                          bool collect_frame_timing = false,
                                                          morseframes::ReferenceFrameReleasePolicy
                                                              release_policy =
                                                                  morseframes::ReferenceFrameReleasePolicy::
                                                                      Eager) {
  return morseframes::build_morse_reference_reduction_input(
      complex,
      morseframes::morse_sequence_strategy_from_name(strategy),
      collect_frame_timing,
      release_policy);
}

template <class ComplexView>
morseframes::MorseSequence build_sequence(const ComplexView& complex, const std::string& strategy) {
  morseframes::FSequenceBuilder<ComplexView> builder(complex);
  switch (morseframes::morse_sequence_strategy_from_name(strategy)) {
    case morseframes::MorseSequenceStrategy::Saturated:
      return builder.build_saturated();
    case morseframes::MorseSequenceStrategy::SameLevelReduction:
      return builder.build_same_level_reduction();
    case morseframes::MorseSequenceStrategy::FMax:
      return builder.build_f_max();
    case morseframes::MorseSequenceStrategy::FMin:
      return builder.build_f_min();
    case morseframes::MorseSequenceStrategy::PlateauGreedy:
      return builder.build_plateau_greedy();
    case morseframes::MorseSequenceStrategy::FloodingMax:
      return builder.build_flooding_max();
    case morseframes::MorseSequenceStrategy::FloodingMin:
      return builder.build_flooding_min();
    case morseframes::MorseSequenceStrategy::FloodingMinMax:
      return builder.build_flooding_minmax();
    case morseframes::MorseSequenceStrategy::FloodingMaxMin:
      return builder.build_flooding_maxmin();
  }
  throw std::invalid_argument("Unknown Morse sequence strategy: " + strategy);
}

struct CaseSpec {
  std::string family;
  std::string name;
  int size = 0;
  int seed = 0;
  double radius = 0.0;
  int levels = 0;
};

struct CaseData {
  GudhiSimplexTree simplex_tree;
  std::size_t vertices = 0;
};

double grid_height(int x, int y, int levels) {
  return static_cast<double>((17 * x + 31 * y + 7 * (x / 3) + 11 * (y / 5)) % levels);
}

CaseData make_grid_case(int side, int levels) {
  GudhiSimplexTree simplex_tree;
  const auto vertex = [side](int x, int y) {
    return static_cast<GudhiSimplexTree::Vertex_handle>(y * side + x);
  };
  const auto height = [levels](int x, int y) {
    return grid_height(x, y, levels);
  };

  for (int y = 0; y < side; ++y) {
    for (int x = 0; x < side; ++x) {
      insert_simplex(simplex_tree, {vertex(x, y)}, height(x, y));
    }
  }

  for (int y = 0; y < side; ++y) {
    for (int x = 0; x < side; ++x) {
      if (x + 1 < side) {
        insert_simplex(simplex_tree,
                       {vertex(x, y), vertex(x + 1, y)},
                       std::max(height(x, y), height(x + 1, y)));
      }
      if (y + 1 < side) {
        insert_simplex(simplex_tree,
                       {vertex(x, y), vertex(x, y + 1)},
                       std::max(height(x, y), height(x, y + 1)));
      }
      if (x + 1 < side && y + 1 < side) {
        insert_simplex(simplex_tree,
                       {vertex(x + 1, y), vertex(x, y + 1)},
                       std::max(height(x + 1, y), height(x, y + 1)));
      }
    }
  }

  for (int y = 0; y + 1 < side; ++y) {
    for (int x = 0; x + 1 < side; ++x) {
      insert_simplex(simplex_tree,
                     {vertex(x, y), vertex(x + 1, y), vertex(x, y + 1)},
                     std::max({height(x, y), height(x + 1, y), height(x, y + 1)}));
      insert_simplex(simplex_tree,
                     {vertex(x + 1, y), vertex(x + 1, y + 1), vertex(x, y + 1)},
                     std::max(
                         {height(x + 1, y), height(x + 1, y + 1), height(x, y + 1)}));
    }
  }

  simplex_tree.initialize_filtration();
  return CaseData{std::move(simplex_tree), static_cast<std::size_t>(side * side)};
}

CaseData make_flag_case(int vertices, double radius, int seed, int levels) {
  struct Point {
    double x = 0.0;
    double y = 0.0;
  };

  std::mt19937 rng(static_cast<std::mt19937::result_type>(seed));
  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  std::vector<Point> points(static_cast<std::size_t>(vertices));
  std::vector<double> vertex_filtration(static_cast<std::size_t>(vertices));
  for (int vertex = 0; vertex < vertices; ++vertex) {
    points[static_cast<std::size_t>(vertex)] = Point{uniform(rng), uniform(rng)};
    vertex_filtration[static_cast<std::size_t>(vertex)] =
        static_cast<double>((vertex * 13 + seed * 5) % levels);
  }

  GudhiSimplexTree simplex_tree;
  for (int vertex = 0; vertex < vertices; ++vertex) {
    insert_simplex(simplex_tree,
                   {static_cast<GudhiSimplexTree::Vertex_handle>(vertex)},
                   vertex_filtration[static_cast<std::size_t>(vertex)]);
  }

  std::vector<std::vector<bool>> adjacent(static_cast<std::size_t>(vertices),
                                          std::vector<bool>(static_cast<std::size_t>(vertices),
                                                            false));
  std::vector<std::vector<double>> edge_filtration(
      static_cast<std::size_t>(vertices),
      std::vector<double>(static_cast<std::size_t>(vertices), 0.0));

  for (int i = 0; i < vertices; ++i) {
    for (int j = i + 1; j < vertices; ++j) {
      const double dx = points[static_cast<std::size_t>(i)].x -
                        points[static_cast<std::size_t>(j)].x;
      const double dy = points[static_cast<std::size_t>(i)].y -
                        points[static_cast<std::size_t>(j)].y;
      const double distance = std::sqrt(dx * dx + dy * dy);
      if (distance > radius) {
        continue;
      }

      const double distance_level = std::floor(4.0 * distance / radius) / 10.0;
      const double filtration =
          std::max(vertex_filtration[static_cast<std::size_t>(i)],
                   vertex_filtration[static_cast<std::size_t>(j)]) +
          distance_level;
      adjacent[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = true;
      adjacent[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = true;
      edge_filtration[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = filtration;
      edge_filtration[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = filtration;

      insert_simplex(simplex_tree,
                     {static_cast<GudhiSimplexTree::Vertex_handle>(i),
                      static_cast<GudhiSimplexTree::Vertex_handle>(j)},
                     filtration);
    }
  }

  for (int i = 0; i < vertices; ++i) {
    for (int j = i + 1; j < vertices; ++j) {
      if (!adjacent[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]) {
        continue;
      }
      for (int k = j + 1; k < vertices; ++k) {
        if (!adjacent[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] ||
            !adjacent[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)]) {
          continue;
        }

        const double filtration =
            std::max({edge_filtration[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)],
                      edge_filtration[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)],
                      edge_filtration[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)]});
        insert_simplex(simplex_tree,
                       {static_cast<GudhiSimplexTree::Vertex_handle>(i),
                        static_cast<GudhiSimplexTree::Vertex_handle>(j),
                        static_cast<GudhiSimplexTree::Vertex_handle>(k)},
                       filtration);
      }
    }
  }

  simplex_tree.initialize_filtration();
  return CaseData{std::move(simplex_tree), static_cast<std::size_t>(vertices)};
}

struct Options {
  int repeats = 3;
  bool quick = false;
  bool large = false;
  bool lean = false;
  bool canonical_order = false;
  bool deferred_release = true;
  std::vector<std::string> strategies = {
      "same-level-reduction",
      "f-max",
      "f-min",
      "plateau-greedy",
  };
};

void print_usage(const char* program) {
  std::cerr
      << "Usage: " << program
      << " [--quick|--large] [--lean] [--canonical-order]"
      << " [--deferred-release|--eager-release] [--repeats N] [--strategy NAME]\n"
      << "\n"
      << "Prints CSV rows comparing GUDHI persistent cohomology, the direct\n"
      << "Gudhi::Simplex_tree Morse view, and import into FilteredSimplicialComplex.\n"
      << "By default the direct view preserves GUDHI order inside each level/dimension;\n"
      << "--canonical-order restores lexicographic tie ordering. Reduction input uses\n"
      << "deferred annotation release by default; --eager-release restores the\n"
      << "lower-live-memory eager policy.\n";
}

Options parse_options(int argc, char** argv) {
  Options options;
  bool explicit_strategy = false;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--help" || argument == "-h") {
      print_usage(argv[0]);
      std::exit(0);
    }
    if (argument == "--quick") {
      if (options.large) {
        throw std::invalid_argument("--quick and --large are mutually exclusive.");
      }
      options.quick = true;
      options.repeats = 1;
      continue;
    }
    if (argument == "--large") {
      if (options.quick) {
        throw std::invalid_argument("--quick and --large are mutually exclusive.");
      }
      options.large = true;
      continue;
    }
    if (argument == "--lean") {
      options.lean = true;
      continue;
    }
    if (argument == "--canonical-order") {
      options.canonical_order = true;
      continue;
    }
    if (argument == "--deferred-release") {
      options.deferred_release = true;
      continue;
    }
    if (argument == "--eager-release") {
      options.deferred_release = false;
      continue;
    }
    if (argument == "--repeats") {
      if (index + 1 >= argc) {
        throw std::invalid_argument("--repeats needs an integer value.");
      }
      options.repeats = std::stoi(argv[++index]);
      if (options.repeats <= 0) {
        throw std::invalid_argument("--repeats must be positive.");
      }
      continue;
    }
    if (argument == "--strategy") {
      if (index + 1 >= argc) {
        throw std::invalid_argument("--strategy needs a strategy name.");
      }
      if (!explicit_strategy) {
        options.strategies.clear();
        explicit_strategy = true;
      }
      options.strategies.push_back(argv[++index]);
      continue;
    }
    throw std::invalid_argument("Unknown argument: " + argument);
  }
  if (options.strategies.empty()) {
    throw std::invalid_argument("At least one strategy is required.");
  }
  return options;
}

void print_header() {
  std::cout
      << "family,name,size,seed,radius,levels,vertices,simplices,filtration_levels,"
      << "strategy,repeat,view_build_ms,view_extract_ms,view_boundary_ms,"
      << "view_boundary_scan_ms,view_boundary_lookup_setup_ms,"
      << "view_boundary_register_ms,view_boundary_faces_ms,"
      << "view_coboundary_ms,view_coboundary_count_ms,"
      << "view_coboundary_reserve_ms,view_coboundary_fill_ms,"
      << "view_order_ms,view_order_size_ms,view_order_bucket_ms,"
      << "view_order_sort_emit_ms,boundary_used_dense_vertex_lookup,"
      << "boundary_used_dense_edge_lookup,boundary_vertex_count,boundary_edge_count,"
      << "boundary_dense_vertex_size,boundary_dense_edge_size,"
      << "boundary_vertex_face_lookups,boundary_edge_face_lookups,"
      << "boundary_general_face_lookups,import_ms,view_reduction_input_ms,"
      << "compact_reduction_input_ms,view_sequence_ms,compact_sequence_ms,"
      << "view_frame_extra_ms,compact_frame_extra_ms,"
      << "view_frame_remaining_cofaces_ms,view_frame_sequence_total_ms,"
      << "view_frame_sequence_core_ms,"
      << "view_sequence_init_ms,view_sequence_candidate_seed_ms,"
      << "view_sequence_candidate_loop_ms,view_sequence_emit_ms,"
      << "view_sequence_callback_ms,view_sequence_replay_ms,"
      << "view_sequence_candidate_pushes,view_sequence_candidate_pops,"
      << "view_sequence_stale_candidate_skips,"
      << "view_sequence_level_mismatch_skips,"
      << "view_sequence_regular_pairs,view_sequence_criticals,"
      << "view_frame_reference_update_ms,"
      << "view_frame_reduction_plan_ms,view_frame_release_cleanup_ms,"
      << "view_frame_working_set_pack_ms,view_frame_local_index_ms,"
      << "view_frame_peak_live_annotations,"
      << "view_frame_peak_live_total_annotation_size,"
      << "view_frame_final_live_annotations,"
      << "view_frame_final_live_total_annotation_size,"
      << "view_frame_released_total_annotation_size,"
      << "compact_frame_remaining_cofaces_ms,compact_frame_sequence_total_ms,"
      << "compact_frame_sequence_core_ms,compact_frame_reference_update_ms,"
      << "compact_frame_reduction_plan_ms,compact_frame_release_cleanup_ms,"
      << "compact_frame_working_set_pack_ms,compact_frame_local_index_ms,"
      << "view_deferred_frame_remaining_cofaces_ms,"
      << "view_deferred_frame_sequence_total_ms,"
      << "view_deferred_frame_sequence_core_ms,"
      << "view_deferred_frame_reference_update_ms,"
      << "view_deferred_frame_reduction_plan_ms,"
      << "view_deferred_frame_release_cleanup_ms,"
      << "view_deferred_frame_working_set_pack_ms,"
      << "view_deferred_frame_local_index_ms,"
      << "view_deferred_frame_peak_live_annotations,"
      << "view_deferred_frame_peak_live_total_annotation_size,"
      << "view_deferred_frame_final_live_annotations,"
      << "view_deferred_frame_final_live_total_annotation_size,"
      << "view_deferred_frame_released_total_annotation_size,"
      << "view_deferred_over_eager_frame_profile,"
      << "direct_reducer_ms,compact_reducer_ms,"
      << "direct_view_frame_total_ms,compact_import_frame_total_ms,"
      << "direct_full_ms,compact_full_ms,gudhi_persistence_ms,"
      << "gudhi_over_direct_full,gudhi_over_compact_full,view_frame_over_compact_frame,"
      << "direct_frame_total_over_import_frame_total,direct_full_over_compact_full,"
      << "view_criticals,compact_criticals,view_regular_pairs,compact_regular_pairs,"
      << "view_working_set,compact_working_set,finite_pairs,essential_intervals,"
      << "gudhi_finite_pairs,gudhi_essential_intervals,"
      << "direct_reducer_boundary_candidates,"
      << "direct_reducer_zero_boundary_criticals,"
      << "direct_reducer_boundary_xors,"
      << "direct_reducer_boundary_input_size,"
      << "direct_reducer_boundary_output_size,"
      << "direct_reducer_pivot_eliminations,"
      << "direct_inverse_initial_nonempty,"
      << "direct_inverse_initial_total_size,"
      << "direct_inverse_remove_scans,"
      << "direct_inverse_remove_applied,"
      << "direct_inverse_xor_scans,"
      << "direct_inverse_xor_applied,"
      << "direct_inverse_xor_changed_labels,"
      << "direct_inverse_xor_input_size,"
      << "direct_inverse_xor_output_size,"
      << "direct_inverse_appends\n";
}

void run_case(const CaseSpec& spec, const Options& options) {
  const auto generated = time_value([&]() {
    if (spec.family == "grid") {
      return make_grid_case(spec.size, spec.levels);
    }
    if (spec.family == "flag") {
      return make_flag_case(spec.size, spec.radius, spec.seed, spec.levels);
    }
    throw std::invalid_argument("Unknown case family: " + spec.family);
  });

  const auto& simplex_tree = generated.value.simplex_tree;
  const auto order_policy =
      options.canonical_order
          ? morseframes::SimplexTreeFiltrationOrder::CanonicalLexicographic
          : morseframes::SimplexTreeFiltrationOrder::PreserveInputWithinDimension;
  const auto release_policy =
      options.deferred_release
          ? morseframes::ReferenceFrameReleasePolicy::Deferred
          : morseframes::ReferenceFrameReleasePolicy::Eager;

  for (int repeat = 0; repeat < options.repeats; ++repeat) {
    const auto view = time_value([&]() {
      return morseframes::SimplexTreeComplexView<GudhiSimplexTree>(simplex_tree, order_policy);
    });
    const auto compact = time_value([&]() {
      return morseframes::filtered_complex_from_simplex_tree(simplex_tree);
    });
    auto gudhi_simplex_tree = simplex_tree;
    const auto gudhi_persistence = time_value([&]() {
      return compute_gudhi_persistence_summary(gudhi_simplex_tree);
    });

    if (view.value.size() != compact.value.size()) {
      throw std::logic_error("Direct view and compact import disagree on simplex count.");
    }
    if (view.value.num_levels() != compact.value.num_levels()) {
      throw std::logic_error("Direct view and compact import disagree on filtration levels.");
    }

    for (const std::string& strategy : options.strategies) {
      const auto& view_build_metrics = view.value.build_metrics();
      auto view_input = time_value([&]() {
        return build_reduction_input(view.value, strategy, false, release_policy);
      });
      auto compact_input = time_value([&]() {
        return build_reduction_input(compact.value, strategy, false, release_policy);
      });
      const std::size_t view_criticals =
          view_input.value.sequence.critical_simplices().size();
      const std::size_t compact_criticals =
          compact_input.value.sequence.critical_simplices().size();
      const std::size_t view_regular_pairs =
          regular_pair_count(view_input.value.sequence);
      const std::size_t compact_regular_pairs =
          regular_pair_count(compact_input.value.sequence);
      const std::size_t view_working_set =
          view_input.value.reduction_plan.working_set.size();
      const std::size_t compact_working_set =
          compact_input.value.reduction_plan.working_set.size();

      double view_sequence_seconds = 0.0;
      double compact_sequence_seconds = 0.0;
      morseframes::MorseReferenceFrameMetrics view_frame_metrics =
          view_input.value.frame_metrics;
      morseframes::MorseReferenceFrameMetrics compact_frame_metrics =
          compact_input.value.frame_metrics;
      morseframes::MorseReferenceFrameMetrics view_deferred_frame_metrics;
      double view_frame_profile_seconds = 0.0;
      double view_deferred_frame_profile_seconds = 0.0;

      if (!options.lean) {
        auto view_sequence = time_value([&]() {
          return build_sequence(view.value, strategy);
        });
        auto compact_sequence = time_value([&]() {
          return build_sequence(compact.value, strategy);
        });
        auto view_frame_profile =
            build_reduction_input(view.value, strategy, true, release_policy);
        auto compact_frame_profile =
            build_reduction_input(compact.value, strategy, true, release_policy);
        auto view_deferred_frame_profile =
            build_reduction_input(view.value,
                                  strategy,
                                  true,
                                  morseframes::ReferenceFrameReleasePolicy::Deferred);

        if (sequence_signature(view.value, view_sequence.value) !=
            sequence_signature(view.value, view_input.value.sequence)) {
          throw std::logic_error(
              "Direct view sequence-only and fused reduction-input sequences disagree.");
        }
        if (sequence_signature(compact.value, compact_sequence.value) !=
            sequence_signature(compact.value, compact_input.value.sequence)) {
          throw std::logic_error(
              "Compact import sequence-only and fused reduction-input sequences disagree.");
        }
        if (sequence_signature(view.value, view_input.value.sequence) !=
            sequence_signature(view.value, view_frame_profile.sequence)) {
          throw std::logic_error(
              "Direct view timed and untimed fused reduction-input sequences disagree.");
        }
        if (sequence_signature(compact.value, compact_input.value.sequence) !=
            sequence_signature(compact.value, compact_frame_profile.sequence)) {
          throw std::logic_error(
              "Compact import timed and untimed fused reduction-input sequences disagree.");
        }
        if (sequence_signature(view.value, view_input.value.sequence) !=
            sequence_signature(view.value, view_deferred_frame_profile.sequence)) {
          throw std::logic_error(
              "Direct view eager and deferred fused reduction-input sequences disagree.");
        }

        view_sequence_seconds = view_sequence.seconds;
        compact_sequence_seconds = compact_sequence.seconds;
        view_frame_metrics = view_frame_profile.frame_metrics;
        compact_frame_metrics = compact_frame_profile.frame_metrics;
        view_deferred_frame_metrics = view_deferred_frame_profile.frame_metrics;
        view_frame_profile_seconds =
            static_cast<double>(frame_profile_nanoseconds(view_frame_metrics)) /
            1000000000.0;
        view_deferred_frame_profile_seconds =
            static_cast<double>(frame_profile_nanoseconds(view_deferred_frame_metrics)) /
            1000000000.0;
      }

      TimedValue<morseframes::MorseReferenceReductionResult> direct_reducer;
      TimedValue<morseframes::MorseReferenceReductionResult> compact_reducer;
      if (options.lean) {
        direct_reducer = time_value([&]() {
          morseframes::MorseReferencePersistenceReducer persistence_reducer(
              view.value,
              view_input.value.sequence,
              std::move(view_input.value.reduction_plan),
              std::move(view_input.value.annotations));
          return persistence_reducer.compute_with_metrics();
        });
        compact_reducer = time_value([&]() {
          morseframes::MorseReferencePersistenceReducer persistence_reducer(
              compact.value,
              compact_input.value.sequence,
              std::move(compact_input.value.reduction_plan),
              std::move(compact_input.value.annotations));
          return persistence_reducer.compute_with_metrics();
        });
      } else {
        auto direct_reducer_input =
            build_reduction_input(view.value, strategy, false, release_policy);
        direct_reducer = time_value([&]() {
          morseframes::MorseReferencePersistenceReducer persistence_reducer(
              view.value,
              direct_reducer_input.sequence,
              std::move(direct_reducer_input.reduction_plan),
              std::move(direct_reducer_input.annotations));
          return persistence_reducer.compute_with_metrics();
        });

        auto compact_reducer_input =
            build_reduction_input(compact.value, strategy, false, release_policy);
        compact_reducer = time_value([&]() {
          morseframes::MorseReferencePersistenceReducer persistence_reducer(
              compact.value,
              compact_reducer_input.sequence,
              std::move(compact_reducer_input.reduction_plan),
              std::move(compact_reducer_input.annotations));
          return persistence_reducer.compute_with_metrics();
        });
      }

      if (finite_barcode_key(direct_reducer.value.diagram) !=
              finite_barcode_key(compact_reducer.value.diagram) ||
          essential_barcode_key(direct_reducer.value.diagram) !=
              essential_barcode_key(compact_reducer.value.diagram)) {
        throw std::logic_error("Direct view and compact import persistence disagree.");
      }

      const double direct_frame_total = view.seconds + view_input.seconds;
      const double compact_frame_total = compact.seconds + compact_input.seconds;
      const double view_frame_extra = options.lean
                                          ? 0.0
                                          : std::max(0.0,
                                                     view_input.seconds -
                                                         view_sequence_seconds);
      const double compact_frame_extra = options.lean
                                             ? 0.0
                                             : std::max(0.0,
                                                        compact_input.seconds -
                                                            compact_sequence_seconds);
      const double direct_full =
          view.seconds + view_input.seconds + direct_reducer.seconds;
      const double compact_full =
          compact.seconds + compact_input.seconds + compact_reducer.seconds;

      std::cout << spec.family << ',' << spec.name << ',' << spec.size << ',' << spec.seed
                << ',' << spec.radius << ',' << spec.levels << ','
                << generated.value.vertices << ',' << view.value.size() << ','
                << view.value.num_levels() << ',' << strategy << ',' << repeat << ','
                << milliseconds(view.seconds) << ','
                << milliseconds_from_nanoseconds(view_build_metrics.extract_nanoseconds) << ','
                << milliseconds_from_nanoseconds(view_build_metrics.boundary_nanoseconds) << ','
                << milliseconds_from_nanoseconds(view_build_metrics.boundary_scan_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_build_metrics.boundary_lookup_setup_nanoseconds)
                << ','
                << milliseconds_from_nanoseconds(
                       view_build_metrics.boundary_register_nanoseconds)
                << ','
                << milliseconds_from_nanoseconds(
                       view_build_metrics.boundary_face_lookup_nanoseconds)
                << ','
                << milliseconds_from_nanoseconds(view_build_metrics.coboundary_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_build_metrics.coboundary_count_nanoseconds)
                << ','
                << milliseconds_from_nanoseconds(
                       view_build_metrics.coboundary_reserve_nanoseconds)
                << ','
                << milliseconds_from_nanoseconds(
                       view_build_metrics.coboundary_fill_nanoseconds)
                << ','
                << milliseconds_from_nanoseconds(view_build_metrics.order_nanoseconds) << ','
                << milliseconds_from_nanoseconds(view_build_metrics.order_size_nanoseconds) << ','
                << milliseconds_from_nanoseconds(view_build_metrics.order_bucket_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_build_metrics.order_sort_emit_nanoseconds)
                << ','
                << (view_build_metrics.boundary_used_dense_vertex_lookup ? 1 : 0) << ','
                << (view_build_metrics.boundary_used_dense_edge_lookup ? 1 : 0) << ','
                << view_build_metrics.boundary_vertex_count << ','
                << view_build_metrics.boundary_edge_count << ','
                << view_build_metrics.boundary_dense_vertex_size << ','
                << view_build_metrics.boundary_dense_edge_size << ','
                << view_build_metrics.boundary_vertex_face_lookups << ','
                << view_build_metrics.boundary_edge_face_lookups << ','
                << view_build_metrics.boundary_general_face_lookups << ','
                << milliseconds(compact.seconds) << ','
                << milliseconds(view_input.seconds) << ','
                << milliseconds(compact_input.seconds) << ','
                << milliseconds(view_sequence_seconds) << ','
                << milliseconds(compact_sequence_seconds) << ','
                << milliseconds(view_frame_extra) << ','
                << milliseconds(compact_frame_extra) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.remaining_cofaces_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.sequence_total_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.sequence_core_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.sequence_init_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.sequence_candidate_seed_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.sequence_candidate_loop_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.sequence_emit_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.sequence_callback_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.sequence_replay_nanoseconds) << ','
                << view_frame_metrics.sequence_candidate_pushes << ','
                << view_frame_metrics.sequence_candidate_pops << ','
                << view_frame_metrics.sequence_stale_candidate_skips << ','
                << view_frame_metrics.sequence_level_mismatch_skips << ','
                << view_frame_metrics.sequence_regular_pairs << ','
                << view_frame_metrics.sequence_criticals << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.reference_update_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.reduction_plan_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.release_cleanup_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.working_set_pack_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_frame_metrics.local_index_nanoseconds) << ','
                << view_frame_metrics.peak_live_nonempty_annotations << ','
                << view_frame_metrics.peak_live_total_annotation_size << ','
                << view_frame_metrics.final_live_nonempty_annotations << ','
                << view_frame_metrics.final_live_total_annotation_size << ','
                << view_frame_metrics.released_total_annotation_size << ','
                << milliseconds_from_nanoseconds(
                       compact_frame_metrics.remaining_cofaces_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       compact_frame_metrics.sequence_total_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       compact_frame_metrics.sequence_core_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       compact_frame_metrics.reference_update_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       compact_frame_metrics.reduction_plan_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       compact_frame_metrics.release_cleanup_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       compact_frame_metrics.working_set_pack_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       compact_frame_metrics.local_index_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_deferred_frame_metrics.remaining_cofaces_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_deferred_frame_metrics.sequence_total_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_deferred_frame_metrics.sequence_core_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_deferred_frame_metrics.reference_update_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_deferred_frame_metrics.reduction_plan_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_deferred_frame_metrics.release_cleanup_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_deferred_frame_metrics.working_set_pack_nanoseconds) << ','
                << milliseconds_from_nanoseconds(
                       view_deferred_frame_metrics.local_index_nanoseconds) << ','
                << view_deferred_frame_metrics.peak_live_nonempty_annotations << ','
                << view_deferred_frame_metrics.peak_live_total_annotation_size << ','
                << view_deferred_frame_metrics.final_live_nonempty_annotations << ','
                << view_deferred_frame_metrics.final_live_total_annotation_size << ','
                << view_deferred_frame_metrics.released_total_annotation_size << ','
                << ratio(view_deferred_frame_profile_seconds,
                         view_frame_profile_seconds) << ','
                << milliseconds(direct_reducer.seconds) << ','
                << milliseconds(compact_reducer.seconds) << ','
                << milliseconds(direct_frame_total) << ','
                << milliseconds(compact_frame_total) << ','
                << milliseconds(direct_full) << ',' << milliseconds(compact_full) << ','
                << milliseconds(gudhi_persistence.seconds) << ','
                << ratio(gudhi_persistence.seconds, direct_full) << ','
                << ratio(gudhi_persistence.seconds, compact_full) << ','
                << ratio(view_input.seconds, compact_input.seconds) << ','
                << ratio(direct_frame_total, compact_frame_total) << ','
                << ratio(direct_full, compact_full) << ','
                << view_criticals << ','
                << compact_criticals << ','
                << view_regular_pairs << ','
                << compact_regular_pairs << ','
                << view_working_set << ','
                << compact_working_set << ','
                << compact_reducer.value.metrics.finite_pairs << ','
                << compact_reducer.value.metrics.essential_intervals << ','
                << gudhi_persistence.value.finite_pairs << ','
                << gudhi_persistence.value.essential_intervals << ','
                << direct_reducer.value.metrics.boundary_annotation_candidate_criticals << ','
                << direct_reducer.value.metrics.boundary_annotation_zero_skipped_criticals << ','
                << direct_reducer.value.metrics.boundary_annotation_xors << ','
                << direct_reducer.value.metrics.boundary_annotation_total_input_size << ','
                << direct_reducer.value.metrics.boundary_annotation_total_output_size << ','
                << direct_reducer.value.metrics.pivot_eliminations << ','
                << direct_reducer.value.metrics.inverse_store.initial_nonempty_annotations << ','
                << direct_reducer.value.metrics.inverse_store.initial_total_annotation_size << ','
                << direct_reducer.value.metrics.inverse_store.remove_candidate_scans << ','
                << direct_reducer.value.metrics.inverse_store.remove_applied << ','
                << direct_reducer.value.metrics.inverse_store.xor_candidate_scans << ','
                << direct_reducer.value.metrics.inverse_store.xor_applied << ','
                << direct_reducer.value.metrics.inverse_store.xor_changed_labels << ','
                << direct_reducer.value.metrics.inverse_store.xor_total_input_size << ','
                << direct_reducer.value.metrics.inverse_store.xor_total_output_size << ','
                << direct_reducer.value.metrics.inverse_store.inverse_list_appends << '\n';
    }
  }
}

std::vector<CaseSpec> case_specs(const Options& options) {
  if (options.large) {
    return {
        CaseSpec{"grid", "grid-64x64-plateau", 64, 0, 0.0, 7},
        CaseSpec{"grid", "grid-96x96-plateau", 96, 0, 0.0, 9},
        CaseSpec{"grid", "grid-128x128-plateau", 128, 0, 0.0, 11},
        CaseSpec{"flag", "flag-240-r014", 240, 271828, 0.14, 7},
    };
  }

  if (options.quick) {
    return {
        CaseSpec{"grid", "grid-24x24-plateau", 24, 0, 0.0, 5},
        CaseSpec{"flag", "flag-90-r018", 90, 1729, 0.18, 5},
    };
  }

  return {
      CaseSpec{"grid", "grid-32x32-plateau", 32, 0, 0.0, 5},
      CaseSpec{"grid", "grid-48x48-plateau", 48, 0, 0.0, 7},
      CaseSpec{"flag", "flag-120-r018", 120, 1729, 0.18, 5},
      CaseSpec{"flag", "flag-160-r016", 160, 314159, 0.16, 7},
  };
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_options(argc, argv);

    std::cout << std::fixed << std::setprecision(6);
    print_header();
    for (const auto& spec : case_specs(options)) {
      run_case(spec, options);
    }
  } catch (const std::exception& error) {
    std::cerr << "benchmark_gudhi_view: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
