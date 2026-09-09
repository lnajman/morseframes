#include <DiscreteGradient.h>
#include <Triangulation.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct InputComplex {
  int dimension{};
  ttk::SimplexId vertex_count{};
  ttk::SimplexId cell_count{};
  std::vector<double> scalars;
  std::vector<ttk::SimplexId> offsets;
  std::vector<float> points;
  std::vector<ttk::LongSimplexId> cells;
};

double seconds(Clock::time_point first, Clock::time_point last) {
  return std::chrono::duration<double>(last - first).count();
}

InputComplex read_complex(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("Cannot open input complex: " + path);
  }

  std::string magic;
  InputComplex result;
  input >> magic >> result.dimension >> result.vertex_count >> result.cell_count;
  if (!input || magic != "morseframes-ttk-v1") {
    throw std::runtime_error("Invalid TTK benchmark input header.");
  }
  if (result.dimension < 1 || result.dimension > 3 || result.vertex_count < 1 ||
      result.cell_count < 1) {
    throw std::runtime_error("TTK benchmark inputs must be nonempty and 1D, 2D, or 3D.");
  }

  result.scalars.resize(static_cast<std::size_t>(result.vertex_count));
  result.offsets.resize(static_cast<std::size_t>(result.vertex_count));
  result.points.assign(static_cast<std::size_t>(3 * result.vertex_count), 0.0F);
  std::vector<std::pair<double, ttk::SimplexId>> scalar_order;
  scalar_order.reserve(static_cast<std::size_t>(result.vertex_count));
  for (ttk::SimplexId vertex = 0; vertex < result.vertex_count; ++vertex) {
    input >> result.scalars[static_cast<std::size_t>(vertex)];
    result.points[static_cast<std::size_t>(3 * vertex)] = static_cast<float>(vertex);
    scalar_order.emplace_back(result.scalars[static_cast<std::size_t>(vertex)], vertex);
  }
  std::sort(scalar_order.begin(), scalar_order.end());
  for (std::size_t rank = 0; rank < scalar_order.size(); ++rank) {
    if (rank > 0 && scalar_order[rank - 1].first == scalar_order[rank].first) {
      throw std::runtime_error("TTK ProcessLowerStars benchmark requires injective values.");
    }
    result.offsets[static_cast<std::size_t>(scalar_order[rank].second)] =
        static_cast<ttk::SimplexId>(rank);
  }

  const auto vertices_per_cell = static_cast<ttk::LongSimplexId>(result.dimension + 1);
  result.cells.reserve(
      static_cast<std::size_t>(result.cell_count) *
      static_cast<std::size_t>(vertices_per_cell + 1));
  for (ttk::SimplexId cell = 0; cell < result.cell_count; ++cell) {
    result.cells.push_back(vertices_per_cell);
    for (ttk::LongSimplexId local = 0; local < vertices_per_cell; ++local) {
      ttk::LongSimplexId vertex{};
      input >> vertex;
      if (vertex < 0 || vertex >= result.vertex_count) {
        throw std::runtime_error("Cell contains an invalid vertex identifier.");
      }
      result.cells.push_back(vertex);
    }
  }
  if (!input) {
    throw std::runtime_error("TTK benchmark input ended before all cells were read.");
  }
  return result;
}

struct Options {
  std::string input_path;
  int workers{1};
  int repeats{5};
  int warmups{1};
};

Options parse_options(int argc, char** argv) {
  Options result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (index + 1 >= argc) {
      throw std::invalid_argument("Missing value after " + option);
    }
    const std::string value = argv[++index];
    if (option == "--input") {
      result.input_path = value;
    } else if (option == "--workers") {
      result.workers = std::stoi(value);
    } else if (option == "--repeats") {
      result.repeats = std::stoi(value);
    } else if (option == "--warmups") {
      result.warmups = std::stoi(value);
    } else {
      throw std::invalid_argument("Unknown option: " + option);
    }
  }
  if (result.input_path.empty() || result.workers < 1 || result.repeats < 1 ||
      result.warmups < 0) {
    throw std::invalid_argument(
        "Usage: benchmark --input FILE --workers N --repeats N --warmups N");
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_options(argc, argv);
    InputComplex input = read_complex(options.input_path);

    const auto setup_start = Clock::now();
    ttk::Triangulation triangulation;
    if (triangulation.setInputPoints(input.vertex_count, input.points.data()) != 0 ||
        triangulation.setInputCells(input.cell_count, input.cells.data()) != 0) {
      throw std::runtime_error("TTK rejected the explicit triangulation.");
    }
    const auto setup_stop = Clock::now();

    ttk::dcg::DiscreteGradient gradient;
    gradient.setDebugLevel(0);
    gradient.setThreadNumber(options.workers);
    gradient.setBackend(ttk::dcg::DiscreteGradient::BACKEND::CLASSIC_BACKEND);
    gradient.setInputScalarField(input.scalars.data(), 1);
    gradient.setInputOffsets(input.offsets.data());

    const auto precondition_start = Clock::now();
    gradient.preconditionTriangulation(&triangulation);
    const auto precondition_stop = Clock::now();

    auto build_gradient = [&]() {
      const int status = gradient.buildGradient<double>(triangulation, true);
      if (status != 0) {
        throw std::runtime_error("TTK ProcessLowerStars gradient construction failed.");
      }
    };
    for (int repeat = 0; repeat < options.warmups; ++repeat) {
      build_gradient();
    }

    double best_gradient_seconds = std::numeric_limits<double>::infinity();
    for (int repeat = 0; repeat < options.repeats; ++repeat) {
      const auto start = Clock::now();
      build_gradient();
      const auto stop = Clock::now();
      best_gradient_seconds = std::min(best_gradient_seconds, seconds(start, stop));
    }

    std::vector<ttk::SimplexId> critical_by_dimension(
        static_cast<std::size_t>(input.dimension + 1), 0);
    ttk::SimplexId critical_total = 0;
    ttk::SimplexId simplex_total = 0;
    for (int dimension = 0; dimension <= input.dimension; ++dimension) {
      const ttk::SimplexId count = gradient.getNumberOfCells(dimension, triangulation);
      simplex_total += count;
      for (ttk::SimplexId cell = 0; cell < count; ++cell) {
        if (gradient.isCellCritical(dimension, cell)) {
          ++critical_by_dimension[static_cast<std::size_t>(dimension)];
          ++critical_total;
        }
      }
    }

    std::cout << std::setprecision(17)
              << "{\"dimension\":" << input.dimension
              << ",\"num_vertices\":" << input.vertex_count
              << ",\"num_simplices\":" << simplex_total
              << ",\"num_critical_simplices\":" << critical_total
              << ",\"critical_simplices_by_dimension\":[";
    for (std::size_t dimension = 0; dimension < critical_by_dimension.size();
         ++dimension) {
      if (dimension != 0) {
        std::cout << ',';
      }
      std::cout << critical_by_dimension[dimension];
    }
    std::cout << "]"
              << ",\"workers\":" << options.workers
              << ",\"repeats\":" << options.repeats
              << ",\"warmups\":" << options.warmups
              << ",\"setup_seconds\":" << seconds(setup_start, setup_stop)
              << ",\"precondition_seconds\":"
              << seconds(precondition_start, precondition_stop)
              << ",\"gradient_seconds\":" << best_gradient_seconds << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
