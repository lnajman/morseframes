#pragma once

#include "morseframes/filtered_complex.hpp"
#include <array>
#include <iterator>
#include <numeric>

namespace morseframes {

struct LowerStarConstructionMetrics {
  double validation_seconds = 0;
  double enumeration_seconds = 0;
  double sort_and_dedup_seconds = 0;
  double insertion_seconds = 0;
  std::size_t generated_faces = 0;
  std::size_t unique_faces_submitted = 0;
};

namespace detail {
struct LowerStarComplexBuilder {
  using Cells = std::vector<std::vector<VertexId>>;
  using Clock = std::chrono::steady_clock;

  template <class Vertices> struct Record {
    Vertices vertices;
    std::size_t first_cell = 0;
  };

  static std::size_t combinations(std::size_t n, std::size_t k) {
    if (k > n) return 0;
    k = std::min(k, n - k);
    std::size_t count = 1;
    for (std::size_t i = 1; i <= k; ++i) {
      auto numerator = n - k + i, denominator = i;
      const auto divisor = std::gcd(numerator, denominator);
      numerator /= divisor; denominator /= divisor;
      count /= denominator; // Remaining denominator divides the prior binomial.
      if (count > std::numeric_limits<std::size_t>::max() / numerator)
        throw std::length_error("Too many faces in lower-star construction.");
      count *= numerator;
    }
    return count;
  }

  template <bool Diagnostic>
  static void record(Clock::time_point& last, double* destination) {
    if constexpr (Diagnostic) {
      const auto now = Clock::now();
      *destination += std::chrono::duration<double>(now - last).count();
      last = now;
    }
  }

  template <class Vertices>
  static void enumerate(const std::vector<VertexId>& cell, std::size_t cell_index,
                        Vertices& face, std::size_t first, std::size_t depth,
                        std::vector<Record<Vertices>>& faces) {
    if (depth == face.size()) {
      faces.push_back({face, cell_index});
      return;
    }
    const auto needed = face.size() - depth;
    for (std::size_t i = first; i + needed <= cell.size(); ++i) {
      face[depth] = cell[i];
      enumerate(cell, cell_index, face, i + 1, depth + 1, faces);
    }
  }

  template <bool Diagnostic, class Vertices>
  static void dimension(FilteredSimplicialComplex& complex,
                        const std::vector<double>& values, const Cells& cells,
                        const Cells& canonical, Vertices face,
                        LowerStarConstructionMetrics* metrics) {
    Clock::time_point last;
    if constexpr (Diagnostic) last = Clock::now();
    std::vector<Record<Vertices>> faces;
    std::size_t count = 0;
    for (const auto& cell : canonical) {
      const auto add = combinations(cell.size(), face.size());
      if (add > faces.max_size() - count)
        throw std::length_error("Lower-star face buffer exceeds its maximum size.");
      count += add;
    }
    faces.reserve(count);
    for (std::size_t i = 0; i < canonical.size(); ++i)
      if (canonical[i].size() >= face.size())
        enumerate(canonical[i], i, face, 0, 0, faces);
    if constexpr (Diagnostic) {
      if (count > std::numeric_limits<std::size_t>::max() - metrics->generated_faces)
        throw std::length_error("Lower-star diagnostic face count overflow.");
      metrics->generated_faces += count;
      record<true>(last, &metrics->enumeration_seconds);
    }
    std::sort(faces.begin(), faces.end(), [](const auto& a, const auto& b) {
      return a.vertices < b.vertices;
    });
    // Equal keys retain their earliest input cell, including signed-zero ties
    // in max-vertex extension. No reliance on an unstable sort's tie order.
    std::size_t unique = 0;
    for (std::size_t i = 0; i < faces.size(); ++i) {
      if (unique && faces[unique - 1].vertices == faces[i].vertices) {
        faces[unique - 1].first_cell = std::min(faces[unique - 1].first_cell, faces[i].first_cell);
      } else {
        if (unique != i) faces[unique] = std::move(faces[i]);
        ++unique;
      }
    }
    faces.resize(unique);
    if constexpr (Diagnostic) {
      metrics->unique_faces_submitted += unique; // Bounded by generated_faces.
      record<true>(last, &metrics->sort_and_dedup_seconds);
    }
    auto hint = complex.pending_.begin();
    for (const auto& item : faces) {
      std::vector<VertexId> vertices(item.vertices.begin(), item.vertices.end());
      double value = -std::numeric_limits<double>::infinity();
      for (VertexId vertex : vertices) value = std::max(value, values[vertex]);
      if (value == 0.0) {
        // Preserve the legacy extension's first encountered zero, not just its
        // numerical value, even when input cell vertices are not sorted.
        value = -std::numeric_limits<double>::infinity();
        for (VertexId vertex : cells[item.first_cell])
          if (std::binary_search(vertices.begin(), vertices.end(), vertex))
            value = std::max(value, values[vertex]);
      }
      // Monotone merge, not a new tree search for every unique face. Existing
      // pending simplices may have other dimensions and are visited once here.
      while (hint != complex.pending_.end() && hint->first < vertices) ++hint;
      if (hint != complex.pending_.end() && hint->first == vertices) {
        if (std::fabs(hint->second - value) > 1e-12)
          throw std::invalid_argument("Duplicate simplex inserted with a different filtration value.");
      } else {
        hint = complex.pending_.emplace_hint(hint, std::move(vertices), value);
      }
      ++hint;
    }
    // Include destruction of temporary records in the insertion/cleanup phase.
    std::vector<Record<Vertices>>().swap(faces);
    if constexpr (Diagnostic) record<true>(last, &metrics->insertion_seconds);
  }

  template <bool Diagnostic>
  static void build(FilteredSimplicialComplex& complex,
                    const std::vector<double>& values, const Cells& cells,
                    LowerStarConstructionMetrics* metrics) {
    Clock::time_point last;
    if constexpr (Diagnostic) { *metrics = {}; last = Clock::now(); }
    for (double value : values)
      if (std::isnan(value)) throw std::invalid_argument("Lower-star vertex values cannot be NaN.");
    Cells canonical = cells;
    std::size_t maximum = 0;
    for (auto& cell : canonical) {
      FilteredSimplicialComplex::canonicalize(cell);
      if (cell.empty()) throw std::invalid_argument("A simplex must contain at least one vertex.");
      FilteredSimplicialComplex::checked_dimension(cell.size() - 1);
      for (VertexId vertex : cell)
        if (vertex >= values.size()) throw std::out_of_range("Missing lower-star vertex value.");
      maximum = std::max(maximum, cell.size());
    }
    if constexpr (Diagnostic) record<true>(last, &metrics->validation_seconds);
    if (maximum) complex.clear_same_level_closure_cache();
    for (std::size_t cardinality = 1; cardinality <= maximum; ++cardinality) {
      switch (cardinality) {
        case 1: dimension<Diagnostic>(complex, values, cells, canonical, std::array<VertexId, 1>{}, metrics); break;
        case 2: dimension<Diagnostic>(complex, values, cells, canonical, std::array<VertexId, 2>{}, metrics); break;
        case 3: dimension<Diagnostic>(complex, values, cells, canonical, std::array<VertexId, 3>{}, metrics); break;
        case 4: dimension<Diagnostic>(complex, values, cells, canonical, std::array<VertexId, 4>{}, metrics); break;
        default: dimension<Diagnostic>(complex, values, cells, canonical, std::vector<VertexId>(cardinality), metrics);
      }
    }
    if constexpr (Diagnostic) last = Clock::now();
    Cells().swap(canonical);
    if constexpr (Diagnostic) record<true>(last, &metrics->insertion_seconds);
  }
};
} // namespace detail

// Insert the closure of cells with max-vertex filtration. Vertex IDs index
// values; use singleton cells to include isolated vertices. Cells may be mixed
// dimensional, duplicated or unordered. Existing equal-weight entries retain
// their values. Call finalize() afterward. Errors during insertion can leave
// partial additions, as with a sequence of add_simplex() calls.
inline void add_lower_star_cells(FilteredSimplicialComplex& complex,
                                 const std::vector<double>& values,
                                 const std::vector<std::vector<VertexId>>& cells) {
  detail::LowerStarComplexBuilder::build<false>(complex, values, cells, nullptr);
}

inline void add_lower_star_cells_with_metrics(
    FilteredSimplicialComplex& complex, const std::vector<double>& values,
    const std::vector<std::vector<VertexId>>& cells, LowerStarConstructionMetrics& metrics) {
  detail::LowerStarComplexBuilder::build<true>(complex, values, cells, &metrics);
}

} // namespace morseframes
