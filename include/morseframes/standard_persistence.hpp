#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "morseframes/field_arithmetic.hpp"
#include "morseframes/filtered_complex.hpp"
#include "morseframes/incidence.hpp"
#include "morseframes/reference_persistence.hpp"

namespace morseframes {

using MatrixColumn = std::vector<std::uint32_t>;

struct ModPEntry {
  std::uint32_t row = 0;
  std::uint32_t coefficient = 0;
};

using ModPColumn = std::vector<ModPEntry>;

inline void xor_columns_in_place(MatrixColumn& lhs, const MatrixColumn& rhs) {
  MatrixColumn result;
  result.reserve(lhs.size() + rhs.size());

  auto left = lhs.begin();
  auto right = rhs.begin();
  while (left != lhs.end() || right != rhs.end()) {
    if (right == rhs.end() || (left != lhs.end() && *left < *right)) {
      result.push_back(*left);
      ++left;
    } else if (left == lhs.end() || *right < *left) {
      result.push_back(*right);
      ++right;
    } else {
      ++left;
      ++right;
    }
  }

  lhs = std::move(result);
}

inline void modp_add_scaled_in_place(ModPColumn& target,
                                     const ModPColumn& source,
                                     std::uint32_t scale,
                                     std::uint32_t modulus) {
  scale %= modulus;
  if (scale == 0 || source.empty()) {
    return;
  }

  ModPColumn result;
  result.reserve(target.size() + source.size());
  std::size_t left = 0;
  std::size_t right = 0;

  while (left < target.size() || right < source.size()) {
    if (right == source.size() ||
        (left < target.size() && target[left].row < source[right].row)) {
      result.push_back(target[left]);
      ++left;
      continue;
    }

    if (left == target.size() || source[right].row < target[left].row) {
      const auto coefficient = modp_multiply(scale, source[right].coefficient, modulus);
      if (coefficient != 0) {
        result.push_back(ModPEntry{source[right].row, coefficient});
      }
      ++right;
      continue;
    }

    const auto scaled = modp_multiply(scale, source[right].coefficient, modulus);
    const auto coefficient = static_cast<std::uint32_t>(
        (static_cast<std::uint64_t>(target[left].coefficient) + scaled) % modulus);
    if (coefficient != 0) {
      result.push_back(ModPEntry{target[left].row, coefficient});
    }
    ++left;
    ++right;
  }

  target = std::move(result);
}

template <class ComplexView>
inline ModPColumn oriented_boundary_column_modp(
    const ComplexView& complex,
    SimplexId simplex,
    const std::vector<std::uint32_t>& order_index,
    std::uint32_t column_index,
    std::uint32_t modulus) {
  static_assert(is_complex_view_v<ComplexView>,
                "oriented_boundary_column_modp requires a Morse complex-view type.");
  const auto& boundary = complex.boundary(simplex);
  ModPColumn column;
  column.reserve(boundary.size());
  for (std::size_t removed_index = 0; removed_index < boundary.size(); ++removed_index) {
    const std::uint32_t row = order_index.at(boundary[removed_index]);
    if (row >= column_index) {
      throw std::logic_error("Filtration order does not place faces before cofaces.");
    }
    column.push_back(ModPEntry{
        row,
        boundary_incidence_coefficient(complex, simplex, removed_index, modulus),
    });
  }
  std::sort(column.begin(), column.end(), [](const ModPEntry& lhs, const ModPEntry& rhs) {
    return lhs.row < rhs.row;
  });
  return column;
}

template <class ComplexView>
inline PersistenceDiagram compute_standard_z2_persistence(
    const ComplexView& complex) {
  static_assert(is_complex_view_v<ComplexView>,
                "compute_standard_z2_persistence requires a Morse complex-view type.");
  const auto& order = complex.filtration_order();
  const std::size_t n = order.size();
  const std::uint32_t invalid = std::numeric_limits<std::uint32_t>::max();

  std::vector<std::uint32_t> order_index(complex.size(), invalid);
  for (std::uint32_t i = 0; i < n; ++i) {
    order_index.at(order[i]) = i;
  }

  std::vector<MatrixColumn> reduced_columns(n);
  std::vector<std::uint32_t> low_to_column(n, invalid);
  std::vector<bool> is_birth(n, false);
  std::vector<bool> is_killed(n, false);

  PersistenceDiagram diagram;

  for (std::uint32_t column_index = 0; column_index < n; ++column_index) {
    const SimplexId sigma = order[column_index];
    MatrixColumn column;
    column.reserve(complex.boundary(sigma).size());
    for (SimplexId face : complex.boundary(sigma)) {
      const std::uint32_t row = order_index.at(face);
      if (row >= column_index) {
        throw std::logic_error("Filtration order does not place faces before cofaces.");
      }
      column.push_back(row);
    }
    std::sort(column.begin(), column.end());

    while (!column.empty()) {
      const std::uint32_t low = column.back();
      if (low_to_column[low] == invalid) {
        break;
      }
      xor_columns_in_place(column, reduced_columns[low_to_column[low]]);
    }

    if (column.empty()) {
      is_birth[column_index] = true;
      continue;
    }

    const std::uint32_t low = column.back();
    low_to_column[low] = column_index;
    reduced_columns[column_index] = std::move(column);
    is_killed[low] = true;

    const SimplexId birth = order[low];
    diagram.finite_pairs.push_back(PersistencePair{
        birth,
        sigma,
        complex.dimension(birth),
        complex.filtration(birth),
        complex.filtration(sigma),
    });
  }

  for (std::uint32_t column_index = 0; column_index < n; ++column_index) {
    if (!is_birth[column_index] || is_killed[column_index]) {
      continue;
    }

    const SimplexId birth = order[column_index];
    diagram.essential.push_back(EssentialInterval{
        birth,
        complex.dimension(birth),
        complex.filtration(birth),
    });
  }

  return diagram;
}

template <class ComplexView>
inline PersistenceDiagram compute_standard_prime_field_persistence(
    const ComplexView& complex,
    std::uint32_t modulus) {
  static_assert(is_complex_view_v<ComplexView>,
                "compute_standard_prime_field_persistence requires a Morse complex-view type.");
  validate_prime_field_characteristic(modulus);
  if (modulus == 2) {
    return compute_standard_z2_persistence(complex);
  }

  const auto& order = complex.filtration_order();
  const std::size_t n = order.size();
  const std::uint32_t invalid = std::numeric_limits<std::uint32_t>::max();

  std::vector<std::uint32_t> order_index(complex.size(), invalid);
  for (std::uint32_t i = 0; i < n; ++i) {
    order_index.at(order[i]) = i;
  }

  std::vector<ModPColumn> reduced_columns(n);
  std::vector<std::uint32_t> low_to_column(n, invalid);
  std::vector<bool> is_birth(n, false);
  std::vector<bool> is_killed(n, false);

  PersistenceDiagram diagram;

  for (std::uint32_t column_index = 0; column_index < n; ++column_index) {
    const SimplexId sigma = order[column_index];
    ModPColumn column = oriented_boundary_column_modp(
        complex, sigma, order_index, column_index, modulus);

    while (!column.empty()) {
      const auto& pivot = column.back();
      const std::uint32_t reducer = low_to_column[pivot.row];
      if (reducer == invalid) {
        break;
      }
      const std::uint32_t reducer_pivot = reduced_columns[reducer].back().coefficient;
      const std::uint32_t scale =
          (modulus - modp_multiply(pivot.coefficient,
                                   modp_inverse(reducer_pivot, modulus),
                                   modulus)) %
          modulus;
      modp_add_scaled_in_place(column, reduced_columns[reducer], scale, modulus);
    }

    if (column.empty()) {
      is_birth[column_index] = true;
      continue;
    }

    const std::uint32_t low = column.back().row;
    low_to_column[low] = column_index;
    reduced_columns[column_index] = std::move(column);
    is_killed[low] = true;

    const SimplexId birth = order[low];
    diagram.finite_pairs.push_back(PersistencePair{
        birth,
        sigma,
        complex.dimension(birth),
        complex.filtration(birth),
        complex.filtration(sigma),
    });
  }

  for (std::uint32_t column_index = 0; column_index < n; ++column_index) {
    if (!is_birth[column_index] || is_killed[column_index]) {
      continue;
    }

    const SimplexId birth = order[column_index];
    diagram.essential.push_back(EssentialInterval{
        birth,
        complex.dimension(birth),
        complex.filtration(birth),
    });
  }

  return diagram;
}

}  // namespace morseframes
