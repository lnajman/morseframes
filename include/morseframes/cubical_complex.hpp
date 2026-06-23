#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "morseframes/filtered_complex.hpp"

namespace morseframes {

template <class Id, std::size_t Capacity>
class FixedIdList {
 public:
  using value_type = Id;
  using const_iterator = const Id*;
  using iterator = Id*;

  void clear() { size_ = 0; }

  void push_back(Id value) {
    if (size_ == Capacity) {
      throw std::overflow_error("FixedIdList capacity exceeded.");
    }
    values_[size_++] = value;
  }

  bool empty() const { return size_ == 0; }
  std::size_t size() const { return size_; }

  const Id& operator[](std::size_t index) const {
    if (index >= size_) {
      throw std::out_of_range("FixedIdList index out of range.");
    }
    return values_[index];
  }

  Id& operator[](std::size_t index) {
    if (index >= size_) {
      throw std::out_of_range("FixedIdList index out of range.");
    }
    return values_[index];
  }

  const_iterator begin() const { return values_.data(); }
  const_iterator end() const { return values_.data() + size_; }
  iterator begin() { return values_.data(); }
  iterator end() { return values_.data() + size_; }

  operator std::vector<Id>() const {
    return std::vector<Id>(begin(), end());
  }

 private:
  std::array<Id, Capacity> values_{};
  std::size_t size_ = 0;
};

template <class Id, std::size_t Capacity>
bool operator<(const FixedIdList<Id, Capacity>& lhs,
               const FixedIdList<Id, Capacity>& rhs) {
  return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

template <class Id, std::size_t Capacity>
bool operator==(const FixedIdList<Id, Capacity>& lhs,
                const std::vector<Id>& rhs) {
  return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template <class Id, std::size_t Capacity>
bool operator==(const std::vector<Id>& lhs,
                const FixedIdList<Id, Capacity>& rhs) {
  return rhs == lhs;
}

enum class CubicalCellType : std::uint8_t {
  Vertex,
  HorizontalEdge,
  VerticalEdge,
  Square,
};

class CubicalGrid2DComplex {
 public:
  using VertexList = FixedIdList<VertexId, 4>;
  using BoundaryList = FixedIdList<SimplexId, 4>;

  CubicalGrid2DComplex(std::size_t vertex_width,
                       std::size_t vertex_height,
                       std::vector<double> vertex_values)
      : vertex_width_(vertex_width),
        vertex_height_(vertex_height),
        vertex_values_(std::move(vertex_values)) {
    if (vertex_width_ == 0 || vertex_height_ == 0) {
      throw std::invalid_argument("A cubical grid must contain at least one vertex.");
    }
    if (vertex_values_.size() != checked_product(vertex_width_, vertex_height_)) {
      throw std::invalid_argument("Cubical grid vertex-value array has the wrong size.");
    }
    build();
  }

  static CubicalGrid2DComplex from_vertex_values(std::size_t vertex_width,
                                                 std::size_t vertex_height,
                                                 std::vector<double> vertex_values) {
    return CubicalGrid2DComplex(vertex_width, vertex_height, std::move(vertex_values));
  }

  std::size_t vertex_width() const { return vertex_width_; }
  std::size_t vertex_height() const { return vertex_height_; }
  std::size_t square_width() const { return vertex_width_ - 1; }
  std::size_t square_height() const { return vertex_height_ - 1; }

  std::size_t size() const { return cells_.size(); }

  CubicalCellType cell_type(SimplexId cell) const { return record(cell).type; }
  std::uint16_t dimension(SimplexId cell) const { return record(cell).dimension; }
  LevelId level(SimplexId cell) const { return record(cell).level; }
  double filtration(SimplexId cell) const { return record(cell).filtration; }

  const VertexList& vertices(SimplexId cell) const { return record(cell).vertices; }
  const BoundaryList& boundary(SimplexId cell) const { return record(cell).boundary; }
  const BoundaryList& coboundary(SimplexId cell) const { return record(cell).coboundary; }

  const std::vector<SimplexId>& filtration_order() const { return filtration_order_; }

  const std::vector<SimplexId>& simplices_of_level(LevelId level) const {
    if (level >= level_buckets_.size()) {
      throw std::out_of_range("Invalid filtration level.");
    }
    return level_buckets_[level];
  }

  const std::vector<double>& level_values() const { return level_values_; }
  std::size_t num_levels() const { return level_values_.size(); }

  std::uint32_t boundary_coefficient(SimplexId cell,
                                     std::size_t boundary_index,
                                     std::uint32_t modulus) const {
    if (modulus < 2) {
      throw std::invalid_argument("Coefficient modulus must be at least 2.");
    }
    const auto& data = record(cell);
    if (boundary_index >= data.boundary.size()) {
      throw std::out_of_range("Boundary coefficient index out of range.");
    }
    return data.boundary_signs[boundary_index] > 0 ? 1u : modulus - 1u;
  }

  SimplexId vertex(std::size_t x, std::size_t y) const {
    check_vertex_coordinates(x, y);
    return checked_simplex_id(y * vertex_width_ + x);
  }

  SimplexId horizontal_edge(std::size_t x, std::size_t y) const {
    if (vertex_width_ < 2 || x + 1 >= vertex_width_ || y >= vertex_height_) {
      throw std::out_of_range("Invalid horizontal edge coordinates.");
    }
    return checked_simplex_id(num_vertices() + y * square_width() + x);
  }

  SimplexId vertical_edge(std::size_t x, std::size_t y) const {
    if (vertex_height_ < 2 || x >= vertex_width_ || y + 1 >= vertex_height_) {
      throw std::out_of_range("Invalid vertical edge coordinates.");
    }
    return checked_simplex_id(num_vertices() + num_horizontal_edges() + y * vertex_width_ + x);
  }

  SimplexId square(std::size_t x, std::size_t y) const {
    if (vertex_width_ < 2 || vertex_height_ < 2 ||
        x + 1 >= vertex_width_ || y + 1 >= vertex_height_) {
      throw std::out_of_range("Invalid square coordinates.");
    }
    return checked_simplex_id(
        num_vertices() + num_horizontal_edges() + num_vertical_edges() +
        y * square_width() + x);
  }

 private:
  struct CellRecord {
    CubicalCellType type = CubicalCellType::Vertex;
    VertexList vertices;
    std::uint16_t dimension = 0;
    LevelId level = 0;
    double filtration = 0.0;
    BoundaryList boundary;
    BoundaryList coboundary;
    std::array<std::int8_t, 4> boundary_signs{};
  };

  static std::size_t checked_product(std::size_t lhs, std::size_t rhs) {
    if (lhs != 0 && rhs > std::numeric_limits<std::size_t>::max() / lhs) {
      throw std::overflow_error("Cubical grid is too large.");
    }
    return lhs * rhs;
  }

  static std::size_t checked_sum(std::size_t lhs, std::size_t rhs) {
    if (rhs > std::numeric_limits<std::size_t>::max() - lhs) {
      throw std::overflow_error("Cubical grid is too large.");
    }
    return lhs + rhs;
  }

  static SimplexId checked_simplex_id(std::size_t value) {
    if (value > std::numeric_limits<SimplexId>::max()) {
      throw std::overflow_error("Too many cubical cells for SimplexId.");
    }
    return static_cast<SimplexId>(value);
  }

  std::size_t num_vertices() const {
    return checked_product(vertex_width_, vertex_height_);
  }

  std::size_t num_horizontal_edges() const {
    return vertex_width_ < 2 ? 0 : checked_product(vertex_width_ - 1, vertex_height_);
  }

  std::size_t num_vertical_edges() const {
    return vertex_height_ < 2 ? 0 : checked_product(vertex_width_, vertex_height_ - 1);
  }

  std::size_t num_squares() const {
    return vertex_width_ < 2 || vertex_height_ < 2
               ? 0
               : checked_product(vertex_width_ - 1, vertex_height_ - 1);
  }

  void check_vertex_coordinates(std::size_t x, std::size_t y) const {
    if (x >= vertex_width_ || y >= vertex_height_) {
      throw std::out_of_range("Invalid vertex coordinates.");
    }
  }

  const CellRecord& record(SimplexId cell) const {
    if (cell >= cells_.size()) {
      throw std::out_of_range("Invalid cubical cell id.");
    }
    return cells_[cell];
  }

  CellRecord& record(SimplexId cell) {
    if (cell >= cells_.size()) {
      throw std::out_of_range("Invalid cubical cell id.");
    }
    return cells_[cell];
  }

  double vertex_value(std::size_t x, std::size_t y) const {
    check_vertex_coordinates(x, y);
    return vertex_values_.at(y * vertex_width_ + x);
  }

  double max_vertex_value(std::initializer_list<VertexId> vertices) const {
    if (vertices.size() == 0) {
      throw std::invalid_argument("A cubical cell must have at least one vertex.");
    }
    double value = -std::numeric_limits<double>::infinity();
    for (VertexId vertex_id : vertices) {
      value = std::max(value, vertex_values_.at(vertex_id));
    }
    return value;
  }

  void assign_cell(SimplexId cell,
                   CubicalCellType type,
                   std::uint16_t dimension,
                   std::initializer_list<VertexId> vertices) {
    auto& data = record(cell);
    data.type = type;
    data.dimension = dimension;
    data.vertices.clear();
    data.boundary.clear();
    data.coboundary.clear();
    data.boundary_signs.fill(0);
    for (VertexId vertex_id : vertices) {
      data.vertices.push_back(vertex_id);
    }
    data.filtration = max_vertex_value(vertices);
  }

  void add_boundary(SimplexId cell, SimplexId face, std::int8_t sign) {
    if (sign != 1 && sign != -1) {
      throw std::invalid_argument("Cubical incidence signs must be +1 or -1.");
    }
    auto& data = record(cell);
    if (record(face).filtration > data.filtration) {
      throw std::logic_error("Cubical filtration is not monotone on faces.");
    }
    const std::size_t index = data.boundary.size();
    data.boundary.push_back(face);
    data.boundary_signs[index] = sign;
    record(face).coboundary.push_back(cell);
  }

  VertexId grid_vertex_id(std::size_t x, std::size_t y) const {
    return static_cast<VertexId>(vertex(x, y));
  }

  void build() {
    const std::size_t total = checked_sum(
        checked_sum(num_vertices(), num_horizontal_edges()),
        checked_sum(num_vertical_edges(), num_squares()));
    checked_simplex_id(total);
    cells_.assign(total, CellRecord{});

    for (std::size_t y = 0; y < vertex_height_; ++y) {
      for (std::size_t x = 0; x < vertex_width_; ++x) {
        assign_cell(vertex(x, y),
                    CubicalCellType::Vertex,
                    0,
                    {grid_vertex_id(x, y)});
      }
    }

    for (std::size_t y = 0; y < vertex_height_; ++y) {
      for (std::size_t x = 0; x + 1 < vertex_width_; ++x) {
        const SimplexId edge = horizontal_edge(x, y);
        assign_cell(edge,
                    CubicalCellType::HorizontalEdge,
                    1,
                    {grid_vertex_id(x, y), grid_vertex_id(x + 1, y)});
        add_boundary(edge, vertex(x + 1, y), 1);
        add_boundary(edge, vertex(x, y), -1);
      }
    }

    for (std::size_t y = 0; y + 1 < vertex_height_; ++y) {
      for (std::size_t x = 0; x < vertex_width_; ++x) {
        const SimplexId edge = vertical_edge(x, y);
        assign_cell(edge,
                    CubicalCellType::VerticalEdge,
                    1,
                    {grid_vertex_id(x, y), grid_vertex_id(x, y + 1)});
        add_boundary(edge, vertex(x, y + 1), 1);
        add_boundary(edge, vertex(x, y), -1);
      }
    }

    for (std::size_t y = 0; y + 1 < vertex_height_; ++y) {
      for (std::size_t x = 0; x + 1 < vertex_width_; ++x) {
        const SimplexId face = square(x, y);
        assign_cell(face,
                    CubicalCellType::Square,
                    2,
                    {grid_vertex_id(x, y),
                     grid_vertex_id(x + 1, y),
                     grid_vertex_id(x, y + 1),
                     grid_vertex_id(x + 1, y + 1)});
        add_boundary(face, vertical_edge(x + 1, y), 1);
        add_boundary(face, vertical_edge(x, y), -1);
        add_boundary(face, horizontal_edge(x, y + 1), -1);
        add_boundary(face, horizontal_edge(x, y), 1);
      }
    }

    build_levels();
    build_orders_and_buckets();
  }

  void build_levels() {
    level_values_.clear();
    level_values_.reserve(cells_.size());
    for (const auto& cell : cells_) {
      level_values_.push_back(cell.filtration);
    }
    std::sort(level_values_.begin(), level_values_.end());
    level_values_.erase(std::unique(level_values_.begin(), level_values_.end()),
                        level_values_.end());

    for (auto& cell : cells_) {
      auto it = std::lower_bound(level_values_.begin(), level_values_.end(), cell.filtration);
      cell.level = static_cast<LevelId>(std::distance(level_values_.begin(), it));
    }
  }

  void build_orders_and_buckets() {
    filtration_order_.resize(cells_.size());
    for (SimplexId cell = 0; cell < cells_.size(); ++cell) {
      filtration_order_[cell] = cell;
    }

    auto cell_less = [this](SimplexId lhs, SimplexId rhs) {
      const auto& a = cells_[lhs];
      const auto& b = cells_[rhs];
      if (a.level != b.level) {
        return a.level < b.level;
      }
      if (a.dimension != b.dimension) {
        return a.dimension < b.dimension;
      }
      if (a.vertices < b.vertices) {
        return true;
      }
      if (b.vertices < a.vertices) {
        return false;
      }
      return lhs < rhs;
    };
    std::sort(filtration_order_.begin(), filtration_order_.end(), cell_less);

    level_buckets_.assign(level_values_.size(), {});
    for (SimplexId cell : filtration_order_) {
      level_buckets_[cells_[cell].level].push_back(cell);
    }
  }

  std::size_t vertex_width_ = 0;
  std::size_t vertex_height_ = 0;
  std::vector<double> vertex_values_;
  std::vector<CellRecord> cells_;
  std::vector<double> level_values_;
  std::vector<std::vector<SimplexId>> level_buckets_;
  std::vector<SimplexId> filtration_order_;
};

}  // namespace morseframes
