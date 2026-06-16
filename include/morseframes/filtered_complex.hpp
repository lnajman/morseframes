#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace morseframes {

using SimplexId = std::uint32_t;
using VertexId = std::uint32_t;
using LevelId = std::uint32_t;

constexpr SimplexId kInvalidSimplex = std::numeric_limits<SimplexId>::max();

struct Simplex {
  std::vector<VertexId> vertices;
  std::uint16_t dimension = 0;
  LevelId level = 0;
  double filtration = 0.0;
  std::vector<SimplexId> boundary;
  std::vector<SimplexId> coboundary;
};

class FilteredSimplicialComplex {
 public:
  void add_simplex(std::vector<VertexId> vertices, double filtration) {
    canonicalize(vertices);
    if (vertices.empty()) {
      throw std::invalid_argument("A simplex must contain at least one vertex.");
    }

    auto [it, inserted] = pending_.emplace(vertices, filtration);
    if (!inserted && std::fabs(it->second - filtration) > 1e-12) {
      throw std::invalid_argument("Duplicate simplex inserted with a different filtration value.");
    }
  }

  void finalize() {
    if (pending_.empty()) {
      throw std::invalid_argument("Cannot finalize an empty complex.");
    }

    simplices_.clear();
    simplex_to_id_.clear();
    level_values_.clear();
    level_buckets_.clear();
    filtration_order_.clear();

    for (const auto& [vertices, filtration] : pending_) {
      (void)filtration;
      const SimplexId id = checked_id(simplices_.size());
      simplex_to_id_.emplace(vertices, id);
      simplices_.push_back(Simplex{});
      simplices_.back().vertices = vertices;
      simplices_.back().dimension = checked_dimension(vertices.size() - 1);
      simplices_.back().filtration = pending_.at(vertices);
    }

    build_levels();
    build_boundaries_and_check_filtration();
    build_coboundaries();
    build_orders_and_buckets();
    finalized_ = true;
  }

  std::size_t size() const { return simplices_.size(); }

  const Simplex& simplex(SimplexId simplex) const {
    check_simplex_id(simplex);
    return simplices_[simplex];
  }

  std::uint16_t dimension(SimplexId simplex) const { return this->simplex(simplex).dimension; }
  LevelId level(SimplexId simplex) const { return this->simplex(simplex).level; }
  double filtration(SimplexId simplex) const { return this->simplex(simplex).filtration; }

  const std::vector<VertexId>& vertices(SimplexId simplex) const {
    return this->simplex(simplex).vertices;
  }

  const std::vector<SimplexId>& boundary(SimplexId simplex) const {
    return this->simplex(simplex).boundary;
  }

  const std::vector<SimplexId>& coboundary(SimplexId simplex) const {
    return this->simplex(simplex).coboundary;
  }

  const std::vector<SimplexId>& filtration_order() const { return filtration_order_; }

  const std::vector<SimplexId>& simplices_of_level(LevelId level) const {
    if (level >= level_buckets_.size()) {
      throw std::out_of_range("Invalid filtration level.");
    }
    return level_buckets_[level];
  }

  const std::vector<double>& level_values() const { return level_values_; }
  std::size_t num_levels() const { return level_values_.size(); }

  SimplexId find_simplex(const std::vector<VertexId>& vertices) const {
    std::vector<VertexId> canonical = vertices;
    canonicalize(canonical);
    auto it = simplex_to_id_.find(canonical);
    if (it == simplex_to_id_.end()) {
      return kInvalidSimplex;
    }
    return it->second;
  }

 private:
  struct VectorLess {
    bool operator()(const std::vector<VertexId>& lhs, const std::vector<VertexId>& rhs) const {
      return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
    }
  };

  static void canonicalize(std::vector<VertexId>& vertices) {
    std::sort(vertices.begin(), vertices.end());
    auto last = std::unique(vertices.begin(), vertices.end());
    if (last != vertices.end()) {
      throw std::invalid_argument("A simplex cannot contain duplicate vertices.");
    }
  }

  static SimplexId checked_id(std::size_t value) {
    if (value > std::numeric_limits<SimplexId>::max()) {
      throw std::overflow_error("Too many simplices for SimplexId.");
    }
    return static_cast<SimplexId>(value);
  }

  static std::uint16_t checked_dimension(std::size_t value) {
    if (value > std::numeric_limits<std::uint16_t>::max()) {
      throw std::overflow_error("Simplex dimension is too large.");
    }
    return static_cast<std::uint16_t>(value);
  }

  void check_simplex_id(SimplexId simplex) const {
    if (simplex >= simplices_.size()) {
      throw std::out_of_range("Invalid simplex id.");
    }
  }

  void build_levels() {
    for (const auto& simplex : simplices_) {
      level_values_.push_back(simplex.filtration);
    }
    std::sort(level_values_.begin(), level_values_.end());
    level_values_.erase(std::unique(level_values_.begin(), level_values_.end()), level_values_.end());

    for (auto& simplex : simplices_) {
      auto it = std::lower_bound(level_values_.begin(), level_values_.end(), simplex.filtration);
      simplex.level = static_cast<LevelId>(std::distance(level_values_.begin(), it));
    }
  }

  void build_boundaries_and_check_filtration() {
    for (SimplexId simplex_id = 0; simplex_id < simplices_.size(); ++simplex_id) {
      auto& simplex = simplices_[simplex_id];
      simplex.boundary.clear();

      if (simplex.vertices.size() == 1) {
        continue;
      }

      for (std::size_t removed = 0; removed < simplex.vertices.size(); ++removed) {
        std::vector<VertexId> face_vertices;
        face_vertices.reserve(simplex.vertices.size() - 1);
        for (std::size_t i = 0; i < simplex.vertices.size(); ++i) {
          if (i != removed) {
            face_vertices.push_back(simplex.vertices[i]);
          }
        }

        auto face_it = simplex_to_id_.find(face_vertices);
        if (face_it == simplex_to_id_.end()) {
          throw std::invalid_argument("Input is not closed under faces.");
        }

        const SimplexId face_id = face_it->second;
        if (simplices_[face_id].filtration > simplex.filtration + 1e-12) {
          throw std::invalid_argument("Filtration is not monotone on faces.");
        }
        simplex.boundary.push_back(face_id);
      }
    }
  }

  void build_coboundaries() {
    for (auto& simplex : simplices_) {
      simplex.coboundary.clear();
    }

    for (SimplexId simplex_id = 0; simplex_id < simplices_.size(); ++simplex_id) {
      for (SimplexId face : simplices_[simplex_id].boundary) {
        simplices_[face].coboundary.push_back(simplex_id);
      }
    }
  }

  void build_orders_and_buckets() {
    filtration_order_.resize(simplices_.size());
    for (SimplexId simplex_id = 0; simplex_id < simplices_.size(); ++simplex_id) {
      filtration_order_[simplex_id] = simplex_id;
    }

    auto simplex_less = [this](SimplexId lhs, SimplexId rhs) {
      const auto& a = simplices_[lhs];
      const auto& b = simplices_[rhs];
      if (a.level != b.level) {
        return a.level < b.level;
      }
      if (a.dimension != b.dimension) {
        return a.dimension < b.dimension;
      }
      return a.vertices < b.vertices;
    };
    std::sort(filtration_order_.begin(), filtration_order_.end(), simplex_less);

    level_buckets_.assign(level_values_.size(), {});
    for (SimplexId simplex_id : filtration_order_) {
      level_buckets_[simplices_[simplex_id].level].push_back(simplex_id);
    }
  }

  bool finalized_ = false;
  std::map<std::vector<VertexId>, double, VectorLess> pending_;
  std::map<std::vector<VertexId>, SimplexId, VectorLess> simplex_to_id_;
  std::vector<Simplex> simplices_;
  std::vector<double> level_values_;
  std::vector<std::vector<SimplexId>> level_buckets_;
  std::vector<SimplexId> filtration_order_;
};

}  // namespace morseframes
