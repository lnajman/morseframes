#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gudhi/Morse_persistence/internal/filtered_complex.h"

namespace Gudhi { namespace morse_persistence { namespace internal {

template <class Id, std::size_t InlineCapacity>
class SmallIdList {
 public:
  using value_type = Id;
  using const_iterator = const Id*;
  using iterator = Id*;

  SmallIdList() = default;
  ~SmallIdList() { delete[] overflow_; }

  SmallIdList(const SmallIdList& other) {
    copy_from(other);
  }

  SmallIdList& operator=(const SmallIdList& other) {
    if (this != &other) {
      SmallIdList copy(other);
      swap(copy);
    }
    return *this;
  }

  SmallIdList(SmallIdList&& other) noexcept {
    move_from(std::move(other));
  }

  SmallIdList& operator=(SmallIdList&& other) noexcept {
    if (this != &other) {
      delete[] overflow_;
      overflow_ = nullptr;
      size_ = 0;
      capacity_ = InlineCapacity;
      move_from(std::move(other));
    }
    return *this;
  }

  void clear() {
    size_ = 0;
  }

  void reserve(std::size_t capacity) {
    if (capacity <= capacity_) {
      return;
    }
    grow_to(capacity);
  }

  void push_back(Id value) {
    if (size_ == capacity_) {
      grow_to(std::max<std::size_t>(capacity_ * 2, InlineCapacity + 1));
    }
    data()[size_++] = value;
  }

  bool empty() const { return size() == 0; }

  std::size_t size() const {
    return size_;
  }

  const Id& operator[](std::size_t index) const { return data()[index]; }
  Id& operator[](std::size_t index) { return data()[index]; }

  const Id& front() const { return (*this)[0]; }
  Id& front() { return (*this)[0]; }

  const_iterator begin() const { return data(); }
  const_iterator end() const { return data() + size(); }
  iterator begin() { return data(); }
  iterator end() { return data() + size(); }

  operator std::vector<Id>() const {
    return std::vector<Id>(begin(), end());
  }

 private:
  void grow_to(std::size_t new_capacity) {
    Id* replacement = new Id[new_capacity];
    std::copy(begin(), end(), replacement);
    delete[] overflow_;
    overflow_ = replacement;
    capacity_ = new_capacity;
  }

  const Id* data() const {
    return overflow_ == nullptr ? inline_.data() : overflow_;
  }

  Id* data() {
    return overflow_ == nullptr ? inline_.data() : overflow_;
  }

  void copy_from(const SmallIdList& other) {
    size_ = other.size_;
    capacity_ = other.size_ > InlineCapacity ? other.size_ : InlineCapacity;
    if (capacity_ > InlineCapacity) {
      overflow_ = new Id[capacity_];
    }
    std::copy(other.begin(), other.end(), data());
  }

  void move_from(SmallIdList&& other) {
    size_ = other.size_;
    if (other.overflow_ != nullptr) {
      overflow_ = other.overflow_;
      capacity_ = other.capacity_;
      other.overflow_ = nullptr;
      other.size_ = 0;
      other.capacity_ = InlineCapacity;
      return;
    }
    capacity_ = InlineCapacity;
    std::copy(other.inline_.begin(), other.inline_.begin() + other.size_, inline_.begin());
    other.size_ = 0;
  }

  void swap(SmallIdList& other) noexcept {
    using std::swap;
    swap(inline_, other.inline_);
    swap(overflow_, other.overflow_);
    swap(size_, other.size_);
    swap(capacity_, other.capacity_);
  }

  std::array<Id, InlineCapacity> inline_{};
  Id* overflow_ = nullptr;
  std::size_t size_ = 0;
  std::size_t capacity_ = InlineCapacity;
};

template <class Id, std::size_t InlineCapacity>
bool operator<(const SmallIdList<Id, InlineCapacity>& lhs,
               const SmallIdList<Id, InlineCapacity>& rhs) {
  return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

template <class Id, std::size_t InlineCapacity>
bool operator==(const SmallIdList<Id, InlineCapacity>& lhs,
                const std::vector<Id>& rhs) {
  return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template <class Id, std::size_t InlineCapacity>
bool operator==(const std::vector<Id>& lhs,
                const SmallIdList<Id, InlineCapacity>& rhs) {
  return rhs == lhs;
}

template <std::size_t InlineCapacity>
using SmallSimplexIdList = SmallIdList<SimplexId, InlineCapacity>;

template <std::size_t InlineCapacity>
using SmallVertexIdList = SmallIdList<VertexId, InlineCapacity>;

enum class DuplicateFiltrationPolicy {
  Min,
  Strict,
};

inline DuplicateFiltrationPolicy duplicate_policy_from_string(const std::string& policy) {
  if (policy == "min") {
    return DuplicateFiltrationPolicy::Min;
  }
  if (policy == "strict") {
    return DuplicateFiltrationPolicy::Strict;
  }
  throw std::invalid_argument("Duplicate policy must be either 'min' or 'strict'.");
}

class SimplexTreeBuilder {
 public:
  explicit SimplexTreeBuilder(
      DuplicateFiltrationPolicy duplicate_policy = DuplicateFiltrationPolicy::Min)
      : duplicate_policy_(duplicate_policy) {}

  explicit SimplexTreeBuilder(const std::string& duplicate_policy)
      : duplicate_policy_(duplicate_policy_from_string(duplicate_policy)) {}

  bool insert(std::vector<VertexId> simplex, double filtration = 0.0) {
    canonicalize(simplex);
    bool changed = false;
    for_each_nonempty_face(simplex, [&](const std::vector<VertexId>& face) {
      changed = insert_canonical_simplex(face, filtration) || changed;
    });
    return changed;
  }

  bool insert_simplex_only(std::vector<VertexId> simplex, double filtration = 0.0) {
    canonicalize(simplex);
    return insert_canonical_simplex(simplex, filtration);
  }

  bool add_simplex(std::vector<VertexId> simplex, double filtration = 0.0) {
    return insert(std::move(simplex), filtration);
  }

  bool contains(std::vector<VertexId> simplex) const {
    canonicalize(simplex);
    const Node* node = find_node(simplex);
    return node != nullptr && node->has_filtration;
  }

  bool find_simplex(std::vector<VertexId> simplex) const {
    return contains(std::move(simplex));
  }

  double filtration(std::vector<VertexId> simplex) const {
    canonicalize(simplex);
    const Node* node = find_node(simplex);
    if (node == nullptr || !node->has_filtration) {
      throw std::out_of_range("Simplex is not present in the builder.");
    }
    return node->filtration;
  }

  double simplex_filtration(std::vector<VertexId> simplex) const {
    return filtration(std::move(simplex));
  }

  std::size_t size() const { return size_; }
  std::size_t num_simplices() const { return size_; }
  std::size_t num_vertices() const { return vertices_.size(); }
  int max_dimension() const { return max_dimension_; }
  bool empty() const { return size_ == 0; }

  std::vector<std::pair<std::vector<VertexId>, double>> simplices() const {
    std::vector<std::pair<std::vector<VertexId>, double>> result;
    result.reserve(size_);
    std::vector<VertexId> prefix;
    collect_simplices(root_, prefix, result);
    return result;
  }

  std::vector<std::pair<std::vector<VertexId>, double>> get_filtration() const {
    auto result = simplices();
    std::stable_sort(
        result.begin(),
        result.end(),
        [](const auto& lhs, const auto& rhs) {
          if (lhs.second != rhs.second) {
            return lhs.second < rhs.second;
          }
          if (lhs.first.size() != rhs.first.size()) {
            return lhs.first.size() < rhs.first.size();
          }
          return lhs.first < rhs.first;
        });
    return result;
  }

  FilteredSimplicialComplex to_filtered_complex(bool finalize = true) const {
    FilteredSimplicialComplex complex;
    for (auto [simplex, filtration_value] : simplices()) {
      complex.add_simplex(std::move(simplex), filtration_value);
    }
    if (finalize) {
      complex.finalize();
    }
    return complex;
  }

  FilteredSimplicialComplex finalize(bool clear_builder = true) {
    FilteredSimplicialComplex complex = to_filtered_complex(true);
    if (clear_builder) {
      clear();
    }
    return complex;
  }

  void clear() {
    root_ = Node{};
    vertices_.clear();
    size_ = 0;
    max_dimension_ = -1;
  }

 private:
  struct Node {
    std::map<VertexId, Node> children;
    bool has_filtration = false;
    double filtration = 0.0;
  };

  static void canonicalize(std::vector<VertexId>& simplex) {
    std::sort(simplex.begin(), simplex.end());
    auto last = std::unique(simplex.begin(), simplex.end());
    if (last != simplex.end()) {
      throw std::invalid_argument("A simplex cannot contain duplicate vertices.");
    }
    if (simplex.empty()) {
      throw std::invalid_argument("A simplex must contain at least one vertex.");
    }
  }

  template <class Callback>
  static void for_each_nonempty_face(const std::vector<VertexId>& simplex, Callback callback) {
    std::vector<VertexId> face;
    face.reserve(simplex.size());

    std::function<void(std::size_t)> visit = [&](std::size_t index) {
      if (index == simplex.size()) {
        if (!face.empty()) {
          callback(face);
        }
        return;
      }

      face.push_back(simplex[index]);
      visit(index + 1);
      face.pop_back();
      visit(index + 1);
    };
    visit(0);
  }

  bool insert_canonical_simplex(const std::vector<VertexId>& simplex, double filtration) {
    Node* node = ensure_node(simplex);
    for (VertexId vertex : simplex) {
      vertices_.insert(vertex);
    }

    if (node->has_filtration) {
      if (duplicate_policy_ == DuplicateFiltrationPolicy::Strict &&
          std::fabs(node->filtration - filtration) > 1e-12) {
        throw std::invalid_argument(
            "Duplicate simplex inserted with a different filtration value.");
      }
      const double old_value = node->filtration;
      node->filtration = std::min(old_value, filtration);
      return std::fabs(old_value - node->filtration) > 1e-12;
    }

    node->has_filtration = true;
    node->filtration = filtration;
    ++size_;
    const auto dimension = static_cast<int>(simplex.size()) - 1;
    if (dimension > max_dimension_) {
      max_dimension_ = dimension;
    }
    return true;
  }

  Node* ensure_node(const std::vector<VertexId>& simplex) {
    Node* node = &root_;
    for (VertexId vertex : simplex) {
      node = &node->children[vertex];
    }
    return node;
  }

  const Node* find_node(const std::vector<VertexId>& simplex) const {
    const Node* node = &root_;
    for (VertexId vertex : simplex) {
      auto child = node->children.find(vertex);
      if (child == node->children.end()) {
        return nullptr;
      }
      node = &child->second;
    }
    return node;
  }

  static void collect_simplices(
      const Node& node,
      std::vector<VertexId>& prefix,
      std::vector<std::pair<std::vector<VertexId>, double>>& result) {
    if (node.has_filtration) {
      result.emplace_back(prefix, node.filtration);
    }

    for (const auto& [vertex, child] : node.children) {
      prefix.push_back(vertex);
      collect_simplices(child, prefix, result);
      prefix.pop_back();
    }
  }

  DuplicateFiltrationPolicy duplicate_policy_ = DuplicateFiltrationPolicy::Min;
  Node root_;
  std::set<VertexId> vertices_;
  std::size_t size_ = 0;
  int max_dimension_ = -1;
};

template <class SimplexTree>
class SimplexTreeComplexView {
 public:
  using SimplexHandle = std::decay_t<decltype(*std::begin(
      std::declval<const SimplexTree&>().filtration_simplex_range()))>;

  struct BuildMetrics {
    std::uint64_t extract_nanoseconds = 0;
    std::uint64_t boundary_nanoseconds = 0;
    std::uint64_t boundary_scan_nanoseconds = 0;
    std::uint64_t boundary_lookup_setup_nanoseconds = 0;
    std::uint64_t boundary_register_nanoseconds = 0;
    std::uint64_t boundary_face_lookup_nanoseconds = 0;
    std::uint64_t coboundary_nanoseconds = 0;
    std::uint64_t coboundary_count_nanoseconds = 0;
    std::uint64_t coboundary_reserve_nanoseconds = 0;
    std::uint64_t coboundary_fill_nanoseconds = 0;
    std::uint64_t order_nanoseconds = 0;
    std::uint64_t order_size_nanoseconds = 0;
    std::uint64_t order_bucket_nanoseconds = 0;
    std::uint64_t order_sort_emit_nanoseconds = 0;
    std::uint64_t total_nanoseconds = 0;
    std::size_t boundary_vertex_count = 0;
    std::size_t boundary_edge_count = 0;
    std::size_t boundary_dense_vertex_size = 0;
    std::size_t boundary_dense_edge_size = 0;
    std::size_t boundary_vertex_face_lookups = 0;
    std::size_t boundary_edge_face_lookups = 0;
    std::size_t boundary_general_face_lookups = 0;
    bool boundary_used_dense_vertex_lookup = false;
    bool boundary_used_dense_edge_lookup = false;
  };

  explicit SimplexTreeComplexView(const SimplexTree& simplex_tree)
      : simplex_tree_(&simplex_tree) {
    build();
  }

  SimplexTreeComplexView(const SimplexTreeComplexView&) = default;
  SimplexTreeComplexView& operator=(const SimplexTreeComplexView&) = default;
  SimplexTreeComplexView(SimplexTreeComplexView&&) noexcept = default;
  SimplexTreeComplexView& operator=(SimplexTreeComplexView&&) noexcept = default;

  std::size_t size() const { return simplices_.size(); }

  std::uint16_t dimension(SimplexId simplex) const { return record(simplex).dimension; }
  LevelId level(SimplexId simplex) const { return record(simplex).level; }
  double filtration(SimplexId simplex) const { return record(simplex).filtration; }

  const SmallVertexIdList<4>& vertices(SimplexId simplex) const {
    return record(simplex).vertices;
  }

  const SmallSimplexIdList<4>& boundary(SimplexId simplex) const {
    return record(simplex).boundary;
  }

  const SmallSimplexIdList<4>& coboundary(SimplexId simplex) const {
    return record(simplex).coboundary;
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

  SimplexId find_simplex(std::vector<VertexId> simplex) const {
    canonicalize(simplex);
    ensure_simplex_lookup();
    auto it = simplex_to_id_.find(simplex);
    if (it == simplex_to_id_.end()) {
      return kInvalidSimplex;
    }
    return it->second;
  }

  SimplexHandle handle(SimplexId simplex) const {
    return record(simplex).handle;
  }

  const BuildMetrics& build_metrics() const { return build_metrics_; }

 private:
  using Clock = std::chrono::steady_clock;

  struct VectorHash {
    std::size_t operator()(const std::vector<VertexId>& simplex) const noexcept {
      std::size_t seed = simplex.size();
      for (VertexId vertex : simplex) {
        seed ^= static_cast<std::size_t>(vertex) + 0x9e3779b97f4a7c15ULL +
                (seed << 6U) + (seed >> 2U);
      }
      return seed;
    }
  };

  class FlatEdgeLookup {
   public:
    void reserve(std::size_t size) {
      std::size_t capacity = 16;
      while (capacity < size * 4) {
        capacity *= 2;
      }
      keys_.assign(capacity, 0);
      values_.assign(capacity, kInvalidSimplex);
      mask_ = capacity - 1;
    }

    bool emplace(std::uint64_t key, SimplexId value) {
      ensure_initialized();
      std::size_t index = hash(key) & mask_;
      while (values_[index] != kInvalidSimplex) {
        if (keys_[index] == key) {
          return false;
        }
        index = (index + 1) & mask_;
      }
      keys_[index] = key;
      values_[index] = value;
      return true;
    }

    SimplexId find(std::uint64_t key) const {
      if (values_.empty()) {
        return kInvalidSimplex;
      }
      std::size_t index = hash(key) & mask_;
      while (values_[index] != kInvalidSimplex) {
        if (keys_[index] == key) {
          return values_[index];
        }
        index = (index + 1) & mask_;
      }
      return kInvalidSimplex;
    }

   private:
    void ensure_initialized() {
      if (values_.empty()) {
        reserve(1);
      }
    }

    static std::size_t hash(std::uint64_t value) {
      std::uint64_t mixed = value ^ (value >> 32U);
      mixed *= 0x9e3779b97f4a7c15ULL;
      mixed ^= mixed >> 32U;
      return static_cast<std::size_t>(mixed);
    }

    std::vector<std::uint64_t> keys_;
    std::vector<SimplexId> values_;
    std::size_t mask_ = 0;
  };

  struct Record {
    SimplexHandle handle;
    SmallVertexIdList<4> vertices;
    std::uint16_t dimension = 0;
    LevelId level = 0;
    double filtration = 0.0;
    SmallSimplexIdList<4> boundary;
    SmallSimplexIdList<4> coboundary;
  };

  static void canonicalize(std::vector<VertexId>& simplex) {
    std::sort(simplex.begin(), simplex.end());
    auto last = std::unique(simplex.begin(), simplex.end());
    if (last != simplex.end()) {
      throw std::invalid_argument("A simplex cannot contain duplicate vertices.");
    }
    if (simplex.empty()) {
      throw std::invalid_argument("A simplex must contain at least one vertex.");
    }
  }

  template <class VertexList>
  static void normalize_simplex_vertices(VertexList& simplex,
                                         bool strictly_increasing,
                                         bool strictly_decreasing) {
    if (simplex.empty()) {
      throw std::invalid_argument("A simplex must contain at least one vertex.");
    }
    if (strictly_increasing) {
      return;
    }
    if (strictly_decreasing) {
      std::reverse(simplex.begin(), simplex.end());
      return;
    }
    std::sort(simplex.begin(), simplex.end());
    auto last = std::unique(simplex.begin(), simplex.end());
    if (last != simplex.end()) {
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

  static std::uint64_t elapsed_nanoseconds(Clock::time_point start,
                                           Clock::time_point stop) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count());
  }

  static std::uint64_t packed_edge_key(VertexId lhs, VertexId rhs) {
    if (rhs < lhs) {
      std::swap(lhs, rhs);
    }
    return (static_cast<std::uint64_t>(lhs) << 32U) | static_cast<std::uint64_t>(rhs);
  }

  template <class Tree>
  static auto size_hint(const Tree& simplex_tree, int)
      -> decltype(static_cast<std::size_t>(simplex_tree.num_simplices())) {
    return static_cast<std::size_t>(simplex_tree.num_simplices());
  }

  template <class Tree>
  static std::size_t size_hint(const Tree&, long) {
    return 0;
  }

  const Record& record(SimplexId simplex) const {
    if (simplex >= simplices_.size()) {
      throw std::out_of_range("Invalid simplex id.");
    }
    return simplices_[simplex];
  }

  Record& record(SimplexId simplex) {
    if (simplex >= simplices_.size()) {
      throw std::out_of_range("Invalid simplex id.");
    }
    return simplices_[simplex];
  }

  bool vertices_less_in_same_dimension(SimplexId lhs,
                                       SimplexId rhs,
                                       std::size_t dimension) const {
    const auto& lhs_vertices = simplices_[lhs].vertices;
    const auto& rhs_vertices = simplices_[rhs].vertices;
    switch (dimension) {
      case 0:
        return lhs_vertices[0] < rhs_vertices[0];
      case 1:
        if (lhs_vertices[0] != rhs_vertices[0]) {
          return lhs_vertices[0] < rhs_vertices[0];
        }
        return lhs_vertices[1] < rhs_vertices[1];
      case 2:
        if (lhs_vertices[0] != rhs_vertices[0]) {
          return lhs_vertices[0] < rhs_vertices[0];
        }
        if (lhs_vertices[1] != rhs_vertices[1]) {
          return lhs_vertices[1] < rhs_vertices[1];
        }
        return lhs_vertices[2] < rhs_vertices[2];
      case 3:
        if (lhs_vertices[0] != rhs_vertices[0]) {
          return lhs_vertices[0] < rhs_vertices[0];
        }
        if (lhs_vertices[1] != rhs_vertices[1]) {
          return lhs_vertices[1] < rhs_vertices[1];
        }
        if (lhs_vertices[2] != rhs_vertices[2]) {
          return lhs_vertices[2] < rhs_vertices[2];
        }
        return lhs_vertices[3] < rhs_vertices[3];
      default:
        return std::lexicographical_compare(lhs_vertices.begin(),
                                            lhs_vertices.end(),
                                            rhs_vertices.begin(),
                                            rhs_vertices.end());
    }
  }

  void ensure_simplex_lookup() const {
    if (simplex_lookup_built_) {
      return;
    }
    simplex_to_id_.clear();
    simplex_to_id_.reserve(simplices_.size());
    for (SimplexId simplex_id = 0; simplex_id < simplices_.size(); ++simplex_id) {
      const auto& vertices = record(simplex_id).vertices;
      std::vector<VertexId> key(vertices.begin(), vertices.end());
      auto [_, inserted] = simplex_to_id_.emplace(std::move(key), simplex_id);
      if (!inserted) {
        throw std::invalid_argument("Simplex tree filtration range contains a duplicate simplex.");
      }
    }
    simplex_lookup_built_ = true;
  }

  void build() {
    build_metrics_ = BuildMetrics{};
    const auto total_start = Clock::now();

    const std::size_t hint = size_hint(*simplex_tree_, 0);
    if (hint > 0) {
      simplices_.reserve(hint);
      level_values_.reserve(hint);
    }

    std::vector<SimplexId> vertex_simplex_ids;
    std::vector<SimplexId> edge_simplex_ids;
    if (hint > 0) {
      vertex_simplex_ids.reserve(hint / 4);
      edge_simplex_ids.reserve(hint / 2);
    }
    VertexId max_vertex = 0;
    bool have_vertex_handle = false;
    bool have_previous_filtration = false;
    double previous_filtration = 0.0;
    LevelId current_level = 0;
    for (auto handle : simplex_tree_->filtration_simplex_range()) {
      SmallVertexIdList<4> vertices;
      bool strictly_increasing = true;
      bool strictly_decreasing = true;
      bool have_previous_vertex = false;
      VertexId previous_vertex = 0;
      for (auto vertex : simplex_tree_->simplex_vertex_range(handle)) {
        const auto current_vertex = static_cast<VertexId>(vertex);
        max_vertex = have_vertex_handle ? std::max(max_vertex, current_vertex) : current_vertex;
        have_vertex_handle = true;
        if (have_previous_vertex) {
          strictly_increasing = strictly_increasing && previous_vertex < current_vertex;
          strictly_decreasing = strictly_decreasing && previous_vertex > current_vertex;
        }
        vertices.push_back(current_vertex);
        previous_vertex = current_vertex;
        have_previous_vertex = true;
      }
      normalize_simplex_vertices(vertices, strictly_increasing, strictly_decreasing);

      const double filtration = simplex_tree_->filtration(handle);
      if (!have_previous_filtration || filtration != previous_filtration) {
        if (have_previous_filtration && filtration < previous_filtration - 1e-12) {
          throw std::invalid_argument(
              "Simplex tree filtration_simplex_range() must be nondecreasing.");
        }
        current_level = static_cast<LevelId>(level_values_.size());
        level_values_.push_back(filtration);
        previous_filtration = filtration;
        have_previous_filtration = true;
      }

      const SimplexId id = checked_id(simplices_.size());
      const auto dimension = checked_dimension(vertices.size() - 1);
      simplices_.push_back(Record{
          handle,
          std::move(vertices),
          dimension,
          current_level,
          filtration,
          {},
          {},
      });
      if (dimension == 0) {
        vertex_simplex_ids.push_back(id);
      } else if (dimension == 1) {
        edge_simplex_ids.push_back(id);
      }
    }
    const auto after_extract = Clock::now();
    build_metrics_.extract_nanoseconds =
        elapsed_nanoseconds(total_start, after_extract);

    if (simplices_.empty()) {
      throw std::invalid_argument("Cannot build a view of an empty simplex tree.");
    }

    build_boundaries_and_check_filtration(vertex_simplex_ids,
                                          edge_simplex_ids,
                                          max_vertex,
                                          have_vertex_handle);
    const auto after_boundary = Clock::now();
    build_metrics_.boundary_nanoseconds =
        elapsed_nanoseconds(after_extract, after_boundary);

    build_coboundaries();
    const auto after_coboundary = Clock::now();
    build_metrics_.coboundary_nanoseconds =
        elapsed_nanoseconds(after_boundary, after_coboundary);

    build_orders_and_buckets();
    const auto after_order = Clock::now();
    build_metrics_.order_nanoseconds =
        elapsed_nanoseconds(after_coboundary, after_order);
    build_metrics_.total_nanoseconds = elapsed_nanoseconds(total_start, after_order);
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

  void build_boundaries_and_check_filtration(
      const std::vector<SimplexId>& vertex_simplex_ids,
      const std::vector<SimplexId>& edge_simplex_ids,
      VertexId max_vertex,
      bool have_vertex_handle) {
    const auto scan_start = Clock::now();
    const std::size_t vertex_count = vertex_simplex_ids.size();
    const std::size_t edge_count = edge_simplex_ids.size();
    const auto after_scan = Clock::now();
    build_metrics_.boundary_scan_nanoseconds =
        elapsed_nanoseconds(scan_start, after_scan);
    build_metrics_.boundary_vertex_count = vertex_count;
    build_metrics_.boundary_edge_count = edge_count;

    const std::size_t dense_vertex_size =
        have_vertex_handle ? static_cast<std::size_t>(max_vertex) + 1 : 0;
    const bool use_dense_vertex_lookup =
        have_vertex_handle &&
        dense_vertex_size <= std::max<std::size_t>(4 * vertex_count, std::size_t{1024});
    const bool dense_edge_size_fits =
        dense_vertex_size > 0 &&
        dense_vertex_size <= std::numeric_limits<std::size_t>::max() / dense_vertex_size;
    const std::size_t dense_edge_size =
        dense_edge_size_fits ? dense_vertex_size * dense_vertex_size : 0;
    const bool use_dense_edge_lookup =
        use_dense_vertex_lookup &&
        dense_edge_size <= std::max<std::size_t>(16 * edge_count, std::size_t{100000}) &&
        dense_edge_size <= std::size_t{100000};
    build_metrics_.boundary_dense_vertex_size = dense_vertex_size;
    build_metrics_.boundary_dense_edge_size = use_dense_edge_lookup ? dense_edge_size : 0;
    build_metrics_.boundary_used_dense_vertex_lookup = use_dense_vertex_lookup;
    build_metrics_.boundary_used_dense_edge_lookup = use_dense_edge_lookup;
    std::vector<SimplexId> dense_vertex_to_id;
    std::vector<SimplexId> dense_edge_to_id;
    std::unordered_map<VertexId, SimplexId> vertex_to_id;
    FlatEdgeLookup edge_to_id;
    if (use_dense_vertex_lookup) {
      dense_vertex_to_id.assign(dense_vertex_size, kInvalidSimplex);
    } else {
      vertex_to_id.reserve(vertex_count);
    }
    if (use_dense_edge_lookup) {
      dense_edge_to_id.assign(dense_edge_size, kInvalidSimplex);
    } else {
      edge_to_id.reserve(edge_count);
    }
    const auto after_lookup_setup = Clock::now();
    build_metrics_.boundary_lookup_setup_nanoseconds =
        elapsed_nanoseconds(after_scan, after_lookup_setup);

    for (SimplexId simplex_id : vertex_simplex_ids) {
      const auto& vertices = record(simplex_id).vertices;
      if (use_dense_vertex_lookup) {
        SimplexId& stored_id = dense_vertex_to_id[vertices[0]];
        if (stored_id != kInvalidSimplex) {
          throw std::invalid_argument(
              "Simplex tree filtration range contains a duplicate vertex.");
        }
        stored_id = simplex_id;
      } else {
        auto [_, inserted] = vertex_to_id.emplace(vertices[0], simplex_id);
        if (!inserted) {
          throw std::invalid_argument(
              "Simplex tree filtration range contains a duplicate vertex.");
        }
      }
    }
    for (SimplexId simplex_id : edge_simplex_ids) {
      const auto& vertices = record(simplex_id).vertices;
      if (use_dense_edge_lookup) {
        VertexId lhs = vertices[0];
        VertexId rhs = vertices[1];
        if (rhs < lhs) {
          std::swap(lhs, rhs);
        }
        const std::size_t edge_index =
            static_cast<std::size_t>(lhs) * dense_vertex_size + static_cast<std::size_t>(rhs);
        SimplexId& stored_id = dense_edge_to_id[edge_index];
        if (stored_id != kInvalidSimplex) {
          throw std::invalid_argument(
              "Simplex tree filtration range contains a duplicate edge.");
        }
        stored_id = simplex_id;
      } else {
        if (!edge_to_id.emplace(packed_edge_key(vertices[0], vertices[1]), simplex_id)) {
          throw std::invalid_argument(
              "Simplex tree filtration range contains a duplicate edge.");
        }
      }
    }
    const auto after_register = Clock::now();
    build_metrics_.boundary_register_nanoseconds =
        elapsed_nanoseconds(after_lookup_setup, after_register);

    auto require_vertex = [&](VertexId vertex) {
      ++build_metrics_.boundary_vertex_face_lookups;
      if (use_dense_vertex_lookup) {
        if (static_cast<std::size_t>(vertex) >= dense_vertex_to_id.size() ||
            dense_vertex_to_id[vertex] == kInvalidSimplex) {
          throw std::invalid_argument("Simplex tree view is not closed under vertex faces.");
        }
        return dense_vertex_to_id[vertex];
      }
      const auto face_it = vertex_to_id.find(vertex);
      if (face_it == vertex_to_id.end()) {
        throw std::invalid_argument("Simplex tree view is not closed under vertex faces.");
      }
      return face_it->second;
    };
    auto require_edge = [&](VertexId lhs, VertexId rhs) {
      ++build_metrics_.boundary_edge_face_lookups;
      if (use_dense_edge_lookup) {
        if (rhs < lhs) {
          std::swap(lhs, rhs);
        }
        if (static_cast<std::size_t>(lhs) >= dense_vertex_size ||
            static_cast<std::size_t>(rhs) >= dense_vertex_size) {
          throw std::invalid_argument("Simplex tree view is not closed under edge faces.");
        }
        const std::size_t edge_index =
            static_cast<std::size_t>(lhs) * dense_vertex_size + static_cast<std::size_t>(rhs);
        if (dense_edge_to_id[edge_index] == kInvalidSimplex) {
          throw std::invalid_argument("Simplex tree view is not closed under edge faces.");
        }
        return dense_edge_to_id[edge_index];
      }
      const SimplexId face_id = edge_to_id.find(packed_edge_key(lhs, rhs));
      if (face_id == kInvalidSimplex) {
        throw std::invalid_argument("Simplex tree view is not closed under edge faces.");
      }
      return face_id;
    };

    std::vector<VertexId> face_vertices;
    for (SimplexId simplex_id = 0; simplex_id < simplices_.size(); ++simplex_id) {
      auto& simplex = record(simplex_id);
      simplex.boundary.clear();
      const auto& vertices = simplex.vertices;
      if (vertices.size() == 1) {
        continue;
      }

      simplex.boundary.reserve(vertices.size());
      auto append_face = [&](SimplexId face_id) {
        if (record(face_id).filtration > simplex.filtration + 1e-12) {
          throw std::invalid_argument("Simplex tree view filtration is not monotone on faces.");
        }
        simplex.boundary.push_back(face_id);
      };

      if (vertices.size() == 2) {
        append_face(require_vertex(vertices[1]));
        append_face(require_vertex(vertices[0]));
        continue;
      }

      if (vertices.size() == 3) {
        append_face(require_edge(vertices[1], vertices[2]));
        append_face(require_edge(vertices[0], vertices[2]));
        append_face(require_edge(vertices[0], vertices[1]));
        continue;
      }

      face_vertices.resize(vertices.size() - 1);
      for (std::size_t removed = 0; removed < vertices.size(); ++removed) {
        std::size_t target = 0;
        for (std::size_t index = 0; index < vertices.size(); ++index) {
          if (index != removed) {
            face_vertices[target++] = vertices[index];
          }
        }

        ensure_simplex_lookup();
        auto face_it = simplex_to_id_.find(face_vertices);
        if (face_it == simplex_to_id_.end()) {
          throw std::invalid_argument("Simplex tree view is not closed under faces.");
        }
        ++build_metrics_.boundary_general_face_lookups;

        const SimplexId face_id = face_it->second;
        append_face(face_id);
      }
    }
    const auto after_faces = Clock::now();
    build_metrics_.boundary_face_lookup_nanoseconds =
        elapsed_nanoseconds(after_register, after_faces);
  }

  void build_coboundaries() {
    const auto count_start = Clock::now();
    for (auto& simplex : simplices_) {
      simplex.coboundary.clear();
    }

    std::vector<std::size_t> coboundary_sizes(simplices_.size(), 0);
    for (const auto& simplex : simplices_) {
      for (SimplexId face : simplex.boundary) {
        ++coboundary_sizes[face];
      }
    }
    const auto after_count = Clock::now();
    build_metrics_.coboundary_count_nanoseconds =
        elapsed_nanoseconds(count_start, after_count);
    for (SimplexId simplex_id = 0; simplex_id < simplices_.size(); ++simplex_id) {
      simplices_[simplex_id].coboundary.reserve(coboundary_sizes[simplex_id]);
    }
    const auto after_reserve = Clock::now();
    build_metrics_.coboundary_reserve_nanoseconds =
        elapsed_nanoseconds(after_count, after_reserve);

    for (SimplexId simplex_id = 0; simplex_id < simplices_.size(); ++simplex_id) {
      for (SimplexId face : simplices_[simplex_id].boundary) {
        record(face).coboundary.push_back(simplex_id);
      }
    }
    const auto after_fill = Clock::now();
    build_metrics_.coboundary_fill_nanoseconds =
        elapsed_nanoseconds(after_reserve, after_fill);
  }

  void build_orders_and_buckets() {
    const auto size_start = Clock::now();
    filtration_order_.clear();
    filtration_order_.reserve(simplices_.size());
    level_buckets_.assign(level_values_.size(), {});
    if (simplices_.empty()) {
      return;
    }

    std::uint16_t max_dimension = 0;
    std::vector<std::size_t> level_sizes(level_values_.size(), 0);
    for (const auto& simplex : simplices_) {
      if (simplex.dimension > max_dimension) {
        max_dimension = simplex.dimension;
      }
      ++level_sizes[simplex.level];
    }
    for (LevelId level = 0; level < level_buckets_.size(); ++level) {
      level_buckets_[level].reserve(level_sizes[level]);
    }
    const auto after_size = Clock::now();
    build_metrics_.order_size_nanoseconds =
        elapsed_nanoseconds(size_start, after_size);

    const std::size_t dimension_count = static_cast<std::size_t>(max_dimension) + 1;
    const std::size_t bucket_count = level_values_.size() * dimension_count;
    std::vector<std::size_t> bucket_sizes(bucket_count, 0);
    auto bucket_index = [dimension_count](LevelId level, std::uint16_t dimension) {
      return static_cast<std::size_t>(level) * dimension_count +
             static_cast<std::size_t>(dimension);
    };

    for (const auto& simplex : simplices_) {
      ++bucket_sizes[bucket_index(simplex.level, simplex.dimension)];
    }

    std::vector<std::size_t> bucket_offsets(bucket_count + 1, 0);
    for (std::size_t index = 0; index < bucket_count; ++index) {
      bucket_offsets[index + 1] = bucket_offsets[index] + bucket_sizes[index];
    }

    std::vector<SimplexId> bucketed_simplices(simplices_.size());
    std::vector<std::size_t> next_offsets = bucket_offsets;
    for (SimplexId simplex_id = 0; simplex_id < simplices_.size(); ++simplex_id) {
      const auto& simplex = record(simplex_id);
      const std::size_t bucket = bucket_index(simplex.level, simplex.dimension);
      bucketed_simplices[next_offsets[bucket]++] = simplex_id;
    }
    const auto after_bucket = Clock::now();
    build_metrics_.order_bucket_nanoseconds =
        elapsed_nanoseconds(after_size, after_bucket);

    auto emit_sorted_bucket = [&](auto begin, auto end, auto&& vertices_less, LevelId level) {
      bool already_sorted = true;
      for (auto it = begin + (begin == end ? 0 : 1); it != end; ++it) {
        if (vertices_less(*it, *(it - 1))) {
          already_sorted = false;
          break;
        }
      }
      if (!already_sorted) {
        std::sort(begin, end, vertices_less);
      }
      filtration_order_.insert(filtration_order_.end(), begin, end);
      level_buckets_[level].insert(level_buckets_[level].end(), begin, end);
    };

    for (LevelId level = 0; level < level_values_.size(); ++level) {
      for (std::size_t dimension = 0; dimension < dimension_count; ++dimension) {
        const std::size_t bucket = static_cast<std::size_t>(level) * dimension_count + dimension;
        auto begin = bucketed_simplices.begin() +
                     static_cast<std::ptrdiff_t>(bucket_offsets[bucket]);
        auto end = bucketed_simplices.begin() +
                   static_cast<std::ptrdiff_t>(bucket_offsets[bucket + 1]);
        switch (dimension) {
          case 0:
            emit_sorted_bucket(
                begin,
                end,
                [this](SimplexId lhs, SimplexId rhs) {
                  return simplices_[lhs].vertices[0] < simplices_[rhs].vertices[0];
                },
                level);
            break;
          case 1:
            emit_sorted_bucket(
                begin,
                end,
                [this](SimplexId lhs, SimplexId rhs) {
                  const auto& lhs_vertices = simplices_[lhs].vertices;
                  const auto& rhs_vertices = simplices_[rhs].vertices;
                  if (lhs_vertices[0] != rhs_vertices[0]) {
                    return lhs_vertices[0] < rhs_vertices[0];
                  }
                  return lhs_vertices[1] < rhs_vertices[1];
                },
                level);
            break;
          case 2:
            emit_sorted_bucket(
                begin,
                end,
                [this](SimplexId lhs, SimplexId rhs) {
                  const auto& lhs_vertices = simplices_[lhs].vertices;
                  const auto& rhs_vertices = simplices_[rhs].vertices;
                  if (lhs_vertices[0] != rhs_vertices[0]) {
                    return lhs_vertices[0] < rhs_vertices[0];
                  }
                  if (lhs_vertices[1] != rhs_vertices[1]) {
                    return lhs_vertices[1] < rhs_vertices[1];
                  }
                  return lhs_vertices[2] < rhs_vertices[2];
                },
                level);
            break;
          case 3:
            emit_sorted_bucket(
                begin,
                end,
                [this](SimplexId lhs, SimplexId rhs) {
                  const auto& lhs_vertices = simplices_[lhs].vertices;
                  const auto& rhs_vertices = simplices_[rhs].vertices;
                  if (lhs_vertices[0] != rhs_vertices[0]) {
                    return lhs_vertices[0] < rhs_vertices[0];
                  }
                  if (lhs_vertices[1] != rhs_vertices[1]) {
                    return lhs_vertices[1] < rhs_vertices[1];
                  }
                  if (lhs_vertices[2] != rhs_vertices[2]) {
                    return lhs_vertices[2] < rhs_vertices[2];
                  }
                  return lhs_vertices[3] < rhs_vertices[3];
                },
                level);
            break;
          default:
            emit_sorted_bucket(
                begin,
                end,
                [this, dimension](SimplexId lhs, SimplexId rhs) {
                  return vertices_less_in_same_dimension(lhs, rhs, dimension);
                },
                level);
            break;
        }
      }
    }
    const auto after_sort_emit = Clock::now();
    build_metrics_.order_sort_emit_nanoseconds =
        elapsed_nanoseconds(after_bucket, after_sort_emit);
  }

  const SimplexTree* simplex_tree_ = nullptr;
  mutable std::unordered_map<std::vector<VertexId>, SimplexId, VectorHash> simplex_to_id_;
  mutable bool simplex_lookup_built_ = false;
  std::vector<Record> simplices_;
  std::vector<double> level_values_;
  std::vector<std::vector<SimplexId>> level_buckets_;
  std::vector<SimplexId> filtration_order_;
  BuildMetrics build_metrics_;
};

template <class SimplexTree>
FilteredSimplicialComplex filtered_complex_from_simplex_tree(const SimplexTree& simplex_tree) {
  FilteredSimplicialComplex complex;
  for (auto simplex_handle : simplex_tree.filtration_simplex_range()) {
    std::vector<VertexId> vertices;
    for (auto vertex : simplex_tree.simplex_vertex_range(simplex_handle)) {
      vertices.push_back(static_cast<VertexId>(vertex));
    }
    complex.add_simplex(std::move(vertices), simplex_tree.filtration(simplex_handle));
  }
  complex.finalize();
  return complex;
}

}  // namespace internal
}  // namespace morse_persistence
}  // namespace Gudhi
