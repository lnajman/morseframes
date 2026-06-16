#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <utility>
#include <vector>

#include "morseframes/field_arithmetic.hpp"

namespace morseframes {

using CriticalId = std::uint32_t;

template <std::size_t InlineCapacity>
class SmallAnnotation {
 public:
  using value_type = CriticalId;
  using iterator = CriticalId*;
  using const_iterator = const CriticalId*;

  SmallAnnotation() = default;

  SmallAnnotation(std::initializer_list<CriticalId> labels) {
    reserve(labels.size());
    for (CriticalId label : labels) {
      push_back(label);
    }
  }

  ~SmallAnnotation() { delete[] overflow_; }

  SmallAnnotation(const SmallAnnotation& other) {
    copy_from(other);
  }

  SmallAnnotation& operator=(const SmallAnnotation& other) {
    if (this != &other) {
      reserve(other.size_);
      std::copy(other.begin(), other.end(), data());
      size_ = other.size_;
    }
    return *this;
  }

  SmallAnnotation(SmallAnnotation&& other) noexcept {
    move_from(std::move(other));
  }

  SmallAnnotation& operator=(SmallAnnotation&& other) noexcept {
    if (this != &other) {
      delete[] overflow_;
      overflow_ = nullptr;
      size_ = 0;
      capacity_ = InlineCapacity;
      move_from(std::move(other));
    }
    return *this;
  }

  void clear() { size_ = 0; }

  void reserve(std::size_t capacity) {
    if (capacity > capacity_) {
      grow_to(capacity);
    }
  }

  void push_back(CriticalId label) {
    if (size_ == capacity_) {
      grow_to(std::max<std::size_t>(capacity_ * 2, InlineCapacity + 1));
    }
    data()[size_++] = label;
  }

  iterator insert(const_iterator position, CriticalId label) {
    const std::size_t index = static_cast<std::size_t>(position - begin());
    if (size_ == capacity_) {
      grow_to(std::max<std::size_t>(capacity_ * 2, InlineCapacity + 1));
    }
    CriticalId* values = data();
    std::move_backward(values + index, values + size_, values + size_ + 1);
    values[index] = label;
    ++size_;
    return data() + index;
  }

  iterator erase(const_iterator position) {
    const std::size_t index = static_cast<std::size_t>(position - begin());
    CriticalId* values = data();
    std::move(values + index + 1, values + size_, values + index);
    --size_;
    return data() + index;
  }

  bool empty() const { return size_ == 0; }
  std::size_t size() const { return size_; }

  const CriticalId& operator[](std::size_t index) const { return data()[index]; }
  CriticalId& operator[](std::size_t index) { return data()[index]; }

  const CriticalId& front() const { return (*this)[0]; }
  CriticalId& front() { return (*this)[0]; }

  const CriticalId& back() const { return (*this)[size_ - 1]; }
  CriticalId& back() { return (*this)[size_ - 1]; }

  const_iterator begin() const { return data(); }
  const_iterator end() const { return data() + size_; }
  iterator begin() { return data(); }
  iterator end() { return data() + size_; }

  void swap(SmallAnnotation& other) noexcept {
    using std::swap;
    swap(inline_, other.inline_);
    swap(overflow_, other.overflow_);
    swap(size_, other.size_);
    swap(capacity_, other.capacity_);
  }

  operator std::vector<CriticalId>() const {
    return std::vector<CriticalId>(begin(), end());
  }

 private:
  void grow_to(std::size_t new_capacity) {
    CriticalId* replacement = new CriticalId[new_capacity];
    std::copy(begin(), end(), replacement);
    delete[] overflow_;
    overflow_ = replacement;
    capacity_ = new_capacity;
  }

  const CriticalId* data() const {
    return overflow_ == nullptr ? inline_.data() : overflow_;
  }

  CriticalId* data() {
    return overflow_ == nullptr ? inline_.data() : overflow_;
  }

  void copy_from(const SmallAnnotation& other) {
    size_ = other.size_;
    capacity_ = other.size_ > InlineCapacity ? other.size_ : InlineCapacity;
    if (capacity_ > InlineCapacity) {
      overflow_ = new CriticalId[capacity_];
    }
    std::copy(other.begin(), other.end(), data());
  }

  void move_from(SmallAnnotation&& other) {
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

  std::array<CriticalId, InlineCapacity> inline_{};
  CriticalId* overflow_ = nullptr;
  std::size_t size_ = 0;
  std::size_t capacity_ = InlineCapacity;
};

template <std::size_t InlineCapacity>
bool operator==(const SmallAnnotation<InlineCapacity>& lhs,
                const SmallAnnotation<InlineCapacity>& rhs) {
  return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template <std::size_t InlineCapacity>
bool operator!=(const SmallAnnotation<InlineCapacity>& lhs,
                const SmallAnnotation<InlineCapacity>& rhs) {
  return !(lhs == rhs);
}

using Annotation = SmallAnnotation<2>;

struct FieldAnnotationEntry {
  CriticalId label = 0;
  std::uint32_t coefficient = 0;
};

using FieldAnnotation = std::vector<FieldAnnotationEntry>;

inline bool contains_label(const Annotation& annotation, CriticalId label) {
  return std::binary_search(annotation.begin(), annotation.end(), label);
}

inline void toggle_label(Annotation& annotation, CriticalId label) {
  auto it = std::lower_bound(annotation.begin(), annotation.end(), label);
  if (it != annotation.end() && *it == label) {
    annotation.erase(it);
  } else {
    annotation.insert(it, label);
  }
}

inline void xor_annotations_in_place(Annotation& lhs, const Annotation& rhs) {
  if (rhs.empty()) {
    return;
  }
  if (lhs.empty()) {
    lhs = rhs;
    return;
  }

  Annotation result;
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

inline void xor_annotations_in_place(Annotation& lhs,
                                     const Annotation& rhs,
                                     Annotation& scratch) {
  if (rhs.empty()) {
    return;
  }
  if (lhs.empty()) {
    lhs = rhs;
    return;
  }

  scratch.clear();
  scratch.reserve(lhs.size() + rhs.size());

  auto left = lhs.begin();
  auto right = rhs.begin();
  while (left != lhs.end() || right != rhs.end()) {
    if (right == rhs.end() || (left != lhs.end() && *left < *right)) {
      scratch.push_back(*left);
      ++left;
    } else if (left == lhs.end() || *right < *left) {
      scratch.push_back(*right);
      ++right;
    } else {
      ++left;
      ++right;
    }
  }

  lhs.swap(scratch);
  scratch.clear();
}

inline void add_scaled_field_annotation_in_place(FieldAnnotation& target,
                                                 const FieldAnnotation& source,
                                                 std::uint32_t scale,
                                                 std::uint32_t modulus) {
  scale %= modulus;
  if (scale == 0 || source.empty()) {
    return;
  }

  FieldAnnotation result;
  result.reserve(target.size() + source.size());
  std::size_t left = 0;
  std::size_t right = 0;

  while (left < target.size() || right < source.size()) {
    if (right == source.size() ||
        (left < target.size() && target[left].label < source[right].label)) {
      result.push_back(target[left]);
      ++left;
      continue;
    }

    if (left == target.size() || source[right].label < target[left].label) {
      const auto coefficient = modp_multiply(scale, source[right].coefficient, modulus);
      if (coefficient != 0) {
        result.push_back(FieldAnnotationEntry{source[right].label, coefficient});
      }
      ++right;
      continue;
    }

    const auto scaled = modp_multiply(scale, source[right].coefficient, modulus);
    const auto coefficient = static_cast<std::uint32_t>(
        (static_cast<std::uint64_t>(target[left].coefficient) + scaled) % modulus);
    if (coefficient != 0) {
      result.push_back(FieldAnnotationEntry{target[left].label, coefficient});
    }
    ++left;
    ++right;
  }

  target = std::move(result);
}

inline void scale_field_annotation_in_place(FieldAnnotation& annotation,
                                            std::uint32_t scale,
                                            std::uint32_t modulus) {
  scale %= modulus;
  if (scale == 0) {
    annotation.clear();
    return;
  }

  FieldAnnotation result;
  result.reserve(annotation.size());
  for (const auto& entry : annotation) {
    const auto coefficient = modp_multiply(entry.coefficient, scale, modulus);
    if (coefficient != 0) {
      result.push_back(FieldAnnotationEntry{entry.label, coefficient});
    }
  }
  annotation = std::move(result);
}

inline std::uint32_t field_annotation_coefficient(const FieldAnnotation& annotation,
                                                  CriticalId label) {
  auto it = std::lower_bound(annotation.begin(),
                             annotation.end(),
                             label,
                             [](const FieldAnnotationEntry& entry, CriticalId value) {
                               return entry.label < value;
                             });
  if (it == annotation.end() || it->label != label) {
    return 0;
  }
  return it->coefficient;
}

inline void remove_label_from_field_annotation(FieldAnnotation& annotation,
                                               CriticalId label) {
  auto it = std::lower_bound(annotation.begin(),
                             annotation.end(),
                             label,
                             [](const FieldAnnotationEntry& entry, CriticalId value) {
                               return entry.label < value;
                             });
  if (it != annotation.end() && it->label == label) {
    annotation.erase(it);
  }
}

}  // namespace morseframes
