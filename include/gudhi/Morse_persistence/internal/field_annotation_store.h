#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "gudhi/Morse_persistence/internal/annotation.h"
#include "gudhi/Morse_persistence/internal/filtered_complex.h"
#include "gudhi/Morse_persistence/internal/inverse_annotation_store.h"

namespace Gudhi { namespace morse_persistence { namespace internal {

class FieldAnnotationStore {
 public:
  FieldAnnotationStore(const std::vector<FieldAnnotation>& full_annotations,
                       const std::vector<SimplexId>& selected_simplices,
                       std::size_t num_labels,
                       std::uint32_t modulus)
      : simplex_to_local_(full_annotations.size(), kInvalidLocal),
        inverse_lists_(num_labels),
        num_labels_(num_labels),
        modulus_(modulus) {
    validate_prime_field_characteristic(modulus_);
    annotations_.reserve(selected_simplices.size());
    annotation_positions_.reserve(selected_simplices.size());

    for (SimplexId simplex : selected_simplices) {
      if (simplex >= full_annotations.size()) {
        throw std::out_of_range("Working set contains an invalid simplex id.");
      }
      if (simplex_to_local_[simplex] != kInvalidLocal) {
        throw std::invalid_argument("Working set contains a duplicate simplex id.");
      }
      const std::size_t local = annotations_.size();
      simplex_to_local_[simplex] = local;
      annotations_.push_back(full_annotations[simplex]);
      initialize_inverse_entries(local, annotations_.back());
    }
  }

  FieldAnnotationStore(std::vector<FieldAnnotation> annotations,
                       const std::vector<SimplexId>& selected_simplices,
                       std::size_t universe_size,
                       std::size_t num_labels,
                       std::uint32_t modulus)
      : annotations_(std::move(annotations)),
        simplex_to_local_(universe_size, kInvalidLocal),
        inverse_lists_(num_labels),
        num_labels_(num_labels),
        modulus_(modulus) {
    validate_prime_field_characteristic(modulus_);
    if (annotations_.size() != selected_simplices.size()) {
      throw std::invalid_argument(
          "Compacted field annotation table size does not match the working set.");
    }

    annotation_positions_.reserve(selected_simplices.size());
    for (std::size_t local = 0; local < selected_simplices.size(); ++local) {
      const SimplexId simplex = selected_simplices[local];
      if (simplex >= universe_size) {
        throw std::out_of_range("Working set contains an invalid simplex id.");
      }
      if (simplex_to_local_[simplex] != kInvalidLocal) {
        throw std::invalid_argument("Working set contains a duplicate simplex id.");
      }
      simplex_to_local_[simplex] = local;
      initialize_inverse_entries(local, annotations_[local]);
    }
  }

  const FieldAnnotation& annotation(SimplexId simplex) const {
    return annotations_.at(local_index(simplex));
  }

  const FieldAnnotation& annotation_by_local_unchecked(std::size_t local) const {
    return annotations_[local];
  }

  std::size_t size() const { return annotations_.size(); }

  const InverseAnnotationStoreMetrics& metrics() const { return metrics_; }

  void remove_label_from_all(CriticalId label) {
    check_label(label);
    auto& candidates = inverse_lists_[label];
    while (!candidates.empty()) {
      const std::size_t local = candidates.back();
      ++metrics_.remove_candidate_scans;
      remove_label_from_annotation(local, label);
      ++metrics_.remove_applied;
    }
  }

  void eliminate_pivot(CriticalId pivot,
                       const FieldAnnotation& update,
                       std::uint32_t pivot_coefficient) {
    check_label(pivot);
    validate_field_annotation(update);
    if (pivot_coefficient == 0) {
      throw std::invalid_argument("Cannot eliminate a zero pivot coefficient.");
    }
    const std::uint32_t update_pivot =
        field_annotation_coefficient(update, pivot);
    if (update_pivot != pivot_coefficient) {
      throw std::invalid_argument("Pivot coefficient does not match the update annotation.");
    }

    const std::uint32_t inverse_pivot = modp_inverse(pivot_coefficient, modulus_);
    auto& candidates = inverse_lists_[pivot];
    while (!candidates.empty()) {
      const std::size_t local = candidates.back();
      const std::uint32_t coefficient =
          field_annotation_coefficient(annotations_.at(local), pivot);
      if (coefficient == 0) {
        throw std::logic_error("Exact inverse list pointed to a missing field label.");
      }
      ++metrics_.xor_candidate_scans;
      ++metrics_.xor_applied;
      const std::uint32_t scale =
          (modulus_ - modp_multiply(coefficient, inverse_pivot, modulus_)) % modulus_;
      add_scaled_into_annotation(local, update, scale);
    }
  }

 private:
  static constexpr std::size_t kInvalidLocal = std::numeric_limits<std::size_t>::max();

  void check_label(CriticalId label) const {
    if (label >= inverse_lists_.size()) {
      throw std::out_of_range("Invalid critical label.");
    }
  }

  std::size_t local_index(SimplexId simplex) const {
    if (simplex >= simplex_to_local_.size() || simplex_to_local_[simplex] == kInvalidLocal) {
      throw std::out_of_range("Simplex annotation is outside the stored working set.");
    }
    return simplex_to_local_[simplex];
  }

  void validate_field_annotation(const FieldAnnotation& annotation) const {
    CriticalId previous_label = 0;
    bool has_previous = false;
    for (const auto& entry : annotation) {
      check_label(entry.label);
      if (entry.coefficient == 0 || entry.coefficient >= modulus_) {
        throw std::invalid_argument("Field annotation has an invalid coefficient.");
      }
      if (has_previous && entry.label <= previous_label) {
        throw std::invalid_argument("Field annotation labels must be strictly sorted.");
      }
      previous_label = entry.label;
      has_previous = true;
    }
  }

  void initialize_inverse_entries(std::size_t local, const FieldAnnotation& annotation) {
    validate_field_annotation(annotation);

    if (!annotation.empty()) {
      ++metrics_.initial_nonempty_annotations;
      metrics_.initial_total_annotation_size += annotation.size();
      if (annotation.size() > metrics_.initial_max_annotation_size) {
        metrics_.initial_max_annotation_size = annotation.size();
      }
    }

    annotation_positions_.push_back({});
    auto& positions = annotation_positions_.back();
    positions.reserve(annotation.size());
    for (const auto& entry : annotation) {
      const CriticalId label = entry.label;
      positions.push_back(inverse_lists_[label].size());
      inverse_lists_[label].push_back(local);
      ++metrics_.initial_inverse_list_entries;
    }
  }

  void remove_label_from_annotation(std::size_t local, CriticalId label) {
    auto& annotation = annotations_.at(local);
    auto it = std::lower_bound(
        annotation.begin(),
        annotation.end(),
        label,
        [](const FieldAnnotationEntry& entry, CriticalId value) {
          return entry.label < value;
        });
    if (it == annotation.end() || it->label != label) {
      throw std::logic_error("Exact inverse list pointed to a missing field label.");
    }
    metrics_.remove_total_annotation_size += annotation.size();
    if (annotation.size() > metrics_.remove_max_annotation_size) {
      metrics_.remove_max_annotation_size = annotation.size();
    }
    const std::size_t label_index = static_cast<std::size_t>(it - annotation.begin());
    remove_inverse_entry(label, local, annotation_positions_.at(local).at(label_index));
    annotation.erase(it);
    annotation_positions_[local].erase(annotation_positions_[local].begin() + label_index);
  }

  void remove_inverse_entry(CriticalId label, std::size_t local, std::size_t position) {
    auto& inverse = inverse_lists_.at(label);
    if (position >= inverse.size() || inverse[position] != local) {
      throw std::logic_error("Field inverse annotation position is inconsistent.");
    }

    const std::size_t moved_local = inverse.back();
    inverse[position] = moved_local;
    inverse.pop_back();

    if (position < inverse.size()) {
      auto& moved_annotation = annotations_.at(moved_local);
      auto it = std::lower_bound(
          moved_annotation.begin(),
          moved_annotation.end(),
          label,
          [](const FieldAnnotationEntry& entry, CriticalId value) {
            return entry.label < value;
          });
      if (it == moved_annotation.end() || it->label != label) {
        throw std::logic_error("Moved field inverse annotation entry is inconsistent.");
      }
      const std::size_t moved_label_index =
          static_cast<std::size_t>(it - moved_annotation.begin());
      annotation_positions_.at(moved_local).at(moved_label_index) = position;
    }
  }

  void add_scaled_into_annotation(std::size_t local,
                                  const FieldAnnotation& update,
                                  std::uint32_t scale) {
    scale %= modulus_;
    if (scale == 0 || update.empty()) {
      return;
    }

    const FieldAnnotation* changed_labels = &update;
    FieldAnnotation copied_update;
    if (&annotations_.at(local) == &update) {
      copied_update = update;
      changed_labels = &copied_update;
    }

    metrics_.xor_changed_labels += changed_labels->size();

    FieldAnnotation old_annotation = std::move(annotations_.at(local));
    std::vector<std::size_t> old_positions = std::move(annotation_positions_.at(local));
    const std::size_t input_size = old_annotation.size() + changed_labels->size();
    metrics_.xor_total_input_size += input_size;
    if (input_size > metrics_.xor_max_input_size) {
      metrics_.xor_max_input_size = input_size;
    }

    FieldAnnotation new_annotation;
    std::vector<std::size_t> new_positions;
    new_annotation.reserve(old_annotation.size() + changed_labels->size());
    new_positions.reserve(old_annotation.size() + changed_labels->size());

    std::size_t left = 0;
    std::size_t right = 0;
    while (left < old_annotation.size() || right < changed_labels->size()) {
      if (right == changed_labels->size() ||
          (left < old_annotation.size() &&
           old_annotation[left].label < (*changed_labels)[right].label)) {
        new_annotation.push_back(old_annotation[left]);
        new_positions.push_back(old_positions[left]);
        ++left;
        continue;
      }

      if (left == old_annotation.size() ||
          (*changed_labels)[right].label < old_annotation[left].label) {
        const CriticalId label = (*changed_labels)[right].label;
        check_label(label);
        const std::uint32_t coefficient =
            modp_multiply(scale, (*changed_labels)[right].coefficient, modulus_);
        if (coefficient != 0) {
          new_annotation.push_back(FieldAnnotationEntry{label, coefficient});
          new_positions.push_back(inverse_lists_[label].size());
          inverse_lists_[label].push_back(local);
          ++metrics_.inverse_list_appends;
          ++metrics_.xor_inserted_labels;
        }
        ++right;
        continue;
      }

      const CriticalId label = old_annotation[left].label;
      const std::uint32_t scaled =
          modp_multiply(scale, (*changed_labels)[right].coefficient, modulus_);
      const std::uint32_t coefficient = static_cast<std::uint32_t>(
          (static_cast<std::uint64_t>(old_annotation[left].coefficient) + scaled) %
          modulus_);
      if (coefficient != 0) {
        new_annotation.push_back(FieldAnnotationEntry{label, coefficient});
        new_positions.push_back(old_positions[left]);
      } else {
        remove_inverse_entry(label, local, old_positions[left]);
        ++metrics_.xor_removed_labels;
      }
      ++left;
      ++right;
    }

    metrics_.xor_total_output_size += new_annotation.size();
    if (new_annotation.size() > metrics_.xor_max_output_size) {
      metrics_.xor_max_output_size = new_annotation.size();
    }

    annotations_[local] = std::move(new_annotation);
    annotation_positions_[local] = std::move(new_positions);
  }

  std::vector<FieldAnnotation> annotations_;
  std::vector<std::vector<std::size_t>> annotation_positions_;
  std::vector<std::size_t> simplex_to_local_;
  std::vector<std::vector<std::size_t>> inverse_lists_;
  std::size_t num_labels_ = 0;
  std::uint32_t modulus_ = 2;
  InverseAnnotationStoreMetrics metrics_;
};

}  // namespace internal
}  // namespace morse_persistence
}  // namespace Gudhi
