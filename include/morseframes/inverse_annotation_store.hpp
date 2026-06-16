#pragma once

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "morseframes/annotation.hpp"
#include "morseframes/filtered_complex.hpp"

namespace morseframes {

struct InverseAnnotationStoreMetrics {
  std::size_t initial_nonempty_annotations = 0;
  std::size_t initial_total_annotation_size = 0;
  std::size_t initial_max_annotation_size = 0;
  std::size_t initial_inverse_list_entries = 0;
  std::size_t remove_candidate_scans = 0;
  std::size_t remove_applied = 0;
  std::size_t remove_total_annotation_size = 0;
  std::size_t remove_max_annotation_size = 0;
  std::size_t xor_candidate_scans = 0;
  std::size_t xor_applied = 0;
  std::size_t xor_changed_labels = 0;
  std::size_t xor_total_input_size = 0;
  std::size_t xor_total_output_size = 0;
  std::size_t xor_max_input_size = 0;
  std::size_t xor_max_output_size = 0;
  std::size_t xor_inserted_labels = 0;
  std::size_t xor_removed_labels = 0;
  std::size_t inverse_list_appends = 0;
};

class InverseAnnotationStore {
 public:
  InverseAnnotationStore(std::vector<Annotation> annotations, std::size_t num_labels)
      : annotations_(std::move(annotations)),
        simplex_to_local_(annotations_.size(), kInvalidLocal),
        inverse_lists_(num_labels) {
    reserve_inverse_lists_from_annotations(annotations_);
    annotation_positions_.reserve(annotations_.size());
    for (std::size_t local = 0; local < annotations_.size(); ++local) {
      simplex_to_local_[local] = checked_local_id(local);
      initialize_inverse_entries(local, annotations_[local]);
    }
  }

  InverseAnnotationStore(const std::vector<Annotation>& full_annotations,
                         const std::vector<SimplexId>& selected_simplices,
                         std::size_t num_labels)
      : simplex_to_local_(full_annotations.size(), kInvalidLocal), inverse_lists_(num_labels) {
    reserve_inverse_lists_from_selected_annotations(full_annotations, selected_simplices);
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
      simplex_to_local_[simplex] = checked_local_id(local);
      annotations_.push_back(full_annotations[simplex]);
      initialize_inverse_entries(local, annotations_.back());
    }
  }

  InverseAnnotationStore(std::vector<Annotation> annotations,
                         const std::vector<SimplexId>& selected_simplices,
                         std::size_t universe_size,
                         std::size_t num_labels)
      : annotations_(std::move(annotations)),
        simplex_to_local_(universe_size, kInvalidLocal),
        inverse_lists_(num_labels) {
    if (annotations_.size() != selected_simplices.size()) {
      throw std::invalid_argument(
          "Compacted annotation table size does not match the working set.");
    }

    reserve_inverse_lists_from_annotations(annotations_);
    annotation_positions_.reserve(annotations_.size());

    for (std::size_t local = 0; local < selected_simplices.size(); ++local) {
      const SimplexId simplex = selected_simplices[local];
      if (simplex >= universe_size) {
        throw std::out_of_range("Working set contains an invalid simplex id.");
      }
      if (simplex_to_local_[simplex] != kInvalidLocal) {
        throw std::invalid_argument("Working set contains a duplicate simplex id.");
      }
      simplex_to_local_[simplex] = checked_local_id(local);
      initialize_inverse_entries(local, annotations_[local]);
    }
  }

  const Annotation& annotation(SimplexId simplex) const {
    return annotations_[local_index(simplex)];
  }

  const Annotation& annotation_unchecked(SimplexId simplex) const {
    return annotations_[static_cast<std::size_t>(simplex_to_local_[simplex])];
  }

  const Annotation& annotation_by_local_unchecked(std::size_t local) const {
    return annotations_[local];
  }

  const std::vector<Annotation>& annotations() const { return annotations_; }
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

  void xor_into_all_containing(CriticalId pivot, const Annotation& update) {
    check_label(pivot);
    if (contains_label(update, pivot)) {
      xor_into_all_containing_known_pivot(pivot, update);
      return;
    }

    const auto candidates = inverse_lists_[pivot];
    metrics_.xor_candidate_scans += candidates.size();
    for (std::size_t local : candidates) {
      if (contains_label(annotations_.at(local), pivot)) {
        ++metrics_.xor_applied;
        xor_into_annotation(local, update);
      }
    }
  }

  void xor_into_all_containing_known_pivot(CriticalId pivot, const Annotation& update) {
    check_label(pivot);
    auto& candidates = inverse_lists_[pivot];
    while (!candidates.empty()) {
      const std::size_t local = candidates.back();
      ++metrics_.xor_candidate_scans;
      ++metrics_.xor_applied;
      xor_into_annotation(local, update);
    }
  }

 private:
  using LocalId = SimplexId;

  static constexpr LocalId kInvalidLocal = std::numeric_limits<LocalId>::max();

  static LocalId checked_local_id(std::size_t local) {
    if (local >= std::numeric_limits<LocalId>::max()) {
      throw std::overflow_error("Too many local annotations.");
    }
    return static_cast<LocalId>(local);
  }

  void check_label(CriticalId label) const {
    if (label >= inverse_lists_.size()) {
      throw std::out_of_range("Invalid critical label.");
    }
  }

  std::size_t local_index(SimplexId simplex) const {
    if (simplex >= simplex_to_local_.size() || simplex_to_local_[simplex] == kInvalidLocal) {
      throw std::out_of_range("Simplex annotation is outside the stored working set.");
    }
    return static_cast<std::size_t>(simplex_to_local_[simplex]);
  }

  void reserve_inverse_lists_from_annotations(const std::vector<Annotation>& annotations) {
    std::vector<std::size_t> counts(inverse_lists_.size(), 0);
    for (const Annotation& annotation : annotations) {
      for (CriticalId label : annotation) {
        check_label(label);
        ++counts[label];
      }
    }
    for (std::size_t label = 0; label < counts.size(); ++label) {
      inverse_lists_[label].reserve(counts[label]);
    }
  }

  void reserve_inverse_lists_from_selected_annotations(
      const std::vector<Annotation>& full_annotations,
      const std::vector<SimplexId>& selected_simplices) {
    std::vector<std::size_t> counts(inverse_lists_.size(), 0);
    for (SimplexId simplex : selected_simplices) {
      if (simplex >= full_annotations.size()) {
        throw std::out_of_range("Working set contains an invalid simplex id.");
      }
      for (CriticalId label : full_annotations[simplex]) {
        check_label(label);
        ++counts[label];
      }
    }
    for (std::size_t label = 0; label < counts.size(); ++label) {
      inverse_lists_[label].reserve(counts[label]);
    }
  }

  void initialize_inverse_entries(std::size_t local, const Annotation& annotation) {
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
    for (CriticalId label : annotation) {
      check_label(label);
      positions.push_back(inverse_lists_[label].size());
      inverse_lists_[label].push_back(checked_local_id(local));
      ++metrics_.initial_inverse_list_entries;
    }
  }

  void remove_label_from_annotation(std::size_t local, CriticalId label) {
    auto& annotation = annotations_[local];
    auto it = std::lower_bound(annotation.begin(), annotation.end(), label);
    if (it == annotation.end() || *it != label) {
      throw std::logic_error("Exact inverse list pointed to a missing label.");
    }
    metrics_.remove_total_annotation_size += annotation.size();
    if (annotation.size() > metrics_.remove_max_annotation_size) {
      metrics_.remove_max_annotation_size = annotation.size();
    }
    const std::size_t label_index = static_cast<std::size_t>(it - annotation.begin());
    remove_inverse_entry(label, local, annotation_positions_[local][label_index]);
    annotation.erase(it);
    annotation_positions_[local].erase(annotation_positions_[local].begin() + label_index);
  }

  void remove_inverse_entry(CriticalId label, std::size_t local, std::size_t position) {
    auto& inverse = inverse_lists_[label];
    const LocalId local_id = checked_local_id(local);
    if (position >= inverse.size() || inverse[position] != local_id) {
      throw std::logic_error("Inverse annotation position is inconsistent.");
    }

    const LocalId moved_local_id = inverse.back();
    inverse[position] = moved_local_id;
    inverse.pop_back();

    if (position < inverse.size()) {
      const std::size_t moved_local = static_cast<std::size_t>(moved_local_id);
      auto& moved_annotation = annotations_[moved_local];
      auto it = std::lower_bound(moved_annotation.begin(), moved_annotation.end(), label);
      if (it == moved_annotation.end() || *it != label) {
        throw std::logic_error("Moved inverse annotation entry is inconsistent.");
      }
      const std::size_t moved_label_index =
          static_cast<std::size_t>(it - moved_annotation.begin());
      annotation_positions_[moved_local][moved_label_index] = position;
    }
  }

  void xor_into_annotation(std::size_t local, const Annotation& update) {
    const Annotation* changed_labels = &update;
    if (&annotations_[local] == &update) {
      scratch_update_ = update;
      changed_labels = &scratch_update_;
    }

    metrics_.xor_changed_labels += changed_labels->size();

    const Annotation& old_annotation = annotations_[local];
    const std::vector<std::size_t>& old_positions = annotation_positions_[local];
    const std::size_t input_size = old_annotation.size() + changed_labels->size();
    metrics_.xor_total_input_size += input_size;
    if (input_size > metrics_.xor_max_input_size) {
      metrics_.xor_max_input_size = input_size;
    }

    scratch_annotation_.clear();
    scratch_positions_.clear();
    scratch_annotation_.reserve(old_annotation.size() + changed_labels->size());
    scratch_positions_.reserve(old_annotation.size() + changed_labels->size());

    std::size_t left = 0;
    std::size_t right = 0;
    while (left < old_annotation.size() || right < changed_labels->size()) {
      if (
          right == changed_labels->size()
          || (left < old_annotation.size() && old_annotation[left] < (*changed_labels)[right])) {
        scratch_annotation_.push_back(old_annotation[left]);
        scratch_positions_.push_back(old_positions[left]);
        ++left;
      } else if (
          left == old_annotation.size()
          || (*changed_labels)[right] < old_annotation[left]) {
        const CriticalId label = (*changed_labels)[right];
        check_label(label);
        scratch_annotation_.push_back(label);
        scratch_positions_.push_back(inverse_lists_[label].size());
        inverse_lists_[label].push_back(checked_local_id(local));
        ++metrics_.inverse_list_appends;
        ++metrics_.xor_inserted_labels;
        ++right;
      } else {
        const CriticalId label = old_annotation[left];
        remove_inverse_entry(label, local, old_positions[left]);
        ++metrics_.xor_removed_labels;
        ++left;
        ++right;
      }
    }

    metrics_.xor_total_output_size += scratch_annotation_.size();
    if (scratch_annotation_.size() > metrics_.xor_max_output_size) {
      metrics_.xor_max_output_size = scratch_annotation_.size();
    }

    annotations_[local].swap(scratch_annotation_);
    annotation_positions_[local].swap(scratch_positions_);
  }

  std::vector<Annotation> annotations_;
  std::vector<std::vector<std::size_t>> annotation_positions_;
  std::vector<LocalId> simplex_to_local_;
  std::vector<std::vector<LocalId>> inverse_lists_;
  Annotation scratch_annotation_;
  Annotation scratch_update_;
  std::vector<std::size_t> scratch_positions_;
  InverseAnnotationStoreMetrics metrics_;
};

}  // namespace morseframes
