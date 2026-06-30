#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "morseframes/annotation.hpp"
#include "morseframes/complex_view.hpp"
#include "morseframes/morse_sequence.hpp"

namespace morseframes {

struct CriticalOrder {
  std::vector<CriticalId> critical_ids;
  std::vector<std::size_t> rank_by_id;
};

template <class ComplexView>
CriticalOrder build_flooding_critical_order(const ComplexView& complex,
                                            const MorseSequence& sequence) {
  static_assert(is_complex_view_v<ComplexView>,
                "build_flooding_critical_order requires a Morse complex-view type.");
  const auto& critical_simplices = sequence.critical_simplices();
  CriticalOrder order;
  order.critical_ids.reserve(critical_simplices.size());
  order.rank_by_id.resize(critical_simplices.size());

  for (CriticalId id = 0; id < critical_simplices.size(); ++id) {
    order.critical_ids.push_back(id);
  }

  std::stable_sort(
      order.critical_ids.begin(),
      order.critical_ids.end(),
      [&](CriticalId lhs, CriticalId rhs) {
        return complex.level(critical_simplices[lhs]) <
               complex.level(critical_simplices[rhs]);
      });

  for (std::size_t rank = 0; rank < order.critical_ids.size(); ++rank) {
    order.rank_by_id[order.critical_ids[rank]] = rank;
  }
  return order;
}

inline std::size_t critical_order_rank(const CriticalOrder& order, CriticalId label) {
  if (label >= order.rank_by_id.size()) {
    throw std::out_of_range("Critical label is outside the critical order.");
  }
  return order.rank_by_id[label];
}

inline CriticalId latest_critical_label(const Annotation& annotation,
                                        const CriticalOrder& order) {
  if (annotation.empty()) {
    throw std::invalid_argument("Cannot choose a pivot from an empty annotation.");
  }
  CriticalId best = annotation.front();
  std::size_t best_rank = critical_order_rank(order, best);
  for (CriticalId label : annotation) {
    const std::size_t rank = critical_order_rank(order, label);
    if (rank > best_rank) {
      best = label;
      best_rank = rank;
    }
  }
  return best;
}

inline CriticalId earliest_critical_label(const Annotation& annotation,
                                          const CriticalOrder& order) {
  if (annotation.empty()) {
    throw std::invalid_argument("Cannot choose a pivot from an empty annotation.");
  }
  CriticalId best = annotation.front();
  std::size_t best_rank = critical_order_rank(order, best);
  for (CriticalId label : annotation) {
    const std::size_t rank = critical_order_rank(order, label);
    if (rank < best_rank) {
      best = label;
      best_rank = rank;
    }
  }
  return best;
}

inline FieldAnnotationEntry latest_critical_entry(const FieldAnnotation& annotation,
                                                  const CriticalOrder& order) {
  if (annotation.empty()) {
    throw std::invalid_argument("Cannot choose a pivot from an empty annotation.");
  }
  FieldAnnotationEntry best = annotation.front();
  std::size_t best_rank = critical_order_rank(order, best.label);
  for (const FieldAnnotationEntry& entry : annotation) {
    const std::size_t rank = critical_order_rank(order, entry.label);
    if (rank > best_rank) {
      best = entry;
      best_rank = rank;
    }
  }
  return best;
}

inline FieldAnnotationEntry earliest_critical_entry(const FieldAnnotation& annotation,
                                                    const CriticalOrder& order) {
  if (annotation.empty()) {
    throw std::invalid_argument("Cannot choose a pivot from an empty annotation.");
  }
  FieldAnnotationEntry best = annotation.front();
  std::size_t best_rank = critical_order_rank(order, best.label);
  for (const FieldAnnotationEntry& entry : annotation) {
    const std::size_t rank = critical_order_rank(order, entry.label);
    if (rank < best_rank) {
      best = entry;
      best_rank = rank;
    }
  }
  return best;
}

}  // namespace morseframes
