#pragma once

#include <cstdint>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

#include "morseframes/annotation.hpp"
#include "morseframes/critical_order.hpp"
#include "morseframes/filtered_complex.hpp"
#include "morseframes/inverse_annotation_store.hpp"
#include "morseframes/morse_sequence.hpp"
#include "morseframes/reference_persistence.hpp"
#include "morseframes/working_sets.hpp"

namespace morseframes {

class MorseCoreferenceComputer {
 public:
  MorseCoreferenceComputer(const FilteredSimplicialComplex& complex, const MorseSequence& sequence)
      : complex_(complex), sequence_(sequence) {}

  std::vector<Annotation> compute_full_coreferences() const {
    std::vector<Annotation> coreferences(complex_.size());

    for (auto step_it = sequence_.steps().rbegin(); step_it != sequence_.steps().rend(); ++step_it) {
      const MorseStep& step = *step_it;
      if (step.type == MorseStepType::Critical) {
        const auto critical_id = checked_critical_id(sequence_.critical_index(step.sigma));
        coreferences[step.sigma] = Annotation{critical_id};
        continue;
      }

      coreferences[step.sigma].clear();
      Annotation upper_coreference;
      for (SimplexId coface : complex_.coboundary(step.sigma)) {
        if (coface != step.tau) {
          xor_annotations_in_place(upper_coreference, coreferences[coface]);
        }
      }
      coreferences[step.tau] = std::move(upper_coreference);
    }

    return coreferences;
  }

 private:
  static CriticalId checked_critical_id(std::int32_t value) {
    if (value < 0) {
      throw std::logic_error("Expected a critical simplex.");
    }
    return static_cast<CriticalId>(value);
  }

  const FilteredSimplicialComplex& complex_;
  const MorseSequence& sequence_;
};

struct MorseCoreferenceFrame {
  MorseSequence sequence;
  std::vector<Annotation> coreferences;
};

struct MorseFieldCoreferenceFrame {
  MorseSequence sequence;
  std::vector<FieldAnnotation> coreferences;
  std::uint32_t modulus = 2;
};

inline std::uint32_t boundary_coefficient_of_face(
    const FilteredSimplicialComplex& complex,
    SimplexId coface,
    SimplexId face,
    std::uint32_t modulus) {
  const auto& boundary = complex.boundary(coface);
  for (std::size_t removed_index = 0; removed_index < boundary.size(); ++removed_index) {
    if (boundary[removed_index] == face) {
      return boundary_coefficient(removed_index, modulus);
    }
  }
  throw std::logic_error("Expected a codimension-one face/coface incidence.");
}

class MorseFieldCoreferenceComputer {
 public:
  MorseFieldCoreferenceComputer(const FilteredSimplicialComplex& complex,
                                const MorseSequence& sequence,
                                std::uint32_t modulus)
      : complex_(complex), sequence_(sequence), modulus_(modulus) {
    validate_prime_field_characteristic(modulus_);
  }

  std::vector<FieldAnnotation> compute_full_coreferences() const {
    std::vector<FieldAnnotation> coreferences(complex_.size());

    for (auto step_it = sequence_.steps().rbegin();
         step_it != sequence_.steps().rend();
         ++step_it) {
      update_coreference_for_step(*step_it, coreferences);
    }

    return coreferences;
  }

 private:
  static CriticalId checked_critical_id(std::int32_t value) {
    if (value < 0) {
      throw std::logic_error("Expected a critical simplex.");
    }
    return static_cast<CriticalId>(value);
  }

  void update_coreference_for_step(const MorseStep& step,
                                   std::vector<FieldAnnotation>& coreferences) const {
    if (step.type == MorseStepType::Critical) {
      coreferences[step.sigma] = FieldAnnotation{
          FieldAnnotationEntry{checked_critical_id(sequence_.critical_index(step.sigma)), 1}};
      return;
    }

    coreferences[step.sigma].clear();
    FieldAnnotation upper_coreference;
    std::uint32_t paired_coface_coefficient = 0;
    bool found_paired_coface = false;
    for (SimplexId coface : complex_.coboundary(step.sigma)) {
      const std::uint32_t coefficient =
          boundary_coefficient_of_face(complex_, coface, step.sigma, modulus_);
      if (coface == step.tau) {
        paired_coface_coefficient = coefficient;
        found_paired_coface = true;
        continue;
      }
      add_scaled_field_annotation_in_place(
          upper_coreference, coreferences[coface], coefficient, modulus_);
    }

    if (!found_paired_coface) {
      throw std::logic_error("Regular pair is not a face/coface pair.");
    }

    const std::uint32_t scale =
        (modulus_ - modp_inverse(paired_coface_coefficient, modulus_)) % modulus_;
    scale_field_annotation_in_place(upper_coreference, scale, modulus_);
    coreferences[step.tau] = std::move(upper_coreference);
  }

  const FilteredSimplicialComplex& complex_;
  const MorseSequence& sequence_;
  std::uint32_t modulus_ = 2;
};

class MorseCoreferenceFrameBuilder {
 public:
  explicit MorseCoreferenceFrameBuilder(const FilteredSimplicialComplex& complex)
      : complex_(complex) {}

  MorseCoreferenceFrame build_coreduction() const {
    return build_same_level_reduction();
  }

  MorseCoreferenceFrame build_same_level_reduction() const {
    const std::size_t n = complex_.size();

    struct CollapsePair {
      SimplexId sigma = kInvalidSimplex;
      SimplexId tau = kInvalidSimplex;
    };

    struct LevelRecord {
      LevelId level = 0;
      std::vector<SimplexId> criticals;
      std::vector<CollapsePair> collapse_pairs;
    };

    std::vector<LevelRecord> records;
    records.reserve(complex_.num_levels());
    std::vector<bool> active(n, false);
    std::vector<std::uint32_t> active_coface_count(n, 0);

    for (LevelId level = 0; level < complex_.num_levels(); ++level) {
      const auto& bucket = complex_.simplices_of_level(level);
      LevelRecord record;
      record.level = level;
      record.criticals.reserve(bucket.size());
      record.collapse_pairs.reserve(bucket.size() / 2);
      std::queue<SimplexId> free_faces;

      for (SimplexId simplex : bucket) {
        active[simplex] = true;
        active_coface_count[simplex] = 0;
      }

      auto count_active_same_level_cofaces = [&](SimplexId simplex) {
        std::uint32_t count = 0;
        for (SimplexId coface : complex_.coboundary(simplex)) {
          if (complex_.level(coface) == level && active[coface]) {
            ++count;
          }
        }
        return count;
      };

      for (SimplexId simplex : bucket) {
        active_coface_count[simplex] = count_active_same_level_cofaces(simplex);
        if (active_coface_count[simplex] == 1) {
          free_faces.push(simplex);
        }
      }

      auto unique_active_same_level_coface = [&](SimplexId simplex) {
        SimplexId result = kInvalidSimplex;
        for (SimplexId coface : complex_.coboundary(simplex)) {
          if (complex_.level(coface) != level || !active[coface]) {
            continue;
          }
          if (result != kInvalidSimplex) {
            return kInvalidSimplex;
          }
          result = coface;
        }
        return result;
      };

      auto decrement_active_faces_of = [&](SimplexId simplex) {
        for (SimplexId face : complex_.boundary(simplex)) {
          if (complex_.level(face) != level || !active[face]) {
            continue;
          }
          if (active_coface_count[face] == 0) {
            throw std::logic_error("Active coface count underflow.");
          }
          --active_coface_count[face];
          if (active_coface_count[face] == 1) {
            free_faces.push(face);
          }
        }
      };

      while (!free_faces.empty()) {
        const SimplexId sigma = free_faces.front();
        free_faces.pop();
        if (!active[sigma] || active_coface_count[sigma] != 1) {
          continue;
        }

        const SimplexId tau = unique_active_same_level_coface(sigma);
        if (tau == kInvalidSimplex) {
          continue;
        }

        record.collapse_pairs.push_back(CollapsePair{sigma, tau});
        active[sigma] = false;
        active[tau] = false;
        decrement_active_faces_of(tau);
        decrement_active_faces_of(sigma);
      }

      for (SimplexId simplex : bucket) {
        if (!active[simplex]) {
          continue;
        }
        active[simplex] = false;
        record.criticals.push_back(simplex);
      }

      records.push_back(std::move(record));
    }

    MorseSequence sequence(n);
    std::vector<bool> inserted(n, false);

    auto emit_critical = [&](SimplexId sigma, LevelId level) {
      if (inserted[sigma]) {
        throw std::logic_error("Tried to insert a simplex twice.");
      }
      for (SimplexId face : complex_.boundary(sigma)) {
        if (!inserted[face]) {
          throw std::logic_error(
              "Same-level reduction critical has a missing boundary face.");
        }
      }
      sequence.add_critical(sigma, level);
      inserted[sigma] = true;
    };

    auto emit_pair = [&](SimplexId sigma, SimplexId tau, LevelId level) {
      if (inserted[sigma] || inserted[tau]) {
        throw std::logic_error("Tried to insert a regular pair twice.");
      }
      for (SimplexId face : complex_.boundary(tau)) {
        if (face != sigma && !inserted[face]) {
          throw std::logic_error("Same-level reduction pair has a missing boundary face.");
        }
      }
      sequence.add_regular_pair(sigma, tau, level);
      inserted[sigma] = true;
      inserted[tau] = true;
    };

    for (const auto& record : records) {
      for (SimplexId critical : record.criticals) {
        emit_critical(critical, record.level);
      }
      for (std::size_t index = record.collapse_pairs.size(); index > 0; --index) {
        const auto& pair = record.collapse_pairs[index - 1];
        emit_pair(pair.sigma, pair.tau, record.level);
      }
    }

    std::vector<Annotation> coreferences(n);

    for (std::size_t level_index = records.size(); level_index > 0; --level_index) {
      const auto& record = records[level_index - 1];

      for (const auto& pair : record.collapse_pairs) {
        coreferences[pair.sigma].clear();
        Annotation upper_coreference;
        for (SimplexId coface : complex_.coboundary(pair.sigma)) {
          if (coface != pair.tau) {
            xor_annotations_in_place(upper_coreference, coreferences[coface]);
          }
        }
        coreferences[pair.tau] = std::move(upper_coreference);
      }

      for (std::size_t index = record.criticals.size(); index > 0; --index) {
        const SimplexId critical = record.criticals[index - 1];
        coreferences[critical] = Annotation{
            checked_critical_id(sequence.critical_index(critical))};
      }
    }

    return MorseCoreferenceFrame{std::move(sequence), std::move(coreferences)};
  }

 private:
  static CriticalId checked_critical_id(std::int32_t value) {
    if (value < 0) {
      throw std::logic_error("Expected a critical simplex.");
    }
    return static_cast<CriticalId>(value);
  }

  const FilteredSimplicialComplex& complex_;
};

class MorseCoreferencePersistenceReducer {
 public:
  MorseCoreferencePersistenceReducer(const FilteredSimplicialComplex& complex,
                                     const MorseSequence& sequence,
                                     const std::vector<Annotation>& coreferences)
      : complex_(complex),
        sequence_(sequence),
        annotations_(coreferences,
                     coreference_working_set(complex, sequence),
                     sequence.critical_simplices().size()) {}

  PersistenceDiagram compute() {
    return compute_with_metrics().diagram;
  }

  MorseReferenceReductionResult compute_with_metrics() {
    const auto& critical_simplices = sequence_.critical_simplices();
    std::vector<bool> active_dual(critical_simplices.size(), false);
    std::vector<bool> killed_dual(critical_simplices.size(), false);
    MorseReferenceReductionResult result;
    result.metrics.working_set_size = annotations_.size();
    result.metrics.critical_count = critical_simplices.size();

    for (auto order_it = critical_order_.critical_ids.rbegin();
         order_it != critical_order_.critical_ids.rend();
         ++order_it) {
      const CriticalId sigma_critical_id = *order_it;
      const SimplexId sigma = critical_simplices[sigma_critical_id];

      Annotation coboundary_annotation;
      for (SimplexId coface : complex_.coboundary(sigma)) {
        const auto& coface_annotation = annotations_.annotation(coface);
        ++result.metrics.boundary_annotation_xors;
        result.metrics.boundary_annotation_total_input_size += coface_annotation.size();
        xor_annotations_in_place(coboundary_annotation, coface_annotation);
      }
      if (coboundary_annotation.size() > result.metrics.boundary_annotation_max_size) {
        result.metrics.boundary_annotation_max_size = coboundary_annotation.size();
      }
      result.metrics.boundary_annotation_total_output_size += coboundary_annotation.size();
      if (coboundary_annotation.size() >
          result.metrics.boundary_annotation_max_output_size) {
        result.metrics.boundary_annotation_max_output_size = coboundary_annotation.size();
      }

      if (coboundary_annotation.empty()) {
        active_dual[sigma_critical_id] = true;
        continue;
      }

      const CriticalId pivot = earliest_critical_label(coboundary_annotation, critical_order_);
      const SimplexId death = critical_simplices.at(pivot);

      result.diagram.finite_pairs.push_back(PersistencePair{
          sigma,
          death,
          complex_.dimension(sigma),
          complex_.filtration(sigma),
          complex_.filtration(death),
      });

      killed_dual[pivot] = true;
      ++result.metrics.pivot_eliminations;

      remove_dual_negative_label(sigma_critical_id);
      eliminate_dual_pivot(pivot, coboundary_annotation);
    }

    for (CriticalId id = 0; id < critical_simplices.size(); ++id) {
      if (!active_dual[id] || killed_dual[id]) {
        continue;
      }
      const SimplexId birth = critical_simplices[id];
      result.diagram.essential.push_back(EssentialInterval{
          birth,
          complex_.dimension(birth),
          complex_.filtration(birth),
      });
    }

    result.metrics.finite_pairs = result.diagram.finite_pairs.size();
    result.metrics.essential_intervals = result.diagram.essential.size();
    result.metrics.inverse_store = annotations_.metrics();
    return result;
  }

 private:
  static CriticalId checked_critical_id(std::int32_t value) {
    if (value < 0) {
      throw std::logic_error("Expected a critical simplex.");
    }
    return static_cast<CriticalId>(value);
  }

  void remove_dual_negative_label(CriticalId label) {
    annotations_.remove_label_from_all(label);
  }

  void eliminate_dual_pivot(CriticalId pivot, const Annotation& coboundary_annotation) {
    annotations_.xor_into_all_containing(pivot, coboundary_annotation);
  }

  const FilteredSimplicialComplex& complex_;
  const MorseSequence& sequence_;
  InverseAnnotationStore annotations_;
  CriticalOrder critical_order_ = build_flooding_critical_order(complex_, sequence_);
};

class MorseFieldCoreferencePersistenceReducer {
 public:
  MorseFieldCoreferencePersistenceReducer(const FilteredSimplicialComplex& complex,
                                          const MorseSequence& sequence,
                                          std::vector<FieldAnnotation> coreferences,
                                          std::uint32_t modulus)
      : complex_(complex),
        sequence_(sequence),
        coreferences_(std::move(coreferences)),
        modulus_(modulus) {
    validate_prime_field_characteristic(modulus_);
    if (coreferences_.size() != complex_.size()) {
      throw std::invalid_argument("Coreference table size does not match the complex size.");
    }
  }

  PersistenceDiagram compute() {
    const auto& critical_simplices = sequence_.critical_simplices();
    std::vector<bool> active_dual(critical_simplices.size(), false);
    std::vector<bool> killed_dual(critical_simplices.size(), false);
    PersistenceDiagram diagram;

    for (auto order_it = critical_order_.critical_ids.rbegin();
         order_it != critical_order_.critical_ids.rend();
         ++order_it) {
      const CriticalId sigma_critical_id = *order_it;
      const SimplexId sigma = critical_simplices[sigma_critical_id];

      FieldAnnotation coboundary_annotation;
      for (SimplexId coface : complex_.coboundary(sigma)) {
        add_scaled_field_annotation_in_place(
            coboundary_annotation,
            coreferences_[coface],
            boundary_coefficient_of_face(complex_, coface, sigma, modulus_),
            modulus_);
      }

      if (coboundary_annotation.empty()) {
        active_dual[sigma_critical_id] = true;
        continue;
      }

      const FieldAnnotationEntry pivot_entry =
          earliest_critical_entry(coboundary_annotation, critical_order_);
      const CriticalId pivot = pivot_entry.label;
      const std::uint32_t pivot_coefficient = pivot_entry.coefficient;
      const SimplexId death = critical_simplices.at(pivot);

      diagram.finite_pairs.push_back(PersistencePair{
          sigma,
          death,
          complex_.dimension(sigma),
          complex_.filtration(sigma),
          complex_.filtration(death),
      });

      killed_dual[pivot] = true;

      for (auto& annotation : coreferences_) {
        remove_label_from_field_annotation(annotation, sigma_critical_id);
      }

      const std::uint32_t inverse_pivot = modp_inverse(pivot_coefficient, modulus_);
      for (auto& annotation : coreferences_) {
        const std::uint32_t coefficient =
            field_annotation_coefficient(annotation, pivot);
        if (coefficient == 0) {
          continue;
        }
        const std::uint32_t scale =
            (modulus_ - modp_multiply(coefficient, inverse_pivot, modulus_)) % modulus_;
        add_scaled_field_annotation_in_place(
            annotation, coboundary_annotation, scale, modulus_);
      }
    }

    for (CriticalId id = 0; id < critical_simplices.size(); ++id) {
      if (!active_dual[id] || killed_dual[id]) {
        continue;
      }
      const SimplexId birth = critical_simplices[id];
      diagram.essential.push_back(EssentialInterval{
          birth,
          complex_.dimension(birth),
          complex_.filtration(birth),
      });
    }

    return diagram;
  }

 private:
  static CriticalId checked_critical_id(std::int32_t value) {
    if (value < 0) {
      throw std::logic_error("Expected a critical simplex.");
    }
    return static_cast<CriticalId>(value);
  }

  const FilteredSimplicialComplex& complex_;
  const MorseSequence& sequence_;
  std::vector<FieldAnnotation> coreferences_;
  std::uint32_t modulus_ = 2;
  CriticalOrder critical_order_ = build_flooding_critical_order(complex_, sequence_);
};

class MorseCompactFieldCoreferencePersistenceReducer {
 public:
  MorseCompactFieldCoreferencePersistenceReducer(
      const FilteredSimplicialComplex& complex,
      const MorseSequence& sequence,
      const std::vector<FieldAnnotation>& coreferences,
      std::uint32_t modulus)
      : complex_(complex),
        sequence_(sequence),
        annotations_(coreferences,
                     coreference_working_set(complex, sequence),
                     sequence.critical_simplices().size(),
                     modulus),
        modulus_(modulus) {
    validate_prime_field_characteristic(modulus_);
  }

  MorseCompactFieldCoreferencePersistenceReducer(
      const FilteredSimplicialComplex& complex,
      const MorseSequence& sequence,
      std::vector<SimplexId> working_set,
      std::vector<FieldAnnotation> annotations,
      std::uint32_t modulus)
      : complex_(complex),
        sequence_(sequence),
        working_set_(std::move(working_set)),
        annotations_(std::move(annotations),
                     working_set_,
                     complex.size(),
                     sequence.critical_simplices().size(),
                     modulus),
        modulus_(modulus) {
    validate_prime_field_characteristic(modulus_);
  }

  PersistenceDiagram compute() {
    const auto& critical_simplices = sequence_.critical_simplices();
    std::vector<bool> active_dual(critical_simplices.size(), false);
    std::vector<bool> killed_dual(critical_simplices.size(), false);
    PersistenceDiagram diagram;

    for (auto order_it = critical_order_.critical_ids.rbegin();
         order_it != critical_order_.critical_ids.rend();
         ++order_it) {
      const CriticalId sigma_critical_id = *order_it;
      const SimplexId sigma = critical_simplices[sigma_critical_id];

      FieldAnnotation coboundary_annotation;
      for (SimplexId coface : complex_.coboundary(sigma)) {
        add_scaled_field_annotation_in_place(
            coboundary_annotation,
            annotations_.annotation(coface),
            boundary_coefficient_of_face(complex_, coface, sigma, modulus_),
            modulus_);
      }

      if (coboundary_annotation.empty()) {
        active_dual[sigma_critical_id] = true;
        continue;
      }

      const FieldAnnotationEntry pivot_entry =
          earliest_critical_entry(coboundary_annotation, critical_order_);
      const CriticalId pivot = pivot_entry.label;
      const std::uint32_t pivot_coefficient = pivot_entry.coefficient;
      const SimplexId death = critical_simplices.at(pivot);

      diagram.finite_pairs.push_back(PersistencePair{
          sigma,
          death,
          complex_.dimension(sigma),
          complex_.filtration(sigma),
          complex_.filtration(death),
      });

      killed_dual[pivot] = true;

      annotations_.remove_label_from_all(sigma_critical_id);
      annotations_.eliminate_pivot(pivot, coboundary_annotation, pivot_coefficient);
    }

    for (CriticalId id = 0; id < critical_simplices.size(); ++id) {
      if (!active_dual[id] || killed_dual[id]) {
        continue;
      }
      const SimplexId birth = critical_simplices[id];
      diagram.essential.push_back(EssentialInterval{
          birth,
          complex_.dimension(birth),
          complex_.filtration(birth),
      });
    }

    return diagram;
  }

 private:
  static CriticalId checked_critical_id(std::int32_t value) {
    if (value < 0) {
      throw std::logic_error("Expected a critical simplex.");
    }
    return static_cast<CriticalId>(value);
  }

  const FilteredSimplicialComplex& complex_;
  const MorseSequence& sequence_;
  std::vector<SimplexId> working_set_;
  FieldAnnotationStore annotations_;
  std::uint32_t modulus_ = 2;
  CriticalOrder critical_order_ = build_flooding_critical_order(complex_, sequence_);
};

inline PersistenceDiagram compute_morse_coreference_persistence(
    const FilteredSimplicialComplex& complex, const MorseSequence& sequence) {
  MorseCoreferenceComputer coreference_computer(complex, sequence);
  auto coreferences = coreference_computer.compute_full_coreferences();
  MorseCoreferencePersistenceReducer reducer(complex, sequence, coreferences);
  return reducer.compute();
}

inline PersistenceDiagram compute_morse_coreference_prime_field_persistence(
    const FilteredSimplicialComplex& complex,
    const MorseSequence& sequence,
    std::uint32_t modulus) {
  auto coreferences =
      MorseFieldCoreferenceComputer(complex, sequence, modulus).compute_full_coreferences();
  auto working_set = coreference_working_set(complex, sequence);
  std::vector<FieldAnnotation> annotations;
  annotations.reserve(working_set.size());
  for (SimplexId simplex : working_set) {
    annotations.push_back(std::move(coreferences.at(simplex)));
  }
  MorseCompactFieldCoreferencePersistenceReducer reducer(
      complex, sequence, std::move(working_set), std::move(annotations), modulus);
  return reducer.compute();
}

}  // namespace morseframes
