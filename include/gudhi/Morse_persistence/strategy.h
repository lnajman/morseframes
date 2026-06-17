#pragma once

/** \file gudhi/Morse_persistence/strategy.h
 *  \brief Strategy names for same-level Morse sequence construction.
 */

#include <cstdint>
#include <stdexcept>
#include <string>

#include "gudhi/Morse_persistence/internal/morse_reference_api.h"

namespace Gudhi {
namespace morse_persistence {

/** \brief Strategy used to choose same-level Morse pairs and critical simplices.
 *
 *  \ingroup morse_persistence
 */
enum class Morse_sequence_strategy : std::uint8_t {
  SAME_LEVEL_REDUCTION,       /**< Collapse same-level free-face pairs until no pair remains. */
  F_MAX,                      /**< After each saturation phase, choose a maximal fillable simplex. */
  F_MIN,                      /**< After each saturation phase, choose a minimal fillable simplex. */
  PLATEAU_GREEDY,             /**< Choose a fillable simplex that unlocks many almost-pairable cofaces. */

  /** Lower-case alias for SAME_LEVEL_REDUCTION. */
  same_level_reduction = SAME_LEVEL_REDUCTION,
  /** Lower-case alias for F_MAX. */
  f_max = F_MAX,
  /** Lower-case alias for F_MIN. */
  f_min = F_MIN,
  /** Lower-case alias for PLATEAU_GREEDY. */
  plateau_greedy = PLATEAU_GREEDY,
};

/** \brief Convert a GUDHI-style strategy to the internal kernel strategy. */
inline morseframes::MorseSequenceStrategy to_kernel_strategy(Morse_sequence_strategy strategy) {
  switch (strategy) {
    case Morse_sequence_strategy::SAME_LEVEL_REDUCTION:
      return morseframes::MorseSequenceStrategy::SameLevelReduction;
    case Morse_sequence_strategy::F_MAX:
      return morseframes::MorseSequenceStrategy::FMax;
    case Morse_sequence_strategy::F_MIN:
      return morseframes::MorseSequenceStrategy::FMin;
    case Morse_sequence_strategy::PLATEAU_GREEDY:
      return morseframes::MorseSequenceStrategy::PlateauGreedy;
  }
  throw std::invalid_argument("Unknown Morse sequence strategy.");
}

/** \brief Convert an internal kernel strategy to the GUDHI-style strategy. */
inline Morse_sequence_strategy from_kernel_strategy(morseframes::MorseSequenceStrategy strategy) {
  switch (strategy) {
    case morseframes::MorseSequenceStrategy::SameLevelReduction:
      return Morse_sequence_strategy::SAME_LEVEL_REDUCTION;
    case morseframes::MorseSequenceStrategy::FMax:
      return Morse_sequence_strategy::F_MAX;
    case morseframes::MorseSequenceStrategy::FMin:
      return Morse_sequence_strategy::F_MIN;
    case morseframes::MorseSequenceStrategy::PlateauGreedy:
      return Morse_sequence_strategy::PLATEAU_GREEDY;
    case morseframes::MorseSequenceStrategy::Saturated:
    case morseframes::MorseSequenceStrategy::FloodingMax:
    case morseframes::MorseSequenceStrategy::FloodingMin:
    case morseframes::MorseSequenceStrategy::FloodingMinMax:
    case morseframes::MorseSequenceStrategy::FloodingMaxMin:
      break;
  }
  throw std::invalid_argument(
      "Internal Morse sequence strategy is not part of the public GUDHI API.");
}

/** \brief Return the stable command-line/string spelling of a strategy. */
inline const char* strategy_name(Morse_sequence_strategy strategy) {
  return morseframes::morse_sequence_strategy_name(to_kernel_strategy(strategy));
}

/** \brief Parse a strategy name.
 *
 *  The accepted names include
 *  `"same-level-reduction"`, `"f-max"`, `"f-min"`, and `"plateau-greedy"`.
 */
inline Morse_sequence_strategy strategy_from_name(const std::string& strategy) {
  if (strategy == "same-level-reduction") {
    return Morse_sequence_strategy::SAME_LEVEL_REDUCTION;
  }
  if (strategy == "f-max") {
    return Morse_sequence_strategy::F_MAX;
  }
  if (strategy == "f-min") {
    return Morse_sequence_strategy::F_MIN;
  }
  if (strategy == "plateau-greedy") {
    return Morse_sequence_strategy::PLATEAU_GREEDY;
  }
  throw std::invalid_argument("Unknown public Morse sequence strategy: " + strategy);
}

}  // namespace morse_persistence
}  // namespace Gudhi
