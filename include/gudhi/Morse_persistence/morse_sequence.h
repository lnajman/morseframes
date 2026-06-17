#pragma once

/** \file gudhi/Morse_persistence/morse_sequence.h
 *  \brief Morse sequence aliases used by the public GUDHI-shaped API.
 */

#include "gudhi/Morse_persistence/internal/morse_sequence.h"

namespace Gudhi {
namespace morse_persistence {

/** \brief Type of one Morse sequence step.
 *
 *  \ingroup morse_persistence
 */
using Morse_step_type = internal::MorseStepType;

/** \brief One step of a same-level Morse sequence.
 *
 *  \ingroup morse_persistence
 *
 *  A step is either a critical simplex or a same-level regular face/coface
 *  pair `(sigma, tau)`.
 */
using Morse_step = internal::MorseStep;

/** \brief Ordered Morse sequence used by the reference-map recurrence.
 *
 *  \ingroup morse_persistence
 */
using Morse_sequence = internal::MorseSequence;

}  // namespace morse_persistence
}  // namespace Gudhi
