#pragma once

/** \file gudhi/Morse_persistence.h
 *  \brief Umbrella header for the GUDHI Morse persistence API.
 *
 *  This header exposes a thin public layer over the current Morse persistence
 *  kernels. The main entry point is
 *  `Gudhi::morse_persistence::compute_morse_persistence`.
 *
 *  \ingroup morse_persistence
 *
 *  Include this header to use the public Morse persistence API.
 */

#include "gudhi/Morse_persistence/complex_view.h"
#include "gudhi/Morse_persistence/diagram.h"
#include "gudhi/Morse_persistence/morse_sequence.h"
#include "gudhi/Morse_persistence/persistence_reducer.h"
#include "gudhi/Morse_persistence/reference_map.h"
#include "gudhi/Morse_persistence/strategy.h"
