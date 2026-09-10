#pragma once

#include <map>
#include <string>
#include "morseframes/morse_sequence.hpp"

// Optional diagnostic children; compile against historical header snapshots too.
inline std::map<std::string, double> pls_phase_profile(
    const morseframes::MorseSequenceBuildMetrics& metrics) {
  std::map<std::string, double> result;
#ifdef MORSEFRAMES_PLS_PHASE_PROFILE_VERSION
#define PLS_PHASE(name) result[#name] = 1e-9 * metrics.process_lower_stars_##name##_nanoseconds
  PLS_PHASE(output_init); PLS_PHASE(vertex_order); PLS_PHASE(executor_init);
  PLS_PHASE(storage_init); PLS_PHASE(owner_keys); PLS_PHASE(partition);
  PLS_PHASE(schedule); PLS_PHASE(execution); PLS_PHASE(cleanup);
  PLS_PHASE(events_index_cleanup); PLS_PHASE(keys_cleanup);
  PLS_PHASE(membership_cleanup); PLS_PHASE(executor_cleanup);
  PLS_PHASE(vertices_cleanup);
#undef PLS_PHASE
#endif
  return result;
}
