#ifndef _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_CYCLES_LOST_L3_MISSES_METRIC_H_
#define _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_CYCLES_LOST_L3_MISSES_METRIC_H_

#ifdef IPCM_AVAILABLE

#include "GenericIpcmMetric.h"

namespace exahype {
namespace profilers {
namespace ipcm {

struct __IpcmCyclesLostL3MissesMetric {
  using method_return_t = double;

  static method_return_t method(const SystemCounterState& before,
                                const SystemCounterState& after) {
    return getCyclesLostDueL3CacheMisses(before, after);
  }

  static const char* method_tag() { return "CyclesLostL3Misses"; }
};

using IpcmCyclesLostL3MissesMetric =
    GenericIpcmMetric<__IpcmCyclesLostL3MissesMetric>;

}  // namespace ipcm
}  // namespace profilers
}  // namespace exahype

#endif  // IPCM_AVAILABLE

#endif  // _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_CYCLES_LOST_L3_MISSES_METRIC_H_
