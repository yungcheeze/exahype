#ifndef _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_CYCLES_LOST_L2_MISSES_METRIC_H_
#define _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_CYCLES_LOST_L2_MISSES_METRIC_H_

#ifdef IPCM_AVAILABLE

#include "GenericIpcmMetric.h"

namespace exahype {
namespace profilers {
namespace ipcm {

struct __IpcmCyclesLostL2MissesMetric {
  using method_return_t = double;

  static method_return_t method(const SystemCounterState& before,
                                const SystemCounterState& after) {
    return getCyclesLostDueL2CacheMisses(before, after);
  }

  static const char* method_tag() { return "CyclesLostL2Misses"; }
};

using IpcmCyclesLostL2MissesMetric =
    GenericIpcmMetric<__IpcmCyclesLostL2MissesMetric>;

}  // namespace ipcm
}  // namespace profilers
}  // namespace exahype

#endif  // IPCM_AVAILABLE

#endif  // _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_CYCLES_LOST_L2_MISSES_METRIC_H_
