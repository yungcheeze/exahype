#ifndef _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_CYCLES_METRIC_H_
#define _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_CYCLES_METRIC_H_

#ifdef IPCM_AVAILABLE

#include "GenericIpcmMetric.h"

// TODO(guera): Guard

namespace exahype {
namespace profilers {
namespace ipcm {

struct __IpcmCyclesMetric {
  using method_return_t = uint64;

  static method_return_t method(const SystemCounterState& before,
                                const SystemCounterState& after) {
    return getCycles(before, after);
  }

  static const char* method_tag() { return "Cycles"; }
};

using IpcmCyclesMetric = GenericIpcmMetric<__IpcmCyclesMetric>;

}  // namespace ipcm
}  // namespace profilers
}  // namespace exahype

#endif  // IPCM_AVAILABLE

#endif  // _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_CYCLES_METRIC_H_
