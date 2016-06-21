#ifndef _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_L2_CACHE_MISSES_METRIC_H_
#define _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_L2_CACHE_MISSES_METRIC_H_

#ifdef IPCM_AVAILABLE

#include "GenericIpcmMetric.h"

namespace exahype {
namespace profilers {
namespace ipcm {

struct __IpcmL2CacheMissesMetric {
  using method_return_t = uint64;

  static method_return_t method(const SystemCounterState& before,
                                const SystemCounterState& after) {
    return getL2CacheMisses(before, after);
  }

  static const char* method_tag() { return "L2CacheMisses"; }
};

using IpcmL2CacheMissesMetric = GenericIpcmMetric<__IpcmL2CacheMissesMetric>;

}  // namespace ipcm
}  // namespace profilers
}  // namespace exahype

#endif  // IPCM_AVAILABLE

#endif  // _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_L2_CACHE_MISSES_METRIC_H_
