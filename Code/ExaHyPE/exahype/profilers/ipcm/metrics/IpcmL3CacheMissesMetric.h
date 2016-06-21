#ifndef _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_L3_CACHE_MISSES_METRIC_H_
#define _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_L3_CACHE_MISSES_METRIC_H_

#ifdef IPCM_AVAILABLE

#include "GenericIpcmMetric.h"

namespace exahype {
namespace profilers {
namespace ipcm {

struct __IpcmL3CacheMissesMetric {
  using method_return_t = uint64;

  static method_return_t method(const SystemCounterState& before,
                                const SystemCounterState& after) {
    return getL3CacheMisses(before, after);
  }

  static const char* method_tag() { return "L3CacheMisses"; }
};

using IpcmL3CacheMissesMetric = GenericIpcmMetric<__IpcmL3CacheMissesMetric>;

}  // namespace ipcm
}  // namespace profilers
}  // namespace exahype

#endif  // IPCM_AVAILABLE

#endif  // _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_L3_CACHE_MISSES_METRIC_H_
