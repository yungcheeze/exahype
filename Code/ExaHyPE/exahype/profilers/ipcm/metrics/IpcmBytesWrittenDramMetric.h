#ifndef _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_BYTES_WRITTEN_DRAM_METRIC_H_
#define _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_BYTES_WRITTEN_DRAM_METRIC_H_

#ifdef IPCM_AVAILABLE

#include "GenericIpcmMetric.h"

namespace exahype {
namespace profilers {
namespace ipcm {

struct __IpcmBytesWrittenDramMetric {
  using method_return_t = uint64;

  static method_return_t method(const SystemCounterState& before,
                                const SystemCounterState& after) {
    return getBytesWrittenToMC(before, after);
  }

  static const char* method_tag() { return "DramWritten_bytes"; }
};

using IpcmBytesWrittenDramMetric =
    GenericIpcmMetric<__IpcmBytesWrittenDramMetric>;

}  // namespace ipcm
}  // namespace profilers
}  // namespace exahype

#endif  // IPCM_AVAILABLE

#endif  // _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_BYTES_WRITTEN_DRAM_METRIC_H_
