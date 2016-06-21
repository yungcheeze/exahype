#ifndef _EXAHYPE_PROFILERS_IPCM_IPCM_PROFILER_H_
#define _EXAHYPE_PROFILERS_IPCM_IPCM_PROFILER_H_

#ifdef IPCM_AVAILABLE

#include <cpucounters.h>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "../Profiler.h"

namespace exahype {
namespace profilers {
namespace ipcm {

class IpcmMetric;

class IpcmProfiler : public Profiler {
 public:
  IpcmProfiler();
  virtual ~IpcmProfiler();

  void addMetric(std::unique_ptr<IpcmMetric> metric);

  void setNumberOfTags(int n) override;
  void registerTag(const std::string& tag) override;
  void start(const std::string& tag) override;
  void stop(const std::string& tag) override;
  void writeToOstream(std::ostream* os) const override;

 private:
  PCM* pcm_ = nullptr;
  std::unordered_map<std::string, SystemCounterState> system_counter_states_;
  std::vector<std::unique_ptr<IpcmMetric>> metrics_;
};

}  // namespace ipcm
}  // namespace profilers
}  // namespace exahype

#endif  // IPCM_AVAILABLE

#endif  // _EXAHYPE_PROFILERS_IPCM_IPCM_PROFILER_H_
