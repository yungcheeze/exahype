#ifndef _EXAHYPE_PROFILERS_PROFILER_LIKWID_LIKWID_PROFILER_H_
#define _EXAHYPE_PROFILERS_PROFILER_LIKWID_LIKWID_PROFILER_H_

#ifdef LIKWID_AVAILABLE

#include <iostream>
#include <likwid.h>
#include <memory>
#include <string>
#include <vector>

#include "../Profiler.h"
#include "modules/LikwidModule.h"

namespace exahype {
namespace profilers {
namespace likwid {

struct LikwidProfilerState {
  CpuInfo* cpi_info_ = nullptr;
  CpuTopology* cpu_topology_ = nullptr;
  NumaTopology* numa_topology_ = nullptr;
  AffinityDomains* affinity_domains_ = nullptr;

  std::vector<int> cpus_;
};

class LikwidProfiler : public Profiler {
 public:
  explicit LikwidProfiler(const std::vector<int>& cpus = {});
  virtual ~LikwidProfiler();

  const LikwidProfilerState& state() const { return state_; }
  void addModule(std::unique_ptr<LikwidModule> module);

  void setNumberOfTags(int n) override;
  void registerTag(const std::string& tag) override;
  void start(const std::string& tag) override;
  void stop(const std::string& tag) override;
  void writeToOstream(std::ostream* os) const override;

 private:
  LikwidProfilerState state_;

  std::vector<std::unique_ptr<LikwidModule>> modules_;
};

}  // namespace likwid
}  // namespace profilers
}  // namespace exahype

#endif  // LIKWID_AVAILABLE

#endif  // _EXAHYPE_PROFILERS_PROFILER_LIKWID_LIKWID_PROFILER_H_
