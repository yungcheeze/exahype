#ifndef _EXAHYPE_PROFILERS_LIKWID_MODULES_LIKWID_TIME_MEASUREMENT_MODULE_H_
#define _EXAHYPE_PROFILERS_LIKWID_MODULES_LIKWID_TIME_MEASUREMENT_MODULE_H_

#ifdef LIKWID_AVAILABLE

#include <iostream>
#include <likwid.h>
#include <string>
#include <unordered_map>

#include "LikwidModule.h"

namespace exahype {
namespace profilers {
namespace likwid {

class LikwidTimeMeasurementModule : public LikwidModule {
 public:
  explicit LikwidTimeMeasurementModule(const LikwidProfilerState& state);
  virtual ~LikwidTimeMeasurementModule();

  void setNumberOfTags(int n) override;
  void registerTag(const std::string& tag) override;
  void start(const std::string& tag) override;
  void stop(const std::string& tag) override;
  void writeToOstream(std::ostream* os) const override;

 private:
  std::unordered_map<std::string, TimerData> timer_data_;
};

}  // namespace likwid
}  // namespace profilers
}  // namespace exahype

#endif  // LIKWID_AVAILABLE

#endif  // _EXAHYPE_PROFILERS_LIKWID_MODULES_LIKWID_TIME_MEASUREMENT_MODULE_H_
