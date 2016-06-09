#ifndef _EXAHYPE_PROFILERS_SIMPLE_CHRONO_ELAPSED_TIME_PROFILER_H_
#define _EXAHYPE_PROFILERS_SIMPLE_CHRONO_ELAPSED_TIME_PROFILER_H_

#include <chrono>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>

#include "../Profiler.h"

namespace exahype {
namespace profilers {
namespace simple {

class ChronoElapsedTimeProfiler : public Profiler {
 public:
  ChronoElapsedTimeProfiler() {}

  virtual ~ChronoElapsedTimeProfiler() { writeToCout(); }

  void setNumberOfTags(int n) override;
  void registerTag(const std::string& tag) override;
  void start(const std::string& tag) override;
  void stop(const std::string& tag) override;
  void writeToOstream(std::ostream* os) const override;

 private:
  std::unordered_map<std::string,
                     std::chrono::time_point<std::chrono::system_clock>>
      time_points_;
  std::unordered_map<
      std::string, std::pair<int, std::chrono::high_resolution_clock::duration>>
      counts_and_durations_;
};

}  // namespace simple
}  // namespace profilers
}  // namespace exahype

#endif  // _EXAHYPE_PROFILERS_SIMPLE_CHRONO_ELAPSED_TIME_PROFILER_H_
