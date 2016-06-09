#ifndef _EXAHYPE_PROFILERS_SIMPLE_NO_OP_PROFILER_H_
#define _EXAHYPE_PROFILERS_SIMPLE_NO_OP_PROFILER_H_

#include <iostream>
#include <string>

#include "../Profiler.h"

namespace exahype {
namespace profilers {
namespace simple {

class NoOpProfiler : public Profiler {
 public:
  NoOpProfiler() {}

  virtual ~NoOpProfiler() { writeToCout(); }

  void setNumberOfTags(int n) override;
  void registerTag(const std::string& tag) override;
  void start(const std::string& tag) override;
  void stop(const std::string& tag) override;
  virtual void writeToOstream(std::ostream* os) const;
};

}  // namespace simple
}  // namespace profilers
}  // namespace exahype

#endif  // _EXAHYPE_PROFILERS_SIMPLE_NO_OP_PROFILER_H_
