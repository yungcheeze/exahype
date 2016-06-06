#ifndef _EXAHYPE_PROFILERS_PROFILER_H_
#define _EXAHYPE_PROFILERS_PROFILER_H_

#include <iostream>
#include <string>

namespace exahype {
namespace profilers {

class Profiler {
 public:
  Profiler() {}

  virtual ~Profiler() {
    // TODO(guera): Remove this
    // TODO(guera): I'm not sure if it is a good idea to (indirectly) call
    // virtual functions inside the destructor.
    writeToCout();
  }

  // Disallow copy and assignment
  Profiler(const Profiler& other) = delete;
  Profiler& operator=(const Profiler& other) = delete;

  virtual void setNumberOfTags(int n) = 0;
  virtual void registerTag(const std::string& tag) = 0;
  virtual void start(const std::string& tag) = 0;
  virtual void stop(const std::string& tag) = 0;
  virtual void writeToOstream(std::ostream* os) const = 0;

  void writeToCout() const;
  void writeToFile(const std::string& path) const;
};

}  // namespace profilers
}  // namespace exahype

#endif  // EXAHYPE_PROFILERS_PROFILER_H_
