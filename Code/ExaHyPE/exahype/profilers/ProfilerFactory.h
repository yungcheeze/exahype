#ifndef _EXAHYPE_PROFILERS_PROFILER_FACTORY_H_
#define _EXAHYPE_PROFILERS_PROFILER_FACTORY_H_

#include <memory>
#include <string>
#include <vector>

namespace exahype {
namespace profilers {
class Profiler;
}  // namespace profilers
}  // namespace exahype

namespace exahype {
namespace profilers {

class ProfilerFactory {
 public:
  virtual ~ProfilerFactory() {}

  static ProfilerFactory& getInstance();

  std::unique_ptr<Profiler> create(const std::string& profiler_name,
                                   const std::vector<std::string>& modules);

 private:
  ProfilerFactory() {}
};

}  // namespace profilers
}  // namespace exahype

#endif  // _EXAHYPE_PROFILERS_PROFILER_FACTORY_H_
