#include "IpcmCountMetric.h"

#ifdef IPCM_AVAILABLE

#include <cassert>

namespace exahype {
namespace profilers {
namespace ipcm {

void IpcmCountMetric::setNumberOfTags(int n) { counts_.reserve(n); }

void IpcmCountMetric::registerTag(const std::string& tag) {
  assert((counts_.count(tag) == 0) &&
         "At least one tag has been registered more than once");
  counts_[tag];
}

void IpcmCountMetric::stop(const std::string& tag,
                           const SystemCounterState& before_state,
                           const SystemCounterState& after_state) {
  counts_[tag]++;
}

void IpcmCountMetric::writeToOstream(std::ostream* os) const {
  for (const auto& count : counts_) {
    *os << "Count " << count.first << " " << count.second << std::endl;
  }
}

}  // namespace ipcm
}  // namespace profilers
}  // namespace exahype

#endif  // IPCM_AVAILABLE
