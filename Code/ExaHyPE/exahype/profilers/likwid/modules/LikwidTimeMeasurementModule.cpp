#include "LikwidTimeMeasurementModule.h"

#ifdef LIKWID_AVAILABLE

#include <cassert>
#include <utility>

namespace exahype {
namespace profilers {
namespace likwid {

LikwidTimeMeasurementModule::LikwidTimeMeasurementModule(
    const LikwidProfilerState& state)
    : LikwidModule(state) {
  timer_init();
}

LikwidTimeMeasurementModule::~LikwidTimeMeasurementModule() {
  timer_finalize();
}

void LikwidTimeMeasurementModule::setNumberOfTags(int n) {
  timer_data_.reserve(n);
}

void LikwidTimeMeasurementModule::registerTag(const std::string& tag) {
  assert((timer_data_.count(tag) == 0) &&
         "At least one tag has been registered twice");
  timer_data_[tag];
}

void LikwidTimeMeasurementModule::start(const std::string& tag) {
  assert((timer_data_.count(tag) == 1) &&
         "At least one tag has not been registered prior to starting the "
         "corresponding measurement");
  timer_start(&timer_data_.at(tag));
}

void LikwidTimeMeasurementModule::stop(const std::string& tag) {
  timer_stop(&timer_data_.at(tag));
}

void LikwidTimeMeasurementModule::writeToOstream(std::ostream* os) const {
  *os << "TimeMeasurementModule" << std::endl;

  for (const auto& m : timer_data_) {
    *os << m.first
        << ": time_sec = " << timer_print(const_cast<TimerData*>(&m.second))
        << std::endl;
    *os << m.first
        << ": cycles = " << timer_printCycles(const_cast<TimerData*>(&m.second))
        << std::endl;
  }
}

}  // namespace likwid
}  // namespace profilers
}  // namespace exahype

#endif  // LIKWID_AVAILABLE
