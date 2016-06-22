/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
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
  for (const auto& m : timer_data_) {
    *os << "TimeMeasurementModule " << m.first << " time_sec "
        << timer_print(const_cast<TimerData*>(&m.second)) << std::endl;
    *os << "TimeMeasurementModule " << m.first << " cycles "
        << timer_printCycles(const_cast<TimerData*>(&m.second)) << std::endl;
  }
}

}  // namespace likwid
}  // namespace profilers
}  // namespace exahype

#endif  // LIKWID_AVAILABLE
