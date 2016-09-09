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

#include "../../ProfilerUtils.h"
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
  aggregates_cycles_seconds_.reserve(n);
}

void LikwidTimeMeasurementModule::registerTag(const std::string& tag) {
  assert((aggregates_cycles_seconds_.count(tag) == 0) &&
         "At least one tag has been registered twice");
  aggregates_cycles_seconds_[tag] = {0, 0.0};
}

void LikwidTimeMeasurementModule::start(const std::string& tag) {
  assert((aggregates_cycles_seconds_.count(tag) == 1) &&
         "At least one tag has not been registered prior to starting the "
         "corresponding measurement");
  timer_start(&timer_data_);
}

void LikwidTimeMeasurementModule::stop(const std::string& tag) {
  timer_stop(&timer_data_);
  utils::escape(&timer_data_);
  auto& cycles_seconds_pair = aggregates_cycles_seconds_.at(tag);
  cycles_seconds_pair.first += timer_printCycles(&timer_data_);
  cycles_seconds_pair.second += timer_print(&timer_data_);
}

void LikwidTimeMeasurementModule::writeToOstream(std::ostream* os) const {
  for (const auto& m : aggregates_cycles_seconds_) {
    *os << "TimeMeasurementModule " << m.first << " cycles " << m.second.first
        << std::endl;
    *os << "TimeMeasurementModule " << m.first << " time_sec "
        << m.second.second << std::endl;
  }
}

}  // namespace likwid
}  // namespace profilers
}  // namespace exahype

#endif  // LIKWID_AVAILABLE
