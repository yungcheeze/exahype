/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released unter the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#include "ChronoElapsedTimeProfiler.h"

namespace exahype {
namespace profilers {
namespace simple {

void ChronoElapsedTimeProfiler::setNumberOfTags(int n) {
  time_points_.reserve(n);
  counts_and_durations_.reserve(n);
}

void ChronoElapsedTimeProfiler::registerTag(const std::string& tag) {
  time_points_[tag];
  counts_and_durations_[tag];
}

void ChronoElapsedTimeProfiler::start(const std::string& tag) {
  time_points_[tag] = std::chrono::high_resolution_clock::now();
}

void ChronoElapsedTimeProfiler::stop(const std::string& tag) {
  auto end = std::chrono::high_resolution_clock::now();
  auto start = time_points_[tag];
  auto& pair = counts_and_durations_[tag];
  pair.first++;                  // count
  pair.second += (end - start);  // total elapsed time
}

void ChronoElapsedTimeProfiler::writeToOstream(std::ostream* os) const {
  *os << "ChronoElapsedTimeProfiler {" << std::endl;

  for (const auto& kv_pair : counts_and_durations_) {
    *os << "  " << kv_pair.first << " {" << std::endl;
    *os << "    count = " << kv_pair.second.first << std::endl;
    *os << "    time_sec = "
        << static_cast<std::chrono::duration<double, std::ratio<1>>>(
               kv_pair.second.second)
               .count()
        << std::endl;
    *os << "  }" << std::endl;
  }

  *os << "}" << std::endl;
}

}  // namespace simple
}  // namespace profilers
}  // namespace exahype
