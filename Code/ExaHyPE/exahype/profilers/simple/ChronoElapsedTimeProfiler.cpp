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

#include "ChronoElapsedTimeProfiler.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <thread>

#include "../ProfilerUtils.h"

namespace {
using namespace exahype::profilers::utils;

const uint64_t kNumberOfSamples1 = 1000000;
const uint64_t kNumberOfSamples2 = 10;

std::unordered_map<std::string,
                   std::chrono::time_point<std::chrono::steady_clock>>
    overhead_time_points;

void initOverheadMeasurementMap() {
  static const int kNumberOfKeys = 12;
  std::array<std::string, kNumberOfKeys> keys = {"boundaryConditions",
                                                 "volumeUnknownsRestriction",
                                                 "riemannSolver",
                                                 "volumeIntegral",
                                                 "surfaceIntegral",
                                                 "solutionUpdate",
                                                 "spaceTimePredictor",
                                                 "stableTimeStepSize",
                                                 "faceUnknownsRestriction",
                                                 "solutionAdjustment",
                                                 "faceUnknownsProlongation",
                                                 "volumeUnknownsProlongation"};
  overhead_time_points.reserve(kNumberOfKeys);
  std::for_each(keys.begin(), keys.end(),
                [](const std::string& key) { overhead_time_points[key]; });
}

void storeInMap() {
  overhead_time_points["solutionUpdate"] = std::chrono::steady_clock::now();
  escape(&overhead_time_points);
}

void getCurrentTime() {
  auto now = std::chrono::steady_clock::now();
  escape(&now);
}

static void estimateOverhead(const std::string& output) {
  if (output == "" || output == "cout") {
    return;
  }

  std::ofstream file;
  file.open(output, std::ios::out | std::ios::trunc);

  file << "steady_clock::period = "
       << static_cast<double>(std::chrono::steady_clock::period::num) /
              std::chrono::steady_clock::period::den
       << "sec" << std::endl;

  // noop
  std::tuple<mean_sec, median_sec, std_sec, min_sec, max_sec> duration_noop =
      meanMedianStdMinMaxOfDurations<kNumberOfSamples1>(
          nTimesDurationOf<kNumberOfSamples1, &noop>());
  file << "duration_noop" << std::endl;
  file << "  mean_sec = " << std::get<0>(duration_noop) << std::endl;
  file << "  median_sec = " << std::get<1>(duration_noop) << std::endl;
  file << "  std_sec = " << std::get<2>(duration_noop) << std::endl;
  file << "  min_sec = " << std::get<3>(duration_noop) << std::endl;
  file << "  max_sec = " << std::get<4>(duration_noop) << std::endl;

  // storeInMap
  initOverheadMeasurementMap();
  std::tuple<mean_sec, median_sec, std_sec, min_sec, max_sec>
      duration_storeInMap = meanMedianStdMinMaxOfDurations<kNumberOfSamples1>(
          nTimesDurationOf<kNumberOfSamples1, &storeInMap>());
  file << "duration_storeInMap" << std::endl;
  file << "  mean_sec = " << std::get<0>(duration_storeInMap) << std::endl;
  file << "  median_sec = " << std::get<1>(duration_storeInMap) << std::endl;
  file << "  std_sec = " << std::get<2>(duration_storeInMap) << std::endl;
  file << "  min_sec = " << std::get<3>(duration_storeInMap) << std::endl;
  file << "  max_sec = " << std::get<4>(duration_storeInMap) << std::endl;

  // getCurrentTime
  std::tuple<mean_sec, median_sec, std_sec, min_sec, max_sec>
      duration_getCurrentTime =
          meanMedianStdMinMaxOfDurations<kNumberOfSamples1>(
              nTimesDurationOf<kNumberOfSamples1, &getCurrentTime>());
  file << "getCurrentTime" << std::endl;
  file << "  mean_sec = " << std::get<0>(duration_getCurrentTime) << std::endl;
  file << "  median_sec = " << std::get<1>(duration_getCurrentTime)
       << std::endl;
  file << "  std_sec = " << std::get<2>(duration_getCurrentTime) << std::endl;
  file << "  min_sec = " << std::get<3>(duration_getCurrentTime) << std::endl;
  file << "  max_sec = " << std::get<4>(duration_getCurrentTime) << std::endl;

  /*
  // sleep 1s
  std::cout << "sleep 1s" << std::endl;
  std::tuple<mean_sec, median_sec, std_sec, min_sec, max_sec> duration_sleep1s =
      meanMedianStdMinMaxOf<kNumberOfSamples2>(
          nTimesDurationOf<kNumberOfSamples2, &sleep<1000>>());
  file << "duration_sleep1s" << std::endl;
  file << "  mean_sec = " << std::get<0>(duration_sleep1s) << std::endl;
  file << "  median_sec = " << std::get<1>(duration_sleep1s) << std::endl;
  file << "  std_sec = " << std::get<2>(duration_sleep1s) << std::endl;
  file << "  min_sec = " << std::get<3>(duration_sleep1s) << std::endl;
  file << "  max_sec = " << std::get<4>(duration_sleep1s) << std::endl;

  // sleep 2s
  std::cout << "sleep 2s" << std::endl;
  std::tuple<mean_sec, median_sec, std_sec, min_sec, max_sec> duration_sleep2s =
      meanMedianStdMinMaxOf<kNumberOfSamples2>(
          nTimesDurationOf<kNumberOfSamples2, &sleep<2000>>());
  file << "duration_sleep2s" << std::endl;
  file << "  mean_sec = " << std::get<0>(duration_sleep2s) << std::endl;
  file << "  median_sec = " << std::get<1>(duration_sleep2s) << std::endl;
  file << "  std_sec = " << std::get<2>(duration_sleep2s) << std::endl;
  file << "  min_sec = " << std::get<3>(duration_sleep2s) << std::endl;
  file << "  max_sec = " << std::get<4>(duration_sleep2s) << std::endl;
  */
}

}  // namespace

namespace exahype {
namespace profilers {
namespace simple {

ChronoElapsedTimeProfiler::ChronoElapsedTimeProfiler(const std::string& output)
    : Profiler(output) {
  estimateOverhead(output);
}

void ChronoElapsedTimeProfiler::setNumberOfTags(int n) {
  time_points_.reserve(n);
  counts_and_durations_.reserve(n);

  //  individual_measurements_ns_.reserve(n);
}

void ChronoElapsedTimeProfiler::registerTag(const std::string& tag) {
  time_points_[tag];
  counts_and_durations_[tag];

  //  individual_measurements_ns_[tag];
}

void ChronoElapsedTimeProfiler::start(const std::string& tag) {
  time_points_[tag] = std::chrono::steady_clock::now();
}

void ChronoElapsedTimeProfiler::stop(const std::string& tag) {
  auto end = std::chrono::steady_clock::now();
  escape(&end);

  auto start = time_points_[tag];
  auto& pair = counts_and_durations_[tag];
  pair.first++;                  // count
  pair.second += (end - start);  // total elapsed time

  //  individual_measurements_ns_[tag].push_back((end - start).count());
}

void ChronoElapsedTimeProfiler::writeToOstream(std::ostream* os) const {
  *os << "ChronoElapsedTimeProfiler" << std::endl;

  for (const auto& kv_pair : counts_and_durations_) {
    *os << "  " << kv_pair.first << std::endl;
    *os << "    count = " << kv_pair.second.first << std::endl;
    *os << "    time_sec = "
        << static_cast<std::chrono::duration<double, std::ratio<1>>>(
               kv_pair.second.second)
               .count()
        << std::endl;
  }

  /*
  *os << "individual measurements" << std::endl;
  for (const auto& kv_pair : individual_measurements_ns_) {
    *os << "  " << kv_pair.first << " = [" << std::endl;
    for (const auto& value : kv_pair.second) {
      *os << "    " << value << std::endl;
    }
    *os << "  ];" << std::endl;
  }
  */
}

}  // namespace simple
}  // namespace profilers
}  // namespace exahype
