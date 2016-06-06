#include "ProfilerFactory.h"

#include <functional>
#include <iostream>
#include <unordered_map>

#include "simple/ChronoElapsedTimeProfiler.h"
#include "simple/NoOpProfiler.h"

namespace {

const std::unordered_map<
    std::string, std::function<std::unique_ptr<exahype::profilers::Profiler>(
                     const std::vector<std::string>&)>>
    map = {
        {"NoOpProfiler",
         [](const std::vector<std::string>& modules) {
           return std::make_unique<exahype::profilers::simple::NoOpProfiler>();
         }},
        {"ChronoElapsedTimeProfiler",
         [](const std::vector<std::string>& modules) {
           return std::make_unique<
               exahype::profilers::simple::ChronoElapsedTimeProfiler>();
         }},
};

}  // namespace

namespace exahype {
namespace profilers {

ProfilerFactory& ProfilerFactory::getInstance() {
  static ProfilerFactory singleton;
  return singleton;
}

std::unique_ptr<Profiler> ProfilerFactory::create(
    const std::string& profiler_name, const std::vector<std::string>& modules) {
  if (map.count(profiler_name)) {  // known profiler
    return map.at(profiler_name)(modules);
  } else {  // unknown profiler
    std::cerr << "ProfilerFactory: Unknown profiler name '" << profiler_name
              << "'. NoOpProfiler created instead." << std::endl;
    return map.at("NoOpProfiler")(modules);
  }
}

}  // namespace profilers
}  // namespace exahype
