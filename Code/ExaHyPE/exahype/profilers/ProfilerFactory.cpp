#include "ProfilerFactory.h"

#include <cstdlib>
#include <functional>
#include <iostream>
#include <unordered_map>

#include "simple/ChronoElapsedTimeProfiler.h"
#include "simple/NoOpProfiler.h"

#ifdef LIKWID_AVAILABLE
#include "likwid/LikwidProfiler.h"
#include "likwid/modules/LikwidTimeMeasurementModule.h"
#endif  // LIKWID_AVAILABLE

namespace {

#ifdef LIKWID_AVAILABLE
const std::unordered_map<
    std::string,
    std::function<std::unique_ptr<exahype::profilers::likwid::LikwidModule>(
        const exahype::profilers::likwid::LikwidProfilerState&)>>
    likwid_module_map = {
        {"LikwidTimeMeasurementModule",
         [](const exahype::profilers::likwid::LikwidProfilerState& state) {
           return std::unique_ptr<
               exahype::profilers::likwid::LikwidTimeMeasurementModule>(
               new exahype::profilers::likwid::LikwidTimeMeasurementModule(
                   state));
         }},
};
#endif  // LIKWID_AVAILABLE

const std::unordered_map<
    std::string, std::function<std::unique_ptr<exahype::profilers::Profiler>(
                     const std::vector<std::string>&)>>
    profiler_map = {
        {"NoOpProfiler",
         [](const std::vector<std::string>& modules) {
           return std::unique_ptr<exahype::profilers::simple::NoOpProfiler>(
               new exahype::profilers::simple::NoOpProfiler());
         }},
        {"ChronoElapsedTimeProfiler",
         [](const std::vector<std::string>& modules) {
           return std::unique_ptr<
               exahype::profilers::simple::ChronoElapsedTimeProfiler>(
               new exahype::profilers::simple::ChronoElapsedTimeProfiler());
         }},
#ifdef LIKWID_AVAILABLE
        {"LikwidProfiler",
         [](const std::vector<std::string>& modules) {
           std::unique_ptr<exahype::profilers::likwid::LikwidProfiler> profiler(
               new exahype::profilers::likwid::LikwidProfiler());
           for (const auto& module : modules) {
             if (likwid_module_map.count(module)) {
               profiler->addModule(
                   likwid_module_map.at(module)(profiler->state()));
             } else {
               std::cerr << "ProfilerFactory: Unknown likwid module name '"
                         << module << "'." << std::endl;
               std::exit(EXIT_FAILURE);
             }
           }
           return profiler;
         }},
#endif  // LIKWID_AVAILABLE
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
  if (profiler_map.count(profiler_name)) {  // known profiler
    return profiler_map.at(profiler_name)(modules);
  } else {  // unknown profiler
    std::cerr << "ProfilerFactory: Unknown profiler name '" << profiler_name
              << "'. NoOpProfiler created instead." << std::endl;
    return profiler_map.at("NoOpProfiler")(modules);
  }
}

}  // namespace profilers
}  // namespace exahype
