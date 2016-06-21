#include "ProfilerFactory.h"

#include <cstdlib>
#include <functional>
#include <iostream>
#include <unordered_map>

#include "simple/ChronoElapsedTimeProfiler.h"
#include "simple/NoOpProfiler.h"

#ifdef LIKWID_AVAILABLE
#include "likwid/LikwidProfiler.h"
#include "likwid/modules/LikwidPerformanceMonitoringModule.h"
#include "likwid/modules/LikwidPowerAndEnergyMonitoringModule.h"
#include "likwid/modules/LikwidTimeMeasurementModule.h"
#endif  // LIKWID_AVAILABLE

#ifdef IPCM_AVAILABLE
#include "ipcm/IpcmProfiler.h"
#include "ipcm/metrics/IpcmBytesReadDramMetric.h"
#include "ipcm/metrics/IpcmBytesWrittenDramMetric.h"
#include "ipcm/metrics/IpcmConsumedJoulesMetric.h"
#include "ipcm/metrics/IpcmCyclesMetric.h"
#endif  // IPCM_AVAILABLE

namespace {

#ifdef LIKWID_AVAILABLE
const std::unordered_map<
    std::string,
    std::function<std::unique_ptr<exahype::profilers::likwid::LikwidModule>(
        const exahype::profilers::likwid::LikwidProfilerState&,
        const std::string&)>>
    likwid_module_map = {
        {"LikwidTimeMeasurementModule",
         [](const exahype::profilers::likwid::LikwidProfilerState& state,
            const std::string& suffix) {
           return std::unique_ptr<
               exahype::profilers::likwid::LikwidTimeMeasurementModule>(
               new exahype::profilers::likwid::LikwidTimeMeasurementModule(
                   state));
         }},
        {
            "LikwidPowerAndEnergyMonitoringModule",
            [](const exahype::profilers::likwid::LikwidProfilerState& state,
               const std::string& suffix) {
              return std::unique_ptr<exahype::profilers::likwid::
                                         LikwidPowerAndEnergyMonitoringModule>(
                  new exahype::profilers::likwid::
                      LikwidPowerAndEnergyMonitoringModule(state));
            },
        },
        {"LikwidPerformanceMonitoringModule",
         [](const exahype::profilers::likwid::LikwidProfilerState& state,
            const std::string& suffix) {
           return std::unique_ptr<
               exahype::profilers::likwid::LikwidPerformanceMonitoringModule>(
               new exahype::profilers::likwid::
                   LikwidPerformanceMonitoringModule(state, suffix));
         }}};
#endif  // LIKWID_AVAILABLE

#ifdef IPCM_AVAILABLE
const std::unordered_map<
    std::string,
    std::function<std::unique_ptr<exahype::profilers::ipcm::IpcmMetric>()>>
    ipcm_metrics_map = {
        {"IpcmCyclesMetric",
         []() {
           return std::unique_ptr<exahype::profilers::ipcm::IpcmMetric>(
               new exahype::profilers::ipcm::IpcmCyclesMetric);
         }},
        {"IpcmBytesReadDramMetric",
         []() {
           return std::unique_ptr<
               exahype::profilers::ipcm::IpcmBytesReadDramMetric>(
               new exahype::profilers::ipcm::IpcmBytesReadDramMetric);
         }},
        {"IpcmBytesWrittenDramMetric",
         []() {
           return std::unique_ptr<
               exahype::profilers::ipcm::IpcmBytesWrittenDramMetric>(
               new exahype::profilers::ipcm::IpcmBytesWrittenDramMetric);
         }},
        {"IpcmConsumedJoulesMetric",
         []() {
           return std::unique_ptr<
               exahype::profilers::ipcm::IpcmConsumedJoulesMetric>(
               new exahype::profilers::ipcm::IpcmConsumedJoulesMetric);
         }

        },
};
#endif  // IPCM_AVAILABLE

const std::unordered_map<
    std::string, std::function<std::unique_ptr<exahype::profilers::Profiler>(
                     const std::vector<std::string>&)>>
    profiler_map = {
        {"NoOpProfiler",
         [](const std::vector<std::string>& metrics) {
           return std::unique_ptr<exahype::profilers::simple::NoOpProfiler>(
               new exahype::profilers::simple::NoOpProfiler());
         }},
        {"ChronoElapsedTimeProfiler",
         [](const std::vector<std::string>& metrics) {
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
             std::string module_identifier, suffix;
             size_t pos_underscore = module.find_first_of('_');
             if (pos_underscore != std::string::npos) {
               module_identifier = module.substr(0, pos_underscore);
               suffix = module.substr(pos_underscore + 1);
             } else {
               module_identifier = module;
               suffix = "";
             }
             if (likwid_module_map.count(module_identifier)) {
               profiler->addModule(likwid_module_map.at(module_identifier)(
                   profiler->state(), suffix));
             } else {
               std::cerr << "ProfilerFactory: Unknown likwid module name '"
                         << module_identifier << "'." << std::endl;
               std::exit(EXIT_FAILURE);
             }
           }
           return profiler;
         }},
#endif  // LIKWID_AVAILABLE
#ifdef IPCM_AVAILABLE
        {"IpcmProfiler",
         [](const std::vector<std::string>& metrics) {
           std::unique_ptr<exahype::profilers::ipcm::IpcmProfiler> profiler(
               new exahype::profilers::ipcm::IpcmProfiler());
           for (const auto& metric : metrics) {
             if (ipcm_metrics_map.count(metric)) {
               profiler->addMetric(ipcm_metrics_map.at(metric)());
             } else {
               std::cerr << "ProfilerFactory: Unknown ipcm metric name '"
                         << metric << "'" << std::endl;
               std::exit(EXIT_FAILURE);
             }
           }
           return profiler;
         }},
#endif  // IPCM_AVAILABLE
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
