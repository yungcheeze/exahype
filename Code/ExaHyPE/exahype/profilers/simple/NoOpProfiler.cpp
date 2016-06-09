#include "NoOpProfiler.h"

void exahype::profilers::simple::NoOpProfiler::setNumberOfTags(int n) {}

void exahype::profilers::simple::NoOpProfiler::registerTag(
    const std::string& tag) {}

void exahype::profilers::simple::NoOpProfiler::start(const std::string& tag) {}

void exahype::profilers::simple::NoOpProfiler::stop(const std::string& tag) {}

void exahype::profilers::simple::NoOpProfiler::writeToOstream(
    std::ostream* os) const {}
