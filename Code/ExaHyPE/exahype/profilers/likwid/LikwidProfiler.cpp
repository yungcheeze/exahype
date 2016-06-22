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
 
#include "LikwidProfiler.h"

#include <algorithm>
#include <cstdlib>

#ifdef LIKWID_AVAILABLE

namespace exahype {
namespace profilers {
namespace likwid {

LikwidProfiler::LikwidProfiler(const std::vector<int>& cpus /* = {} */) {
  int err;

  setenv("LIKWID_FORCE", "1", 1);

  topology_init();
  state_.cpi_info_ = get_cpuInfo();
  state_.cpu_topology_ = get_cpuTopology();

  err = numa_init();
  if (err != 0) {
    std::cerr << "LikwidProfiler: numa_init() returned " << err << " != 0"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  state_.numa_topology_ = get_numaTopology();
  affinity_init();
  state_.affinity_domains_ = get_affinityDomains();

  // Initialize CPU ID vector
  if (cpus.empty()) {
    // Add all processors to list
    for (unsigned int i = 0; i < state_.cpu_topology_->numHWThreads; i++) {
      state_.cpus_.push_back(state_.cpu_topology_->threadPool[i].apicId);
    }
  } else {
    // Only add CPU IDs specified by constructor argument
    state_.cpus_ = cpus;
  }

  // TODO(guera): Remove this hack
  affinity_pinThread(0);
}

LikwidProfiler::~LikwidProfiler() {
  modules_.clear();  // ensure correct order of deinitialization

  affinity_finalize();
  numa_finalize();
  topology_finalize();
}

void LikwidProfiler::addModule(std::unique_ptr<LikwidModule> module) {
  modules_.emplace_back(std::move(module));
}

void LikwidProfiler::setNumberOfTags(int n) {
  for (auto& module : modules_) {
    module->setNumberOfTags(n);
  }
}

void LikwidProfiler::registerTag(const std::string& tag) {
  for (auto& module : modules_) {
    module->registerTag(tag);
  }
}

void LikwidProfiler::start(const std::string& tag) {
  for (auto& module : modules_) {
    module->start(tag);
  }
}

void LikwidProfiler::stop(const std::string& tag) {
  for (auto& module : modules_) {
    module->stop(tag);
  }
}

void LikwidProfiler::writeToOstream(std::ostream* os) const {
  for (const auto& module : modules_) {
    module->writeToOstream(os);
  }
}

}  // namespace likwid
}  // namespace profilers
}  // namespace exahype

#endif  // LIKWID_AVAILABLE
