#include "LikwidPowerAndEnergyMonitoringModule.h"

#ifdef LIKWID_AVAILABLE

#include "../LikwidProfiler.h"
#include <cassert>
#include <cstdlib>

namespace exahype {
namespace profilers {
namespace likwid {

constexpr PowerType LikwidPowerAndEnergyMonitoringModule::kProfiledPowerTypes
    [LikwidPowerAndEnergyMonitoringModule::kNumberOfProfiledPowerTypes];

LikwidPowerAndEnergyMonitoringModule::LikwidPowerAndEnergyMonitoringModule(
    const LikwidProfilerState& state)
    : LikwidModule(state) {
  int has_rapl = power_init(0);  // test for CPU 0
  if (has_rapl == 0) {
    std::cerr << "LikwidPowerAndEnergyMonitoringModule: has_rapl == 0"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const AffinityDomain* affinity_domains_begin =
      state_.affinity_domains_->domains;
  const AffinityDomain* affinity_domains_end =
      affinity_domains_begin +
      state_.affinity_domains_->numberOfAffinityDomains;

  // Find all sockets
  number_of_sockets_ = 0;
  for (const AffinityDomain* current_affinity_domain = affinity_domains_begin;
       current_affinity_domain != affinity_domains_end;
       current_affinity_domain++) {
    const std::string tag = std::string(current_affinity_domain->tag->data,
                                        current_affinity_domain->tag->data +
                                            current_affinity_domain->tag->slen);

    // Currently possible values: N (node), SX (socket/package X), CX (cache
    // domain X) and MX (memory domain X)
    size_t pos = tag.find("S");
    if (pos != std::string::npos) {  // current affinity domain is socket
      number_of_sockets_++;

      int socket_id = std::stoi(tag.substr(pos + 1));  // Parse X in SX
      assert(socket_id == static_cast<int>(cpuid_of_sockets_.size()) &&
             "Sockets are not numbered consecutively");

      int cpu_id = current_affinity_domain->processorList[0];
      cpuid_of_sockets_.push_back(cpu_id);

      HPMaddThread(cpu_id);
      // power_init(cpu_id);
    }
  }

  power_info_ = get_powerInfo();
}

LikwidPowerAndEnergyMonitoringModule::~LikwidPowerAndEnergyMonitoringModule() {
  power_finalize();
  power_info_ = nullptr;
}

void LikwidPowerAndEnergyMonitoringModule::setNumberOfTags(int n) {
  power_data_.reserve(n);
  aggregates_.reserve(n);
}

void LikwidPowerAndEnergyMonitoringModule::registerTag(const std::string& tag) {
  assert((power_data_.count(tag) == 0) &&
         "At least one tag has been registered twice");
  power_data_[tag].resize(number_of_sockets_);

  auto& vec_array_value = aggregates_[tag];
  vec_array_value.resize(number_of_sockets_);

  // For all sockets...
  for (auto& array_value : vec_array_value) {
    // initialize array of aggregates to zero
    std::fill(array_value.begin(), array_value.end(), 0.0);
  }
}

void LikwidPowerAndEnergyMonitoringModule::start(const std::string& tag) {
  assert(power_data_.count(tag) && "Unregistered tag encountered");

  // For all sockets...
  auto& vec_array_power_data = power_data_[tag];

  for (int i = 0; i < number_of_sockets_; i++) {
    // for all profiled power types...
    for (int j = 0; j < kNumberOfProfiledPowerTypes; j++) {
      // start power measurement
      power_start(&vec_array_power_data[i][j], cpuid_of_sockets_[i],
                  kProfiledPowerTypes[i]);
    }
  }
}

void LikwidPowerAndEnergyMonitoringModule::stop(const std::string& tag) {
  assert(power_data_.count(tag) && "Unregistered tag encountered");

  // Stop all measurements first
  auto& vec_array_power_data = power_data_[tag];

  // For all sockets
  for (int i = 0; i < number_of_sockets_; i++) {
    // for all profiled power types
    for (int j = 0; j < kNumberOfProfiledPowerTypes; j++) {
      // stop power measurement
      power_stop(&vec_array_power_data[i][j], cpuid_of_sockets_[i],
                 kProfiledPowerTypes[j]);
    }
  }

  // then update all aggregates
  auto& vec_array_value = aggregates_[tag];

  // For all sockets
  for (int i = 0; i < number_of_sockets_; i++) {
    // for all power types
    for (int j = 0; j < kNumberOfProfiledPowerTypes; j++) {
      vec_array_value[i][j] += power_printEnergy(&vec_array_power_data[i][j]);
    }
  }
}

void LikwidPowerAndEnergyMonitoringModule::writeToOstream(
    std::ostream* os) const {
  // For all tags
  for (const auto& tag_vec_array_value_pair : aggregates_) {
    // for all sockets
    for (int i = 0;
         i < static_cast<int>(tag_vec_array_value_pair.second.size()); i++) {
      // for all profiled power types
      for (int j = 0; j < kNumberOfProfiledPowerTypes; j++) {
        *os << "PowerAndEnergyMonitoringModule "
            << tag_vec_array_value_pair.first << " socket" << i << " "
            << powerTypeToString(kProfiledPowerTypes[j]) << " "
            << tag_vec_array_value_pair.second[i][j] << std::endl;
      }
    }
  }
}

}  // namespace likwid
}  // namespace profilers
}  // namespace exahype

#endif  // LIKWID_AVAILABLE
