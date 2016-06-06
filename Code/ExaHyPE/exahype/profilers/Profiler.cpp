#include "Profiler.h"

#include <fstream>

namespace exahype {
namespace profilers {

void Profiler::writeToCout() const { writeToOstream(&std::cout); }

void Profiler::writeToFile(const std::string& path) const {
  std::ofstream ofs(path);

  if (ofs.is_open()) {
    writeToOstream(&ofs);
  } else {
    std::cerr << "Profiler: Could not write to file '" << path << "'"
              << std::endl;
  }
}

}  // namespace profilers
}  // namespace exahype
