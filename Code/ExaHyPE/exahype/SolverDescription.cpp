#include "exahype/SolverDescription.h"


std::vector<exahype::SolverDescription>  exahype::SolverDescription::ExistingSolvers;


exahype::SolverDescription::SolverDescription(const std::string& identifier, int kernelNumber):
  _identifier(identifier),
  _kernelNumber(kernelNumber) {
}
