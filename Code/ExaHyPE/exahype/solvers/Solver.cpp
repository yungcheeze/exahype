#include "exahype/solvers/Solver.h"


//std::vector<exahype::solvers::Solver>  exahype::solvers::Solver::RegisteredSolvers;
std::vector<exahype::solvers::Solver*>  exahype::solvers::RegisteredSolvers;


exahype::solvers::Solver::Solver(const std::string& identifier, int kernelNumber):
  _identifier(identifier),
  _kernelNumber(kernelNumber) {
}
