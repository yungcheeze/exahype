#include "exahype/solvers/Solver.h"


//std::vector<exahype::solvers::Solver>  exahype::solvers::Solver::RegisteredSolvers;
std::vector<exahype::solvers::Solver*>  exahype::solvers::RegisteredSolvers;


exahype::solvers::Solver::Solver(const std::string& identifier, Type type, int kernelNumber, int numberOfVariables, int nodesPerCoordinateAxis):
  _identifier(identifier),
  _type(type),
  _kernelNumber(kernelNumber),
  _numberOfVariables(numberOfVariables),
  _nodesPerCoordinateAxis(nodesPerCoordinateAxis) {
}


exahype::solvers::Solver::Type exahype::solvers::Solver::getType() const {
  return _type;
}


int exahype::solvers::Solver::getNumberOfVariables() const {
  return _numberOfVariables;
}


int exahype::solvers::Solver::getNodesPerCoordinateAxis() const {
  return _nodesPerCoordinateAxis;
}

