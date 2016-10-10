/*
 * SingleSolverCoupling.cpp
 *
 *  Created on: 6 Oct 2016
 *      Author: dominic
 */

#include "SingleSolverCoupling.h"

#include "tarch/Assertions.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "kernels/limiter/generic/Limiter.h"


exahype::solvers::SingleSolverCoupling::SingleSolverCoupling(
    int solverNumber)
: CellWiseCoupling(0.0,0.0),
  _solverNumber(solverNumber)
{
  assertion2(_solverNumber >= 0 && _solverNumber<static_cast<int>(exahype::solvers::RegisteredSolvers.size()),
      _solverNumber,exahype::solvers::RegisteredSolvers.size());
}

void exahype::solvers::SingleSolverCoupling::coupleFirstTime(
    const int cellDescriptionsIndex,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[_solverNumber];
  int element = solver->tryGetElement(cellDescriptionsIndex,_solverNumber);

  if (element!=exahype::solvers::Solver::NotFound) {
    std::cout << "Hallo!";
    solver->synchroniseTimeStepping(cellDescriptionsIndex,element);

    solver->setInitialConditions(
        cellDescriptionsIndex,element,
        fineGridVertices,fineGridVerticesEnumerator);
  }
}

void exahype::solvers::SingleSolverCoupling::couple(
    const int cellDescriptionsIndex,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[_solverNumber];
  int element = solver->tryGetElement(cellDescriptionsIndex,_solverNumber);

  if (element!=exahype::solvers::Solver::NotFound) {
    solver->updateSolution(cellDescriptionsIndex,element,fineGridVertices,fineGridVerticesEnumerator);
  }
}
