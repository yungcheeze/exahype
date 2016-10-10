/*
 * ADERDGAPosterioriSubcellLimiter.cpp
 *
 *  Created on: 6 Oct 2016
 *      Author: dominic
 */

#include "ADERDGAPosterioriSubcellLimiter.h"

#include "tarch/Assertions.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "kernels/limiter/generic/Limiter.h"


exahype::solvers::ADERDGAPosterioriSubcellLimiter::ADERDGAPosterioriSubcellLimiter(
    int aderdgSolverNumber,int finiteVolumesSolverNumber)
      : CellWiseCoupling(0.0,0.0),
        _aderdgSolverNumber(aderdgSolverNumber),
        _finiteVolumesSolverNumber(finiteVolumesSolverNumber)
{
  assertion2(_aderdgSolverNumber >= 0 && _aderdgSolverNumber<static_cast<int>(exahype::solvers::RegisteredSolvers.size()),
      _aderdgSolverNumber,exahype::solvers::RegisteredSolvers.size());
  assertion2(_finiteVolumesSolverNumber >= 0 && _finiteVolumesSolverNumber<static_cast<int>(exahype::solvers::RegisteredSolvers.size()),
      _finiteVolumesSolverNumber,exahype::solvers::RegisteredSolvers.size());
  assertion1(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]->getType()==exahype::solvers::Solver::Type::ADER_DG,
      exahype::solvers::Solver::toString(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]->getType()));
  assertion1(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]->getType()==exahype::solvers::Solver::Type::FiniteVolumes,
        exahype::solvers::Solver::toString(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]->getType()));
}

void exahype::solvers::ADERDGAPosterioriSubcellLimiter::coupleFirstTime(
    const int cellDescriptionsIndex,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  exahype::solvers::ADERDGSolver* aderdgSolver =
      static_cast<exahype::solvers::ADERDGSolver*>(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]);
  exahype::solvers::FiniteVolumesSolver* finiteVolumesSolver =
      static_cast<exahype::solvers::FiniteVolumesSolver*>(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]);
  int aderdgElement        = aderdgSolver->tryGetElement(cellDescriptionsIndex,_aderdgSolverNumber);
  int finiteVolumesElement = finiteVolumesSolver->tryGetElement(cellDescriptionsIndex,_finiteVolumesSolverNumber);

  if (aderdgElement!=exahype::solvers::Solver::NotFound) {
    auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
        cellDescriptionsIndex)[aderdgElement];
    assertion2(finiteVolumesElement!=exahype::solvers::Solver::NotFound,cellDescriptionsIndex,_aderdgSolverNumber);
    auto& finiteVolumesCellDescription = exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(
        cellDescriptionsIndex)[finiteVolumesElement];

    // set initial conditions on both solvers
    // TODO(Dominic): Does it make sense to project initial conditions onto
    // FV solution and then on ADER-DG solution?
    // Limiting might consider initial condition with oscillations
    // as the correct one.
    aderdgSolver->synchroniseTimeStepping(cellDescriptionsIndex,aderdgElement);
    aderdgSolver->setInitialConditions(
        cellDescriptionsIndex,aderdgElement,
        fineGridVertices,fineGridVerticesEnumerator);

    finiteVolumesSolver->synchroniseTimeStepping(cellDescriptionsIndex,finiteVolumesElement);
    finiteVolumesSolver->setInitialConditions(
        cellDescriptionsIndex,finiteVolumesElement,
        fineGridVertices,fineGridVerticesEnumerator);

    double* aderdgSolution  = exahype::solvers::ADERDGSolver::DataHeap::getInstance().
        getData(aderdgCellDescription.getSolution()).data();
    double* minOfNeighbours =  exahype::solvers::ADERDGSolver::DataHeap::getInstance().
        getData(aderdgCellDescription.getSolutionMin()).data();
    double* maxOfNeighbours =  exahype::solvers::ADERDGSolver::DataHeap::getInstance().
        getData(aderdgCellDescription.getSolutionMax()).data();

    bool initialSolutionNeedsLimiting =
        kernels::limiter::generic::c::isTroubledCell(
            aderdgSolution,aderdgSolver->getNumberOfVariables(),aderdgSolver->getNodesPerCoordinateAxis(),
            minOfNeighbours,maxOfNeighbours);
    if (initialSolutionNeedsLimiting) { // This means that the limiting of the previous iteration was not diffusive enough
      finiteVolumesSolver->updateSolution(cellDescriptionsIndex,finiteVolumesElement,fineGridVertices,fineGridVerticesEnumerator);

      // TODO(Dominic): Project on old solution.
    } else {
      aderdgSolver->updateSolution(cellDescriptionsIndex,finiteVolumesElement,fineGridVertices,fineGridVerticesEnumerator);
    }

    // TODO(Dominic): Project on FV subgrid.
  }
}

void exahype::solvers::ADERDGAPosterioriSubcellLimiter::couple(
    const int cellDescriptionsIndex,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  exahype::solvers::ADERDGSolver* aderdgSolver =
      static_cast<exahype::solvers::ADERDGSolver*>(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]);
  exahype::solvers::FiniteVolumesSolver* finiteVolumesSolver =
        static_cast<exahype::solvers::FiniteVolumesSolver*>(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]);
  int aderdgElement        = aderdgSolver->tryGetElement(cellDescriptionsIndex,_aderdgSolverNumber);
  int finiteVolumesElement = finiteVolumesSolver->tryGetElement(cellDescriptionsIndex,_finiteVolumesSolverNumber);

  if (aderdgElement!=exahype::solvers::Solver::NotFound) {
    auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
        cellDescriptionsIndex)[aderdgElement];
    assertion2(finiteVolumesElement!=exahype::solvers::Solver::NotFound,cellDescriptionsIndex,_aderdgSolverNumber);
    auto& finiteVolumesCellDescription = exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(
        cellDescriptionsIndex)[finiteVolumesElement];

    // TODO(DOMINIC): We already have the old FV solution available.

    double* aderdgSolution  = exahype::solvers::ADERDGSolver::DataHeap::getInstance().
        getData(aderdgCellDescription.getSolution()).data();
    double* minOfNeighbours =  exahype::solvers::ADERDGSolver::DataHeap::getInstance().
        getData(aderdgCellDescription.getSolutionMin()).data();
    double* maxOfNeighbours =  exahype::solvers::ADERDGSolver::DataHeap::getInstance().
        getData(aderdgCellDescription.getSolutionMax()).data();

    bool oldSolutionNeedsLimiting =
        kernels::limiter::generic::c::isTroubledCell(
        aderdgSolution,aderdgSolver->getNumberOfVariables(),aderdgSolver->getNodesPerCoordinateAxis(),
        minOfNeighbours,maxOfNeighbours);
    if (oldSolutionNeedsLimiting) { // This means that the limiting of the previous iteration was not diffusive enough
      double* finiteVolumesSolution = exahype::solvers::FiniteVolumesSolver::DataHeap::getInstance().
              getData(finiteVolumesCellDescription.getSolution()).data();

      finiteVolumesSolver->updateSolution(cellDescriptionsIndex,finiteVolumesElement,fineGridVertices,fineGridVerticesEnumerator);

      // TODO(Dominic): Project on old solution.
    } else {
      aderdgSolver->updateSolution(cellDescriptionsIndex,finiteVolumesElement,fineGridVertices,fineGridVerticesEnumerator);
    }

    // TODO(Dominic): Project on FV subgrid.
  }
}

//void exahype::solvers::ADERDGAPosterioriSubcellLimiter::coupleSolversAfterSolutionUpdate(
//    const int cellDescriptionsIndex,
//    exahype::Vertex* const fineGridVertices,
//    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
//  exahype::solvers::ADERDGSolver* aderdgSolver =
//      static_cast<exahype::solvers::ADERDGSolver*>(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]);
//  int aderdgElement = aderdgSolver->tryGetElement(cellDescriptionsIndex,_aderdgSolverNumber);
//
//  if (aderdgElement!=exahype::solvers::Solver::NotFound) {
//    exahype::solvers::FiniteVolumesSolver* finiteVolumesSolver =
//        static_cast<exahype::solvers::FiniteVolumesSolver*>(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]);
//    int finiteVolumesElement = finiteVolumesSolver->tryGetElement(cellDescriptionsIndex,_finiteVolumesSolverNumber);
//    assertion2(finiteVolumesElement!=exahype::solvers::Solver::NotFound,cellDescriptionsIndex,_finiteVolumesSolverNumber);
//
//    auto& aderdgCellDescription        = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
//        cellDescriptionsIndex)[aderdgElement];
//    auto& finiteVolumesCellDescription = exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(
//        cellDescriptionsIndex)[finiteVolumesElement];
//
//    if (_limiterIsActive) {
//
//      // TODO(JM):
//      // In the second method, you project the fv solution onto the aderdg solution again.
//      // You get all you need from the cell descriptions
//
//      // TODO(Dominic): Project FV solution on DG solution space
//    } else {
//      // TODO(Dominic): Project ADER-DG solution on FV solution space
//    }
//  }
//}
