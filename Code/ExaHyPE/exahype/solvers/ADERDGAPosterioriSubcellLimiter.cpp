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

    finiteVolumesSolver->synchroniseTimeStepping(cellDescriptionsIndex,finiteVolumesElement);
    finiteVolumesSolver->setInitialConditions(
        cellDescriptionsIndex,finiteVolumesElement,
        fineGridVertices,fineGridVerticesEnumerator);

    double* finiteVolumesSolution = exahype::solvers::FiniteVolumesSolver::DataHeap::getInstance().
            getData(finiteVolumesCellDescription.getSolution()).data();
    double* aderdgSolution  = exahype::solvers::ADERDGSolver::DataHeap::getInstance().
                getData(aderdgCellDescription.getSolution()).data();
    kernels::limiter::generic::c::updateSubcellWithLimiterData(
        finiteVolumesSolution,
        finiteVolumesSolver->getNumberOfVariables(),
        finiteVolumesSolver->getNodesPerCoordinateAxis(),
        aderdgSolver->getNodesPerCoordinateAxis(),
        aderdgSolution);
    aderdgSolver->synchroniseTimeStepping(cellDescriptionsIndex,aderdgElement);
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

    double* aderdgSolution  = exahype::solvers::ADERDGSolver::DataHeap::getInstance().
        getData(aderdgCellDescription.getSolution()).data();
    double* minOfNeighbours =  exahype::solvers::ADERDGSolver::DataHeap::getInstance().
        getData(aderdgCellDescription.getSolutionMin()).data();
    double* maxOfNeighbours =  exahype::solvers::ADERDGSolver::DataHeap::getInstance().
        getData(aderdgCellDescription.getSolutionMax()).data();

    double* finiteVolumesSolution = exahype::solvers::FiniteVolumesSolver::DataHeap::getInstance().
            getData(finiteVolumesCellDescription.getSolution()).data();

    // 1. We first check if the old solution still needs limiting.
    // 1.1 If so, we update the FV solution and project the result back
    // onto the ADER-DG solution space.
    bool aderdgSolutionNeedsLimiting =
        kernels::limiter::generic::c::isTroubledCell(
        aderdgSolution,aderdgSolver->getNumberOfVariables(),aderdgSolver->getNodesPerCoordinateAxis(),
        minOfNeighbours,maxOfNeighbours);
    if (aderdgSolutionNeedsLimiting) {
      // assertion ( finite volumes solver solution is correct ).

      finiteVolumesSolver->updateSolution(cellDescriptionsIndex,finiteVolumesElement,fineGridVertices,fineGridVerticesEnumerator);

      kernels::limiter::generic::c::updateSubcellWithLimiterData(
              finiteVolumesSolution,
              finiteVolumesSolver->getNumberOfVariables(),
              finiteVolumesSolver->getNodesPerCoordinateAxis(),
              aderdgSolver->getNodesPerCoordinateAxis(),
              aderdgSolution);

      // TODO(Dominic): What about the min max in this case? Keep or recompute?
    } else {
      // 2. If not so, we first update the solution of the ADER-DG solver, and then
      // check if we need to perform limiting to the new ADER-DG solution.
      // 2.1 If so, we update the FV solution and project the result back
      // onto the ADER-DG solution space.
      aderdgSolver->updateSolution(cellDescriptionsIndex,finiteVolumesElement,fineGridVertices,fineGridVerticesEnumerator);

      aderdgSolutionNeedsLimiting =
          kernels::limiter::generic::c::isTroubledCell(
              aderdgSolution,aderdgSolver->getNumberOfVariables(),aderdgSolver->getNodesPerCoordinateAxis(),
              minOfNeighbours,maxOfNeighbours);
      if (aderdgSolutionNeedsLimiting) {
        finiteVolumesSolver->updateSolution(cellDescriptionsIndex,finiteVolumesElement,fineGridVertices,fineGridVerticesEnumerator);

        kernels::limiter::generic::c::updateSubcellWithLimiterData(
            finiteVolumesSolution,
            finiteVolumesSolver->getNumberOfVariables(),
            finiteVolumesSolver->getNodesPerCoordinateAxis(),
            aderdgSolver->getNodesPerCoordinateAxis(),
            aderdgSolution);
      } else {
        // 3. If everything is fine with the ADER-DG solution, we project it back
        // onto the FV solution space to prepare limiting of the ADER-DG
        // solution values in this cell and in the neighbouring cells.
//        kernels::limiter::generic::c::getFVMData()
      }
    }
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
