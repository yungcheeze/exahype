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

tarch::logging::Log exahype::solvers::ADERDGAPosterioriSubcellLimiter::_log("exahype::solvers::ADERDGAPosterioriSubcellLimiter");

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

    // Set the initial conditions of the finite volumes solver.
    finiteVolumesSolver->synchroniseTimeStepping(cellDescriptionsIndex,finiteVolumesElement);
    finiteVolumesSolver->setInitialConditions(
        cellDescriptionsIndex,finiteVolumesElement,
        fineGridVertices,fineGridVerticesEnumerator);

    // Set the initial conditions of the ADER-DG solver by projecting
    // the solution values of the FV solver onto the DG solution space.
    aderdgSolver->synchroniseTimeStepping(cellDescriptionsIndex,aderdgElement);
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

    // This writes the min max values to face 0;
    double* solutionMin = DataHeap::getInstance().getData(aderdgCellDescription.getSolutionMin()).data();
    double* solutionMax = DataHeap::getInstance().getData(aderdgCellDescription.getSolutionMax()).data();
    kernels::limiter::generic::c::findCellLocalLimMinAndMax(
        finiteVolumesSolution,
        aderdgSolver->getNumberOfVariables(),
        aderdgSolver->getNodesPerCoordinateAxis(),
        solutionMin,
        solutionMax);
    // This writes the face 0 min max values to the other faces
    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
      std::copy(solutionMin, solutionMin+aderdgSolver->getNumberOfVariables(),
          solutionMin+i*aderdgSolver->getNumberOfVariables());
      std::copy(solutionMax, solutionMax+aderdgSolver->getNumberOfVariables(),
          solutionMax+i*aderdgSolver->getNumberOfVariables());
    }
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
    aderdgCellDescription.setLimiterStatus(exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Ok);

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
      aderdgCellDescription.setLimiterStatus(exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Troubled);
      logInfo("couple(...)","ADER-DG solution in cell "<<
          aderdgCellDescription.getOffset()<<"x"<<aderdgCellDescription.getSize()<<" still needs limiting.");

      // assertion ( finite volumes solver solution is correct ). todo

      finiteVolumesSolver->updateSolution(cellDescriptionsIndex,finiteVolumesElement,fineGridVertices,fineGridVerticesEnumerator);

      kernels::limiter::generic::c::updateSubcellWithLimiterData(
              finiteVolumesSolution,
              finiteVolumesSolver->getNumberOfVariables(),
              finiteVolumesSolver->getNodesPerCoordinateAxis(),
              aderdgSolver->getNodesPerCoordinateAxis(),
              aderdgSolution);

      // TODO(Dominic): We have to provide ADER-DG boundary data to the neighbours based on the
      // subcells. We need some boundary extrapolation from the troubled cell here.
      // I do now understand why Michael has an extra layer of FV cells around the
      // troubled cells.

      // TODO(Dominic): What about the min max in this case? Keep or recompute?
    } else {
      logInfo("couple(...)","Updating ADER-DG solution in cell "<<
          aderdgCellDescription.getOffset()<<"x"<<aderdgCellDescription.getSize()<<".");

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
        aderdgCellDescription.setLimiterStatus(exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Troubled);
        logInfo("couple(...)","New ADER-DG solution in cell " <<
            aderdgCellDescription.getOffset()<<"x"<<aderdgCellDescription.getSize()<<" needs limiting.");
        finiteVolumesSolver->updateSolution(cellDescriptionsIndex,finiteVolumesElement,fineGridVertices,fineGridVerticesEnumerator);

        kernels::limiter::generic::c::updateSubcellWithLimiterData(
            finiteVolumesSolution,
            finiteVolumesSolver->getNumberOfVariables(),
            finiteVolumesSolver->getNodesPerCoordinateAxis(),
            aderdgSolver->getNodesPerCoordinateAxis(),
            aderdgSolution);

        // TODO(Dominic): We have to provide ADER-DG boundary data to the neighbours based on the
        // subcells. We need some boundary extrapolation from the troubled cell here.
        // I do now understand why Michael has an extra layer of FV cells around the
        // troubled cells.
      } else {
        logInfo("couple(...)","New ADER-DG solution in cell " <<
                    aderdgCellDescription.getOffset()<<"x"<<aderdgCellDescription.getSize()<<" does not need limiting.");
        // 3. If everything is fine with the ADER-DG solution, we project it back
        // onto the FV solution space to prepare limiting of the ADER-DG
        // solution values in this cell and in the neighbouring cells.
        kernels::limiter::generic::c::getFVMData(
            aderdgSolution,
            aderdgSolver->getNumberOfVariables(),
            aderdgSolver->getNodesPerCoordinateAxis(),
            finiteVolumesSolution);

        // This writes the min max values into an array for face 0;
        double* solutionMin = DataHeap::getInstance().getData(aderdgCellDescription.getSolutionMin()).data();
        double* solutionMax = DataHeap::getInstance().getData(aderdgCellDescription.getSolutionMax()).data();
        kernels::limiter::generic::c::findCellLocalMinAndMax(
            aderdgSolution,
            finiteVolumesSolution,
            aderdgSolver->getNumberOfVariables(),
            aderdgSolver->getNodesPerCoordinateAxis(),
            solutionMin,
            solutionMax);
        // This writes the face 0 min max values to arrays for the other faces
        for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
          std::copy(solutionMin, solutionMin+aderdgSolver->getNumberOfVariables(),
              solutionMin+i*aderdgSolver->getNumberOfVariables());
          std::copy(solutionMax, solutionMax+aderdgSolver->getNumberOfVariables(),
              solutionMax+i*aderdgSolver->getNumberOfVariables());
        }
      }
    }
  }
}
