/*
 * LimitingADERDGSolver.cpp
 *
 *  Created on: 26 Oct 2016
 *      Author: dominic
 */

#include "LimitingADERDGSolver.h"

namespace exahype {
namespace solvers {


tarch::multicore::BooleanSemaphore exahype::solvers::LimitingADERDGSolver::_semaphoreForLimiterHeapAccess;

} /* namespace solvers */
} /* namespace exahype */

bool exahype::solvers::LimitingADERDGSolver::isValidCellDescriptionIndex(
    const int cellDescriptionsIndex) const  {
  return _solver->isValidCellDescriptionIndex(cellDescriptionsIndex);
}

int exahype::solvers::LimitingADERDGSolver::tryGetElement(
    const int cellDescriptionsIndex,
    const int solverNumber) const  {
  return _solver->tryGetElement(cellDescriptionsIndex,solverNumber);
}

///////////////////////////////////
// MODIFY CELL DESCRIPTION
///////////////////////////////////
bool exahype::solvers::LimitingADERDGSolver::enterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber)  {
  return _solver->enterCell(
      fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
      coarseGridVertices,coarseGridVerticesEnumerator,coarseGridCell,
      fineGridPositionOfCell,solverNumber);
}

bool exahype::solvers::LimitingADERDGSolver::leaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber)  {
  return _solver->leaveCell(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
          coarseGridVertices,coarseGridVerticesEnumerator,coarseGridCell,
          fineGridPositionOfCell,solverNumber);
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////
double exahype::solvers::LimitingADERDGSolver::startNewTimeStep(
    const int cellDescriptionsIndex,
    const int element,
    double*   tempEigenvalues)  {
  exahype::solvers::ADERDGSolver::CellDescription& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);
  exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch;
  double admissibleTimeStepSize = std::numeric_limits<double>::max();

  switch (solverPatch.getLimiterStatus()) {
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Ok:
      admissibleTimeStepSize =
          _solver->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
      break;
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      admissibleTimeStepSize =
                  _solver->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
      int finiteVolumesElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(finiteVolumesElement!=exahype::solvers::Solver::NotFound);
      limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,element);
      limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
      limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
      break;
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Troubled:
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsTroubledCell:
      admissibleTimeStepSize =
          _limiter->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
      int finiteVolumesElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(finiteVolumesElement!=exahype::solvers::Solver::NotFound);
      limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,element);
      solverPatch.setCorrectorTimeStamp(limiterPatch.getTimeStamp());
      solverPatch.setCorrectorTimeStepSize(limiterPatch.getTimeStepSize());
      solverPatch.setPredictorTimeStamp(limiterPatch.getTimeStamp());
      solverPatch.setPredictorTimeStepSize(limiterPatch.getTimeStepSize());
      break;
  }

  return admissibleTimeStepSize;
}

void exahype::solvers::LimitingADERDGSolver::setInitialConditions(
    const int cellDescriptionsIndex,
    const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  _solver->setInitialConditions(
      cellDescriptionsIndex,element,
      fineGridVertices,fineGridVerticesEnumerator);
}

/**
 * This method assumes the ADERDG solver's cell-local limiter status has
 * already been determined.
 */
void exahype::solvers::LimitingADERDGSolver::updateSolution(
    const int cellDescriptionsIndex,
    const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator)  {
  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);
  exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch;
  int finiteVolumesElement;

  double* solverSolution    = nullptr;
  double* update      = nullptr;
  double* solutionMin = nullptr;
  double* solutionMax = nullptr;

  double* limiterSolution = nullptr;

  // 1. Update the solution in the cells
  switch (solverPatch.getLimiterStatus()) {
  case SolverPatch::LimiterStatus::Ok:
    _solver->updateSolution(
        cellDescriptionsIndex,element,
        fineGridVertices,fineGridVerticesEnumerator);
    break;
  case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
    _solver->updateSolution(
        cellDescriptionsIndex,element,
        fineGridVertices,fineGridVerticesEnumerator);
    finiteVolumesElement =
        _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    assertion(finiteVolumesElement!=exahype::solvers::Solver::NotFound);
    limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,element);

    solverSolution = DataHeap::getInstance().getData(
        solverPatch.getSolution()).data();
    limiterSolution = DataHeap::getInstance().getData(
        limiterPatch.getSolution()).data();

    // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
    kernels::limiter::generic::c::projectOnFVLimiterSpace(
        solverSolution,_solver->getNumberOfVariables(),
        _solver->getNodesPerCoordinateAxis(),limiterSolution);
    break;
  case SolverPatch::LimiterStatus::Troubled:
  case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
    _limiter->updateSolution(
        cellDescriptionsIndex,element,
        fineGridVertices,fineGridVerticesEnumerator);
    finiteVolumesElement =
        _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    assertion(finiteVolumesElement!=exahype::solvers::Solver::NotFound);
    limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,element);

    // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
    kernels::limiter::generic::c::projectOnDGSpace(limiterSolution,_solver->getNumberOfVariables(),
        _solver->getNodesPerCoordinateAxis(),solverSolution);
    break;
  }
}

/**
 * Determine the limiter status of a cell description after a solution update.
 * Check which cell is troubled and which cell holds a valid solution.
 * If the limiter subdomain changes, i.e., if a cell changes from holding a
 * valid solution (Ok) to troubled (Troubled) or vice versa,
 * this function returns true.
 *
 * Further reinitialise the merged limiter status on every
 * face with either Ok or Troubled.
 */
bool exahype::solvers::LimitingADERDGSolver::determineLimiterStatusAfterSolutionUpdate(
    const int cellDescriptionsIndex,
    const int element) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  bool limiterDomainHasChanged=false;

  switch (solverPatch.getLimiterStatus()) {
  case SolverPatch::LimiterStatus::Ok:
  case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
  case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
    if (solutionIsTroubled(solverPatch)) {
      solverPatch.setLimiterStatus(SolverPatch::LimiterStatus::Troubled);
      limiterDomainHasChanged = true;

      for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
        solverPatch.setMergedLimiterStatus(i,SolverPatch::LimiterStatus::Troubled);
      }
    } else {
      for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
        solverPatch.setMergedLimiterStatus(i,SolverPatch::LimiterStatus::Ok);
      }
    }
    break;
  case SolverPatch::LimiterStatus::Troubled:
    if (!solutionIsTroubled(solverPatch)) {
      solverPatch.setLimiterStatus(SolverPatch::LimiterStatus::Ok);
      limiterDomainHasChanged = true;

      for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
        solverPatch.setMergedLimiterStatus(i,SolverPatch::LimiterStatus::Ok);
      }
    } else {
      for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
        solverPatch.setMergedLimiterStatus(i,SolverPatch::LimiterStatus::Troubled);
      }
    }
    break;
  default:
    break;
  }

  return limiterDomainHasChanged;
}

/**
 * Checks if updated solution
 * of the ADER-DG solver is valid
 * or if it contains unphysical oscillations.
 */
bool exahype::solvers::LimitingADERDGSolver::solutionIsTroubled(SolverPatch& solverPatch) {
  double* solution    = nullptr;
  double* update      = nullptr;
  double* solutionMin = nullptr;
  double* solutionMax = nullptr;

  solution = DataHeap::getInstance().getData(
      solverPatch.getSolution()).data();
  update = DataHeap::getInstance().getData(
      solverPatch.getUpdate()).data();
  solutionMin = DataHeap::getInstance().getData(
      solverPatch.getSolutionMin()).data();
  solutionMax = DataHeap::getInstance().getData(
      solverPatch.getSolutionMax()).data();

  return kernels::limiter::generic::c::isTroubledCell(
      solution,_solver->getNumberOfVariables(),
      _solver->getNodesPerCoordinateAxis(),solutionMin,solutionMax);
}

/**
 * Assumes that the status of the cells is either Ok or Troubled
 * at the start of the limiter spreading.
 */
void exahype::solvers::LimitingADERDGSolver::unifyMergedLimiterStatus(
    SolverPatch& solverPatch) {
  // 1. Determine new limiter status.
  SolverPatch::LimiterStatus limiterStatus = SolverPatch::LimiterStatus::Ok;
  for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
    switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok:
      switch (solverPatch.getMergedLimiterStatus(i)) {
      case SolverPatch::LimiterStatus::Troubled:
        limiterStatus = SolverPatch::LimiterStatus::Troubled;
        break;
      case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
        limiterStatus = SolverPatch::LimiterStatus::NeighbourIsTroubledCell;
        break;
      case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
        limiterStatus = SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell;
        break;
      default:
        break;
      }
      break;
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      switch (solverPatch.getMergedLimiterStatus(i)) {
      case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
        limiterStatus = SolverPatch::LimiterStatus::NeighbourIsTroubledCell;
        break;
      default:
        break;
      }
      break;
      default:
        break;
    }
  }

  // 2. Finally, set the merged limiter status on all faces to
  // the determined value.
  // We do not update the final limiter status value since
  // it contains historic information.
  for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
    solverPatch.setMergedLimiterStatus(i,limiterStatus);
  }
}

/**
 * TODO(Dominic):
 *
 * This method is invoked after the limiter status spreading. We need to rollback here
 * if a cell needs to use the FV limiter. After we have finished the rollback in all changed cells,
 * we can recompute them.
 *
 * Adapters:
 * LimitingADERDGSolver LimiterStatusSpreading
 * LimitingADERDGSolver Reinitialisation
 * LimitingADERDGSolver Recomputation
 */
void exahype::solvers::LimitingADERDGSolver::reinitialiseSolvers(
    SolverPatch& solverPatch,
    const exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  SolverPatch::LimiterStatus previousLimiterStatus = solverPatch.getLimiterStatus();
  unifyMergedLimiterStatus(solverPatch);
  SolverPatch::LimiterStatus& limiterStatus = solverPatch.getMergedLimiterStatus(0);

  int limiterElement =
      _limiter->tryGetElement(
          fineGridCell.getCellDescriptionsIndex(),solverPatch.getSolverNumber());
  double* solverSolution  = nullptr;
  double* limiterSolution = nullptr;

  switch (limiterStatus) {
  case SolverPatch::LimiterStatus::Ok:
    if (limiterElement!=exahype::solvers::Solver::NotFound) {
      assertion(previousLimiterStatus==SolverPatch::LimiterStatus::Troubled ||
          previousLimiterStatus==SolverPatch::LimiterStatus::NeighbourIsTroubledCell ||
          previousLimiterStatus==SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell);

      LimiterPatch& limiterPatch =
          LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[limiterElement];
      limiterPatch.setType(LimiterPatch::Type::Erased);
      _limiter->ensureNoUnnecessaryMemoryIsAllocated(limiterPatch);
      LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).erase(
          LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).begin()+limiterElement);
    }
    break;
  case SolverPatch::LimiterStatus::Troubled:
  case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
  case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
    switch (previousLimiterStatus) {
    case SolverPatch::LimiterStatus::Troubled:  // TODO(Dominic): Add to docu: Here we work with a valid old FVM solution
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      _limiter->rollbackSolution(
          fineGridCell.getCellDescriptionsIndex(),limiterElement,
          fineGridVertices,fineGridVerticesEnumerator); // TODO(Dominic): Add update or old solution field such that we can perform a rollback
      break;
    case SolverPatch::LimiterStatus::Ok: // TODO(Dominic): Add to docu: Here we work with a valid old ADER-DG solution
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      if (limiterElement==exahype::solvers::Solver::NotFound) {
        assertion(previousLimiterStatus==SolverPatch::LimiterStatus::Ok);

        // TODO(Dominic): Use solver patch's cell type
        // TODO(Dominic): This is some sort of mesh refinement.
        // In AMR settings, we need to add ancestor and descendant cells
        // after we add a limiter cell description (this will be fun).
        fineGridCell.addNewCellDescription(
            solverPatch.getSolverNumber(),
            LimiterPatch::Type::Cell,
            solverPatch.getLevel(),
            solverPatch.getParentIndex(),
            solverPatch.getSize(),
            solverPatch.getOffset());
      }

      _solver->rollbackSolution(
          fineGridCell.getCellDescriptionsIndex(),limiterElement,
          fineGridVertices,fineGridVerticesEnumerator);
      solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution());

      LimiterPatch& limiterPatch =
          LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[limiterElement];
      _limiter->ensureNecessaryMemoryIsAllocated(limiterPatch);
      limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution());
      kernels::limiter::generic::c::projectOnFVLimiterSpace(
          solverSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),limiterSolution); // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      break;
    }
    break;
  default:
    break;
  }
}

void exahype::solvers::LimitingADERDGSolver::recomputeSolution(
        const int cellDescriptionsIndex,
         const int element,
        exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus previousLimiterStatus = solverPatch.getLimiterStatus();
  solverPatch.setLimiterStatus(solverPatch.getMergedLimiterStatus(0));

  int limiterElement =
      _limiter->tryGetElement(
          cellDescriptionsIndex,solverPatch.getSolverNumber());

  switch (solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Ok:
      assertion(limiterElement==exahype::solvers::Solver::NotFound);
      // do nothing
      break;
    case SolverPatch::LimiterStatus::NewlyNeighbourIsNeighbourOfTroubledCell:
      _solver->updateSolution(cellDescriptionsIndex,limiterElement,
                              fineGridVertices,fineGridVerticesEnumerator);
      break;
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->updateSolution(cellDescriptionsIndex,limiterElement,
                               fineGridVertices,fineGridVerticesEnumerator);
      break;
    case SolverPatch::LimiterStatus::Troubled:
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->updateSolution(cellDescriptionsIndex,limiterElement,
                               fineGridVertices,fineGridVerticesEnumerator);
      break;
  }
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeNeighbours(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    double**                                  tempFaceUnknownsArrays,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  SolverPatch& solverPatch1 = _solver->getCellDescription(cellDescriptionsIndex1,element1);
  SolverPatch& solverPatch2 = _solver->getCellDescription(cellDescriptionsIndex2,element2);

  if (solverPatch1.getLimiterStatus()==)
}
