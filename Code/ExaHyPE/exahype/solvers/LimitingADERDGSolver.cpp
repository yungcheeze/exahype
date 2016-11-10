/*
 * LimitingADERDGSolver.cpp
 *
 *  Created on: 26 Oct 2016
 *      Author: dominic
 */

#include "LimitingADERDGSolver.h"

namespace exahype {
namespace solvers {

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
  double admissibleTimeStepSize = std::numeric_limits<double>::max();

  switch (solverPatch.getLimiterStatus()) {
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Ok: {
      admissibleTimeStepSize =
          _solver->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
      break;
    }
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsNeighbourOfTroubledCell: {
      admissibleTimeStepSize =
                  _solver->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
      int finiteVolumesElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(finiteVolumesElement!=exahype::solvers::Solver::NotFound);
      exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,element);
      limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
      limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
      break;
    }
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Troubled:
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsTroubledCell: {
      admissibleTimeStepSize =
          _limiter->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
      int finiteVolumesElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(finiteVolumesElement!=exahype::solvers::Solver::NotFound);
      exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,element);
      solverPatch.setCorrectorTimeStamp(limiterPatch.getTimeStamp());
      solverPatch.setCorrectorTimeStepSize(limiterPatch.getTimeStepSize());
      solverPatch.setPredictorTimeStamp(limiterPatch.getTimeStamp());
      solverPatch.setPredictorTimeStepSize(limiterPatch.getTimeStepSize());
      break;
    }
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
    double** tempStateSizedArrays,
    double** tempUnknowns,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator)  {
  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);
  int finiteVolumesElement;

  double* solverSolution    = nullptr;
  double* limiterSolution = nullptr;

  // 1. Update the solution in the cells
  switch (solverPatch.getLimiterStatus()) {
  case SolverPatch::LimiterStatus::Ok:
    _solver->updateSolution(
        cellDescriptionsIndex,element,
        tempStateSizedArrays,
        tempUnknowns,
        fineGridVertices,fineGridVerticesEnumerator);
    break;
  case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell: {
      _solver->updateSolution(
          cellDescriptionsIndex,element,
          tempStateSizedArrays,
          tempUnknowns,
          fineGridVertices,fineGridVerticesEnumerator);
      finiteVolumesElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(finiteVolumesElement!=exahype::solvers::Solver::NotFound);
      exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,element);

      solverSolution = DataHeap::getInstance().getData(
          solverPatch.getSolution()).data();
      limiterSolution = DataHeap::getInstance().getData(
          limiterPatch.getSolution()).data();

      // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      kernels::limiter::generic::c::projectOnFVLimiterSpace(
          solverSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),limiterSolution);
      break;
    }
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell: {
      _limiter->updateSolution(
          cellDescriptionsIndex,element,
          tempStateSizedArrays,
          tempUnknowns,
          fineGridVertices,fineGridVerticesEnumerator);
      finiteVolumesElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(finiteVolumesElement!=exahype::solvers::Solver::NotFound);
      exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,element);

      // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      kernels::limiter::generic::c::projectOnDGSpace(limiterSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),solverSolution);
      break;
    }
  }
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
 * Determine the limiter status of a cell description after a solution update.
 * Check which cell is troubled and which cell holds a valid solution.
 * If the limiter subdomain changes, i.e., if a cell changes from holding a
 * valid solution (Ok) to troubled (Troubled) or vice versa,
 * this function returns true.
 *
 * Further reinitialise the merged limiter status on every
 * face with either Ok or Troubled.
 *
 * Do not override the cell-based limiter status of the solver only
 * the merged ones on the faces. We need to hold this value a
 * little longer.
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


void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourLimiterStatus(
    SolverPatch& solverPatch,
    const int faceIndex,
    const SolverPatch::LimiterStatus& neighbourLimiterStatus) const {
  SolverPatch::LimiterStatus limiterStatus = determineLimiterStatus(solverPatch);

  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok:
      switch (neighbourLimiterStatus) {
        case SolverPatch::LimiterStatus::Troubled:
          solverPatch.setMergedLimiterStatus(faceIndex,SolverPatch::LimiterStatus::NeighbourIsTroubledCell);
          break;
        case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
          solverPatch.setMergedLimiterStatus(faceIndex,SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell);
          break;
        default:
          break;
      }
      break;
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      switch (neighbourLimiterStatus) {
        case SolverPatch::LimiterStatus::Troubled:
          solverPatch.setMergedLimiterStatus(faceIndex,SolverPatch::LimiterStatus::NeighbourIsTroubledCell);
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }
}

/**
 * Iterate over the merged limiter statuses per face and
 * determine a unique value.
 */
exahype::solvers::LimitingADERDGSolver::SolverPatch::LimiterStatus
exahype::solvers::LimitingADERDGSolver::determineLimiterStatus(
    SolverPatch& solverPatch) const {
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
  return limiterStatus;
}


void exahype::solvers::LimitingADERDGSolver::reinitialiseSolvers(
    SolverPatch& solverPatch,
    exahype::Cell& fineGridCell,
    double** tempStateSizedArrays,
    double** tempUnknowns,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  SolverPatch::LimiterStatus oldLimiterStatus = solverPatch.getLimiterStatus();
  SolverPatch::LimiterStatus limiterStatus = determineLimiterStatus(solverPatch);

  int limiterElement = _limiter->tryGetElement(
          fineGridCell.getCellDescriptionsIndex(),solverPatch.getSolverNumber());
  LimiterPatch* limiterPatch = nullptr;

  double* solverSolution  = nullptr;
  double* limiterSolution = nullptr;

  switch (limiterStatus) {
  case SolverPatch::LimiterStatus::Ok:
    if (limiterElement!=exahype::solvers::Solver::NotFound) {
      assertion(oldLimiterStatus==SolverPatch::LimiterStatus::Troubled ||
          oldLimiterStatus==SolverPatch::LimiterStatus::NeighbourIsTroubledCell ||
          oldLimiterStatus==SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell);

      limiterPatch =
          &LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[limiterElement];
      limiterPatch->setType(LimiterPatch::Type::Erased);
      _limiter->ensureNoUnnecessaryMemoryIsAllocated(*limiterPatch);
      LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).erase(
          LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).begin()+limiterElement);
    }
    break;
  case SolverPatch::LimiterStatus::Troubled:
  case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
  case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
    // TODO(Dominic): Add to docu: We need to rollback the solution for this case since we need to supply the
    // NeighbourIsTroubledCell neighbours with solution values from the old time step.
    switch (oldLimiterStatus) {
    case SolverPatch::LimiterStatus::Troubled:  // TODO(Dominic): Add to docu: Here we work with a valid old FVM solution
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      _limiter->rollbackSolution(
          fineGridCell.getCellDescriptionsIndex(),limiterElement,
          fineGridVertices,fineGridVerticesEnumerator);
      break;
    case SolverPatch::LimiterStatus::Ok: // TODO(Dominic): Add to docu: Here we work with a valid old ADER-DG solution
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      if (limiterElement==exahype::solvers::Solver::NotFound) {
        assertion(oldLimiterStatus==SolverPatch::LimiterStatus::Ok);

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
      solverSolution = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

      limiterPatch =
          &LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[limiterElement];
      _limiter->ensureNecessaryMemoryIsAllocated(*limiterPatch);
      limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
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

/**
 * Recompute the solution in cells that have been subject to a limiter status change
 * This method is invoked after the solver reinitialisation
 * (see exahype::solvers::LimitingADERDGSolver::reinitialiseSolvers).
 *
 * It evolves the solution of the solver and limiter in the reinitialised cells to the
 * correct time stamp.
 *
 * We perform the following actions based on the
 * new limiter status:
 *
 * |New Status | Action                                                                                                                                      |
 * ----------------------------------------------------------------------------------------------------------------------------------------------------------|
 * |O          | Do nothing. Solver solution has been evolved correctly before.                                                                              |
 * |T/NT       | Evolve FV solver project result onto the ADER-DG space.                                                                                     |
 * |NNT        | Evolve solver and project its solution onto the limiter solution space. (We had to do a rollback beforehand in the reinitialisation phase.) |
 *
 * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
 *
 * We do not overwrite the old limiter status set in this method.
 * We compute the new limiter status based on the merged limiter statuses associated
 * with the faces.
 *
 * <h2>Overlapping status spreading and reinitialisation with solution reconputation</h2>
 * We can recompute the new solution in cells with status Troubled after one iteration
 * since old solution values from direct neighbours are available then.
 *
 * We can recompute the
 *
 * TODO(Dominic)
 * Adapters:
 * LimitingADERDGSolver LimiterStatusSpreading
 * LimitingADERDGSolver Reinitialisation
 * LimitingADERDGSolver Recomputation
 */
void exahype::solvers::LimitingADERDGSolver::recomputeSolution(
        const int cellDescriptionsIndex, const int element,
         double** tempStateSizedArrays,
         double** tempUnknowns,
        exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus limiterStatus = determineLimiterStatus(solverPatch);

  int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  LimiterPatch* limiterPatch = nullptr;

  double* solverSolution  = nullptr;
  double* limiterSolution = nullptr;

  switch (solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Ok:
      assertion(limiterElement==exahype::solvers::Solver::NotFound);
      // do nothing
      break;
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      _solver->updateSolution(cellDescriptionsIndex,limiterElement,
                              tempStateSizedArrays,tempUnknowns,
                              fineGridVertices,fineGridVerticesEnumerator);

      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      limiterPatch = &_limiter->getCellDescription(cellDescriptionsIndex,element);
      solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
      limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();

      // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      kernels::limiter::generic::c::projectOnFVLimiterSpace(
          solverSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),limiterSolution);
      break;
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->updateSolution(cellDescriptionsIndex,limiterElement,
                               tempStateSizedArrays,tempUnknowns,
                               fineGridVertices,fineGridVerticesEnumerator);

      limiterPatch = &_limiter->getCellDescription(cellDescriptionsIndex,element);
      limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
      solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

      // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      kernels::limiter::generic::c::projectOnDGSpace(
          limiterSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),solverSolution);
      break;
    case SolverPatch::LimiterStatus::Troubled:
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->updateSolution(cellDescriptionsIndex,limiterElement,
                               tempStateSizedArrays,tempUnknowns,
                               fineGridVertices,fineGridVerticesEnumerator);

      limiterPatch = &_limiter->getCellDescription(cellDescriptionsIndex,element);
      limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
      solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

      // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      kernels::limiter::generic::c::projectOnDGSpace(
          limiterSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),solverSolution);
      break;
  }
}

void exahype::solvers::LimitingADERDGSolver::updateLimiterStatus(int cellDescriptionsIndex, int element) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  solverPatch.setLimiterStatus(determineLimiterStatus(solverPatch));

  for (int i=0; i<DIMENSIONS_TIMES_TWO; ++i) {
    solverPatch.setMergedLimiterStatus(i,solverPatch.getLimiterStatus());
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

//  if (solverPatch1.getLimiterStatus()==)
}
