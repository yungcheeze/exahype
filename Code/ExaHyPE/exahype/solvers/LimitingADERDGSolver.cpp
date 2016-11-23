/*
 * LimitingADERDGSolver.cpp
 *
 *  Created on: 26 Oct 2016
 *      Author: dominic
 */

#include "LimitingADERDGSolver.h"

#include "tarch/multicore/Lock.h"
#include "kernels/limiter/generic/Limiter.h"

namespace exahype {
namespace solvers {

} /* namespace solvers */
} /* namespace exahype */

tarch::logging::Log exahype::solvers::LimitingADERDGSolver::_log("exahype::solvers::LimitingADERDGSolver");

bool exahype::solvers::LimitingADERDGSolver::isValidCellDescriptionIndex(
    const int cellDescriptionsIndex) const  {
  return _solver->isValidCellDescriptionIndex(cellDescriptionsIndex);
}

int exahype::solvers::LimitingADERDGSolver::tryGetElement(
    const int cellDescriptionsIndex,
    const int solverNumber) const  {
  return _solver->tryGetElement(cellDescriptionsIndex,solverNumber);
}

int exahype::solvers::LimitingADERDGSolver::tryGetLimiterElement(
    const int cellDescriptionsIndex,
    const int solverNumber) const  {
  return _limiter->tryGetElement(cellDescriptionsIndex,solverNumber);
}

exahype::solvers::Solver::SubcellPosition exahype::solvers::LimitingADERDGSolver::computeSubcellPositionOfCellOrAncestor(
      const int cellDescriptionsIndex,
      const int element) {
  return _solver->computeSubcellPositionOfCellOrAncestor(cellDescriptionsIndex,element);
}

exahype::solvers::LimitingADERDGSolver::LimitingADERDGSolver(
    const std::string& identifier,
    std::unique_ptr<exahype::solvers::ADERDGSolver> solver,
    std::unique_ptr<exahype::solvers::FiniteVolumesSolver> limiter)
    :
    exahype::solvers::Solver(identifier, Solver::Type::LimitingADERDG, solver->getNumberOfVariables(),
        solver->getNumberOfParameters(), solver->getNodesPerCoordinateAxis(), solver->getMaximumMeshSize(),
          solver->getTimeStepping()),
          _solver(std::move(solver)),
          _limiter(std::move(limiter)),
          _limiterDomainHasChanged(false)
{
  assertion(_solver->getNumberOfParameters() == 0);
}

double exahype::solvers::LimitingADERDGSolver::getMinTimeStamp() const {
  return _solver->getMinTimeStamp();
}

double exahype::solvers::LimitingADERDGSolver::getMinTimeStepSize() const {
  return _solver->getMinTimeStepSize();
}

double exahype::solvers::LimitingADERDGSolver::getMinNextTimeStepSize() const {
  return _solver->getMinNextTimeStepSize();
}

void exahype::solvers::LimitingADERDGSolver::updateMinNextTimeStepSize(double value) {
  _solver->updateMinNextTimeStepSize(value);
  _limiter->updateMinNextTimeStepSize(value);
}

void exahype::solvers::LimitingADERDGSolver::initInitialTimeStamp(double value) {
  _solver->initInitialTimeStamp(value);
  _limiter->initInitialTimeStamp(value);
}

void exahype::solvers::LimitingADERDGSolver::synchroniseTimeStepping(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->synchroniseTimeStepping(cellDescriptionsIndex,element);

  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);

  const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    assertion(solverPatch.getLimiterStatus()==SolverPatch::LimiterStatus::Troubled
              || solverPatch.getLimiterStatus()==SolverPatch::LimiterStatus::NeighbourIsTroubledCell
              || solverPatch.getLimiterStatus()==SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell);
    _limiter->synchroniseTimeStepping(cellDescriptionsIndex,limiterElement);
  }
}

void exahype::solvers::LimitingADERDGSolver::startNewTimeStep() {
  double minNextTimeStepSize =
            std::min( _solver->getMinNextPredictorTimeStepSize(), _limiter->getMinNextTimeStepSize() );

  switch (_timeStepping) {
    case TimeStepping::Global:
      _solver->updateMinNextPredictorTimeStepSize(minNextTimeStepSize);
      _limiter->updateMinNextTimeStepSize(minNextTimeStepSize);

      _solver->startNewTimeStep();
      _limiter->startNewTimeStep();
      break;
    case TimeStepping::GlobalFixed:
      _solver->startNewTimeStep();
      _limiter->startNewTimeStep();
      break;
  } // TODO(Dominic): Switch-case probably not necessary
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep() {
  _solver->rollbackToPreviousTimeStep();
  _limiter->rollbackToPreviousTimeStep();
}

void exahype::solvers::LimitingADERDGSolver::reconstructStandardTimeSteppingDataAfterRollback() {
  _solver->reconstructStandardTimeSteppingDataAfterRollback();

  _limiter->setMinTimeStamp(_solver->getMinCorrectorTimeStamp());
  _limiter->setMinTimeStepSize(_solver->getMinCorrectorTimeStepSize());
  _limiter->setMinNextTimeStepSize(std::numeric_limits<double>::max());
  _limiter->setPreviousMinTimeStepSize(std::numeric_limits<double>::max());
}

void exahype::solvers::LimitingADERDGSolver::reinitialiseTimeStepData() {
  _solver->reinitialiseTimeStepData();
  _limiter->reinitialiseTimeStepData();
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
      int limiterElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

      limiterPatch.setPreviousTimeStepSize(solverPatch.getPreviousCorrectorTimeStepSize());
      limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
      limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
      break;
    }
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Troubled:
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsTroubledCell: {
      admissibleTimeStepSize =
          _limiter->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
      int limiterElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      LimiterPatch& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

      solverPatch.setPreviousCorrectorTimeStepSize(limiterPatch.getPreviousTimeStepSize());
      solverPatch.setCorrectorTimeStamp(limiterPatch.getTimeStamp());
      solverPatch.setCorrectorTimeStepSize(limiterPatch.getTimeStepSize());

      solverPatch.setPredictorTimeStamp(limiterPatch.getTimeStamp()+limiterPatch.getTimeStepSize());
      solverPatch.setPredictorTimeStepSize(limiterPatch.getTimeStepSize()); // TODO(Dominic): Reassess this for Standard time stepping.
      break;
    }
  }

  return admissibleTimeStepSize;
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep(
    const int cellDescriptionsIndex,
    const int solverElement) {
  _solver->rollbackToPreviousTimeStep(cellDescriptionsIndex,solverElement);

  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,solverElement);
  if (solverPatch.getLimiterStatus()==SolverPatch::Troubled
      || solverPatch.getLimiterStatus()==SolverPatch::NeighbourIsTroubledCell
      || solverPatch.getLimiterStatus()==SolverPatch::NeighbourIsNeighbourOfTroubledCell) {
    int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    assertion(limiterElement!=exahype::solvers::Solver::NotFound);
    _limiter->rollbackToPreviousTimeStep(cellDescriptionsIndex,limiterElement);
  }
}

void exahype::solvers::LimitingADERDGSolver::reconstructStandardTimeSteppingDataAfterRollback(
    const int cellDescriptionsIndex,
    const int solverElement) const {
  _solver->reconstructStandardTimeSteppingDataAfterRollback(cellDescriptionsIndex,solverElement);

  SolverPatch& solverPatch =
        _solver->getCellDescription(cellDescriptionsIndex,solverElement);
  if (solverPatch.getLimiterStatus()==SolverPatch::Troubled
      || solverPatch.getLimiterStatus()==SolverPatch::NeighbourIsTroubledCell
      || solverPatch.getLimiterStatus()==SolverPatch::NeighbourIsNeighbourOfTroubledCell) {
    int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    assertion(limiterElement!=exahype::solvers::Solver::NotFound);
    LimiterPatch& limiterPatch =
        _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

    limiterPatch.setPreviousTimeStepSize(solverPatch.getPreviousCorrectorTimeStepSize());
    limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
    limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
  }
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

void exahype::solvers::LimitingADERDGSolver::initialiseLimiter(
    const int cellDescriptionsIndex,
    const int element,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus limiterStatus = determineLimiterStatus(solverPatch);

  int limiterElement         = exahype::solvers::Solver::NotFound;
  LimiterPatch* limiterPatch = nullptr;

  double* solverSolution  = nullptr;
  double* limiterSolution = nullptr;

  switch (limiterStatus) {
  case SolverPatch::LimiterStatus::Ok:
    // do nothing
    break;
  case SolverPatch::LimiterStatus::Troubled:
  case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
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

    limiterElement = tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    limiterPatch =
        &LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[limiterElement];
    _limiter->ensureNecessaryMemoryIsAllocated(*limiterPatch);

    limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
    _limiter->setInitialConditions(cellDescriptionsIndex,limiterElement,fineGridVertices,fineGridVerticesEnumerator);

    // Finite Volumes -> ADER-DG
    solverSolution = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
    kernels::limiter::generic::c::projectOnDGSpace(
       limiterSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),
          _limiter->getGhostLayerWidth(),
          solverSolution); // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
    break;
  case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
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

    limiterElement = tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    limiterPatch =
        &LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[limiterElement];
    _limiter->ensureNecessaryMemoryIsAllocated(*limiterPatch);
    limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();

    // ADER-DG -> Finite Volumes
    solverSolution = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
    kernels::limiter::generic::c::projectOnFVLimiterSpace(
        solverSolution,_solver->getNumberOfVariables(),
        _solver->getNodesPerCoordinateAxis(),
        _limiter->getGhostLayerWidth(),
        limiterSolution); // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
    break;
  }
}

/**
 * This method assumes the ADERDG solver's cell-local limiter status has
 * already been determined.
 */
void exahype::solvers::LimitingADERDGSolver::updateSolution(
    const int cellDescriptionsIndex,
    const int element,
    double** tempStateSizedVectors,
    double** tempUnknowns,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator)  {
  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());

  double* solverSolution  = nullptr;
  double* limiterSolution = nullptr;

  // 1. Update the solution in the cells
  switch (solverPatch.getLimiterStatus()) {
  case SolverPatch::LimiterStatus::Ok:
    _solver->updateSolution(
        cellDescriptionsIndex,element,
        tempStateSizedVectors,
        tempUnknowns,
        fineGridVertices,fineGridVerticesEnumerator);
    break;
  case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell: {
    _solver->updateSolution(
        cellDescriptionsIndex,element,
        tempStateSizedVectors,
        tempUnknowns,
        fineGridVertices,fineGridVerticesEnumerator);
    assertion(limiterElement!=exahype::solvers::Solver::NotFound);
    exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
        _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

    solverSolution = DataHeap::getInstance().getData(
        solverPatch.getSolution()).data();
    limiterSolution = DataHeap::getInstance().getData(
        limiterPatch.getSolution()).data();

    // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
    kernels::limiter::generic::c::projectOnFVLimiterSpace(
        solverSolution,_solver->getNumberOfVariables(),
        _solver->getNodesPerCoordinateAxis(),
        _limiter->getGhostLayerWidth(),
        limiterSolution);
    } break;
  case SolverPatch::LimiterStatus::Troubled:
  case SolverPatch::LimiterStatus::NeighbourIsTroubledCell: {
    _limiter->updateSolution(
        cellDescriptionsIndex,element,
        tempStateSizedVectors,
        tempUnknowns,
        fineGridVertices,fineGridVerticesEnumerator);
    assertion(limiterElement!=exahype::solvers::Solver::NotFound);
    exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
        _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

    solverSolution = DataHeap::getInstance().getData(
        solverPatch.getSolution()).data();
    limiterSolution = DataHeap::getInstance().getData(
        limiterPatch.getSolution()).data();

    // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
    kernels::limiter::generic::c::projectOnDGSpace(limiterSolution,_solver->getNumberOfVariables(),
        _solver->getNodesPerCoordinateAxis(),
        _limiter->getGhostLayerWidth(),
        solverSolution);
    } break;
  }
}

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


bool exahype::solvers::LimitingADERDGSolver::solutionIsTroubled(SolverPatch& solverPatch) {
  double* solution = DataHeap::getInstance().getData(
      solverPatch.getSolution()).data();
  double* solutionMin = DataHeap::getInstance().getData(
      solverPatch.getSolutionMin()).data();
  double* solutionMax = DataHeap::getInstance().getData(
      solverPatch.getSolutionMax()).data();

  return kernels::limiter::generic::c::isTroubledCell(
      solution,_solver->getNumberOfVariables(),
      _solver->getNodesPerCoordinateAxis(),solutionMin,solutionMax);
}

void exahype::solvers::LimitingADERDGSolver::determineMinAndMax(
    const int cellDescriptionsIndex,
    const int element) {
  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());

  switch(solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      determineSolverMinAndMax(solverPatch);
      break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());
      exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      determineLimiterMinAndMax(solverPatch,limiterPatch);
      break;
  }
}

void exahype::solvers::LimitingADERDGSolver::determineSolverMinAndMax(SolverPatch& solverPatch) {
  double* solution = DataHeap::getInstance().getData(
      solverPatch.getSolution()).data();
  double* solutionMin = DataHeap::getInstance().getData(
      solverPatch.getSolutionMin()).data();
  double* solutionMax = DataHeap::getInstance().getData(
      solverPatch.getSolutionMax()).data();

  // Write the result to the face with index "0"
  kernels::limiter::generic::c::findCellLocalMinAndMax(
      solution,_numberOfVariables,_nodesPerCoordinateAxis,solutionMin,solutionMax);

  // Copy the result on the other faces as well
  for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
    std::copy(
        solutionMin,
        solutionMin+_numberOfVariables, // past-the-end element
        solutionMin+i*_numberOfVariables);
    std::copy(
        solutionMax,
        solutionMax+_numberOfVariables, // past-the-end element
        solutionMax+i*_numberOfVariables);
  }
}

void exahype::solvers::LimitingADERDGSolver::determineLimiterMinAndMax(SolverPatch& solverPatch,LimiterPatch& limiterPatch) {
  double* limiterSolution = DataHeap::getInstance().getData(
      limiterPatch.getSolution()).data();
  double* solutionMin = DataHeap::getInstance().getData(
      solverPatch.getSolutionMin()).data();
  double* solutionMax = DataHeap::getInstance().getData(
      solverPatch.getSolutionMax()).data();

  // Write the result to the face with index "0"
  kernels::limiter::generic::c::findCellLocalLimiterMinAndMax(
      limiterSolution,_numberOfVariables,_nodesPerCoordinateAxis,_limiter->getGhostLayerWidth(),solutionMin,solutionMax);

  // Copy the result on the other faces as well
  for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
    std::copy(
        solutionMin,
        solutionMin+_numberOfVariables, // past-the-end element
        solutionMin+i*_numberOfVariables);
    std::copy(
        solutionMax,
        solutionMax+_numberOfVariables, // past-the-end element
        solutionMax+i*_numberOfVariables);
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
 *
 * Do not override the cell-based limiter status of the solver only
 * the merged ones on the faces. We need to hold this value a
 * little longer.
 */
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

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourLimiterStatus(
    SolverPatch& solverPatch,
    const int faceIndex,
    const SolverPatch::LimiterStatus& neighbourLimiterStatus,
    const SolverPatch::LimiterStatus& neighbourOfNeighbourLimiterStatus) const {
  SolverPatch::LimiterStatus limiterStatus = determineLimiterStatus(solverPatch);

  mergeWithNeighbourLimiterStatus(solverPatch,faceIndex,neighbourLimiterStatus);

  if (limiterStatus==SolverPatch::LimiterStatus::Ok) {
    if (neighbourLimiterStatus==SolverPatch::LimiterStatus::Ok
        && neighbourOfNeighbourLimiterStatus==SolverPatch::LimiterStatus::Troubled) {
      solverPatch.setMergedLimiterStatus(faceIndex,SolverPatch::LimiterStatus::NeighbourIsTroubledCell);
    }
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
    const int cellDescriptionsIndex,
    const int element,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  SolverPatch::LimiterStatus previousLimiterStatus = solverPatch.getLimiterStatus();
  SolverPatch::LimiterStatus limiterStatus         = determineLimiterStatus(solverPatch);

  int limiterElement = _limiter->tryGetElement(
          fineGridCell.getCellDescriptionsIndex(),solverPatch.getSolverNumber());
  LimiterPatch* limiterPatch = nullptr;

  double* solverSolution  = nullptr;
  double* limiterSolution = nullptr;

  switch (limiterStatus) {
  case SolverPatch::LimiterStatus::Ok:
    if (limiterElement!=exahype::solvers::Solver::NotFound) {
      assertion(previousLimiterStatus==SolverPatch::LimiterStatus::Troubled ||
          previousLimiterStatus==SolverPatch::LimiterStatus::NeighbourIsTroubledCell ||
          previousLimiterStatus==SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell);

      limiterPatch = &_limiter->getCellDescription(fineGridCell.getCellDescriptionsIndex(),limiterElement);
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
    switch (previousLimiterStatus) {
    case SolverPatch::LimiterStatus::Troubled:  // TODO(Dominic): Add to docu: Here we work with a valid old FVM solution
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->rollbackSolution(
          fineGridCell.getCellDescriptionsIndex(),limiterElement,
          fineGridVertices,fineGridVerticesEnumerator);
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

        limiterElement = _limiter->tryGetElement(
            fineGridCell.getCellDescriptionsIndex(),solverPatch.getSolverNumber());
      }
      _solver->rollbackSolution(
          fineGridCell.getCellDescriptionsIndex(),element,
          fineGridVertices,fineGridVerticesEnumerator);
      solverSolution = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

      limiterPatch =
          &LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[limiterElement];
      _limiter->ensureNecessaryMemoryIsAllocated(*limiterPatch);
      limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
      kernels::limiter::generic::c::projectOnFVLimiterSpace(
          solverSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),
          _limiter->getGhostLayerWidth(),
          limiterSolution); // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.

      // Copy time step data from the solver patch
      limiterPatch->setPreviousTimeStepSize(solverPatch.getPreviousCorrectorTimeStepSize());
      limiterPatch->setTimeStamp(solverPatch.getCorrectorTimeStamp());
      limiterPatch->setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
      break;
    }
    break;
  default:
    break;
  }
}

void exahype::solvers::LimitingADERDGSolver::recomputeSolution(
        const int cellDescriptionsIndex, const int element,
         double** tempStateSizedArrays,
         double** tempUnknowns,
        exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus previousLimiterStatus = solverPatch.getLimiterStatus();
  SolverPatch::LimiterStatus limiterStatus         = determineLimiterStatus(solverPatch);

  int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  LimiterPatch* limiterPatch = nullptr;

  double* solverSolution  = nullptr;
  double* limiterSolution = nullptr;

  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok:
      assertion(limiterElement==exahype::solvers::Solver::NotFound);
      // do nothing
      break;
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell: // TODO(Dominic): Add to docu: Here, were just went back one time step to supply the NT neighbours with old limiter unknowns.
      switch (previousLimiterStatus) {
        case SolverPatch::LimiterStatus::Ok:
        case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
          _solver->updateSolution(cellDescriptionsIndex,limiterElement,
                                  tempStateSizedArrays,tempUnknowns,
                                  fineGridVertices,fineGridVerticesEnumerator);

          assertion(limiterElement!=exahype::solvers::Solver::NotFound);
          limiterPatch = &_limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
          solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
          limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();

          // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
          kernels::limiter::generic::c::projectOnFVLimiterSpace(
              solverSolution,_solver->getNumberOfVariables(),
              _solver->getNodesPerCoordinateAxis(),
              _limiter->getGhostLayerWidth(),
              limiterSolution);
          break;
        case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
        case SolverPatch::LimiterStatus::Troubled:
          _limiter->swapSolutionAndPreviousSolution(cellDescriptionsIndex,limiterElement);

          limiterPatch = &_limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
          limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
          solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

          // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
          kernels::limiter::generic::c::projectOnDGSpace(
              limiterSolution,_solver->getNumberOfVariables(),
              _solver->getNodesPerCoordinateAxis(),
              _limiter->getGhostLayerWidth(),
              solverSolution);
          break;
      }
      break;
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);

      limiterPatch = &_limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      assertionEquals1(limiterPatch->getPreviousTimeStepSize(),solverPatch.getPreviousCorrectorTimeStepSize(),cellDescriptionsIndex);
      assertionEquals1(limiterPatch->getTimeStamp(),solverPatch.getCorrectorTimeStamp(),cellDescriptionsIndex);
      assertionEquals1(limiterPatch->getTimeStepSize(),solverPatch.getCorrectorTimeStepSize(),cellDescriptionsIndex);

      // 1. Evolve solution to desired  time step again
      _limiter->updateSolution(cellDescriptionsIndex,limiterElement,
                               tempStateSizedArrays,tempUnknowns,
                               fineGridVertices,fineGridVerticesEnumerator);

      // 2. Project FV solution on ADER-DG space
      limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
      solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

      // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      kernels::limiter::generic::c::projectOnDGSpace(
          limiterSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),
          _limiter->getGhostLayerWidth(),
          solverSolution);
      break;
    case SolverPatch::LimiterStatus::Troubled:
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->updateSolution(cellDescriptionsIndex,limiterElement,
                               tempStateSizedArrays,tempUnknowns,
                               fineGridVertices,fineGridVerticesEnumerator);

      limiterPatch    = &_limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
      solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

      // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      kernels::limiter::generic::c::projectOnDGSpace(
          limiterSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),
          _limiter->getGhostLayerWidth(),
          solverSolution);
      break;
  }
}

void exahype::solvers::LimitingADERDGSolver::prepareNextNeighbourMerging(
    const int cellDescriptionsIndex,const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  _solver->prepareNextNeighbourMerging(
      cellDescriptionsIndex,element,
      fineGridVertices,fineGridVerticesEnumerator);

  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->prepareNextNeighbourMerging(
        cellDescriptionsIndex,limiterElement,
        fineGridVertices,fineGridVerticesEnumerator);
  }
}

void exahype::solvers::LimitingADERDGSolver::updateMergedLimiterStatus(
    const int cellDescriptionsIndex,const int element) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus mergedLimiterStatus = determineLimiterStatus(solverPatch);

  for (int i=0; i<DIMENSIONS_TIMES_TWO; ++i) {
    solverPatch.setMergedLimiterStatus(i,mergedLimiterStatus);
  }
}

void exahype::solvers::LimitingADERDGSolver::updateLimiterStatus(
    const int cellDescriptionsIndex,const int element) const  {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  solverPatch.setLimiterStatus(determineLimiterStatus(solverPatch));
}

void exahype::solvers::LimitingADERDGSolver::preProcess(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->preProcess(cellDescriptionsIndex,element);
  _limiter->preProcess(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::postProcess(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->postProcess(cellDescriptionsIndex,element);
  _limiter->postProcess(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::prolongateDataAndPrepareDataRestriction(
    const int cellDescriptionsIndex,
    const int element) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement = exahype::solvers::Solver::NotFound;

  switch(solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Ok:
      _solver->prolongateDataAndPrepareDataRestriction(cellDescriptionsIndex,element);
      break;
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      _solver->prolongateDataAndPrepareDataRestriction(cellDescriptionsIndex,element);

      limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      _limiter->prolongateDataAndPrepareDataRestriction(cellDescriptionsIndex,limiterElement);
      break;
    case SolverPatch::LimiterStatus::Troubled:
      limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      _limiter->prolongateDataAndPrepareDataRestriction(cellDescriptionsIndex,limiterElement);
      break;
  }
}

void exahype::solvers::LimitingADERDGSolver::restrictData(
    const int cellDescriptionsIndex,
    const int element,
    const int parentCellDescriptionsIndex,
    const int parentElement,
    const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement       = exahype::solvers::Solver::NotFound;
  int parentLimiterElement = exahype::solvers::Solver::NotFound;

  switch(solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Ok:
      _solver->restrictData(cellDescriptionsIndex,element,parentCellDescriptionsIndex,parentElement,subcellIndex);
      break;
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      _solver->restrictData(cellDescriptionsIndex,element,parentCellDescriptionsIndex,parentElement,subcellIndex);

      limiterElement       = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      parentLimiterElement = tryGetLimiterElement(parentCellDescriptionsIndex,solverPatch.getSolverNumber());
      _limiter->restrictData(cellDescriptionsIndex,limiterElement,parentCellDescriptionsIndex,parentLimiterElement,subcellIndex);
      break;
    case SolverPatch::LimiterStatus::Troubled:
      limiterElement       = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      parentLimiterElement = tryGetLimiterElement(parentCellDescriptionsIndex,solverPatch.getSolverNumber());
      _limiter->restrictData(cellDescriptionsIndex,limiterElement,parentCellDescriptionsIndex,parentLimiterElement,subcellIndex);
      break;
  }
}


///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeLimiterStatusOfNeighbours(
      const int                                 SolverPatchsIndex1,
      const int                                 element1,
      const int                                 SolverPatchsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  const int normalDirection = tarch::la::equalsReturnIndex(pos1, pos2);
  assertion(normalDirection >= 0 && normalDirection < DIMENSIONS);
  const int faceIndex1 = 2 * normalDirection +
      (pos2(normalDirection) > pos1(normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
  const int faceIndex2 = 2 * normalDirection +
      (pos1(normalDirection) > pos2(normalDirection) ? 1 : 0);   // !!! Be aware of the ">" !!!

  int SolverPatchsIndexLeft  = SolverPatchsIndex1;
  int elementLeft            = element1;
  int faceIndexLeft          = faceIndex1;

  int SolverPatchsIndexRight = SolverPatchsIndex2;
  int elementRight           = element2;
  int faceIndexRight         = faceIndex2;

  if (pos1(normalDirection) > pos2(normalDirection)) {
    SolverPatchsIndexLeft  = SolverPatchsIndex2;
    elementLeft            = element2;
    faceIndexLeft          = faceIndex2;

    SolverPatchsIndexRight = SolverPatchsIndex1;
    elementRight           = element1;
    faceIndexRight         = faceIndex1;
  }

  SolverPatch& solverPatchLeft  = _solver->getCellDescription(SolverPatchsIndexLeft,elementLeft);
  SolverPatch& solverPatchRight = _solver->getCellDescription(SolverPatchsIndexRight,elementRight);

  // We need to copy the limiter status since the routines below modify
  // the limiter status on the cell descriptions.
  const SolverPatch::LimiterStatus& limiterStatusLeft  = solverPatchLeft.getMergedLimiterStatus(faceIndexLeft);
  const SolverPatch::LimiterStatus& limiterStatusRight = solverPatchRight.getMergedLimiterStatus(faceIndexRight);
  mergeWithNeighbourLimiterStatus(solverPatchLeft,faceIndexLeft,limiterStatusRight);
  mergeWithNeighbourLimiterStatus(solverPatchRight,faceIndexRight,limiterStatusLeft);
}

void exahype::solvers::LimitingADERDGSolver::mergeLimiterStatusOfNeighboursOfNeighbours(
      const int                                 SolverPatchsIndex1,
      const int                                 element1,
      const int                                 SolverPatchsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  const int normalDirection = tarch::la::equalsReturnIndex(pos1, pos2);
  assertion(normalDirection >= 0 && normalDirection < DIMENSIONS);
  const int faceIndex1 = 2 * normalDirection +
      (pos2(normalDirection) > pos1(normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
  const int faceIndex2 = 2 * normalDirection +
      (pos1(normalDirection) > pos2(normalDirection) ? 1 : 0);   // !!! Be aware of the ">" !!!

  int SolverPatchsIndexLeft  = SolverPatchsIndex1;
  int elementLeft            = element1;
  int faceIndexLeft          = faceIndex1;

  int SolverPatchsIndexRight = SolverPatchsIndex2;
  int elementRight           = element2;
  int faceIndexRight         = faceIndex2;

  if (pos1(normalDirection) > pos2(normalDirection)) {
    SolverPatchsIndexLeft  = SolverPatchsIndex2;
    elementLeft            = element2;
    faceIndexLeft          = faceIndex2;

    SolverPatchsIndexRight = SolverPatchsIndex1;
    elementRight           = element1;
    faceIndexRight         = faceIndex1;
  }

  SolverPatch& solverPatchLeft  = _solver->getCellDescription(SolverPatchsIndexLeft,elementLeft);
  SolverPatch& solverPatchRight = _solver->getCellDescription(SolverPatchsIndexRight,elementRight);

  // We need to copy the limiter status since the routines below modify
  // the limiter status on the cell descriptions.
  const SolverPatch::LimiterStatus& limiterStatusLeft       = solverPatchLeft.getMergedLimiterStatus(faceIndexLeft);
  const SolverPatch::LimiterStatus& limiterStatusLeftLeft   = solverPatchLeft.getMergedLimiterStatus(faceIndexRight);
  const SolverPatch::LimiterStatus& limiterStatusRight      = solverPatchRight.getMergedLimiterStatus(faceIndexRight);
  const SolverPatch::LimiterStatus& limiterStatusRightRight = solverPatchRight.getMergedLimiterStatus(faceIndexLeft);
  mergeWithNeighbourLimiterStatus(solverPatchLeft,faceIndexLeft,limiterStatusRight, limiterStatusRightRight);
  mergeWithNeighbourLimiterStatus(solverPatchRight,faceIndexRight,limiterStatusLeft, limiterStatusLeftLeft);
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch& pLeft,
  SolverPatch& pRight,
  const int faceIndexLeft,
  const int faceIndexRight
) const {
  if (pLeft.getType()==SolverPatch::Cell ||
      pRight.getType()==SolverPatch::Cell) {
    assertion( pLeft.getSolverNumber() == pRight.getSolverNumber() );
    const int numberOfVariables = getNumberOfVariables();
    double* minLeft  = DataHeap::getInstance().getData( pLeft.getSolutionMin()  ).data()  + faceIndexLeft  * numberOfVariables;
    double* minRight = DataHeap::getInstance().getData( pRight.getSolutionMin()  ).data() + faceIndexRight * numberOfVariables;
    double* maxLeft  = DataHeap::getInstance().getData( pLeft.getSolutionMax()  ).data()  + faceIndexLeft  * numberOfVariables;
    double* maxRight = DataHeap::getInstance().getData( pRight.getSolutionMax()  ).data() + faceIndexRight * numberOfVariables;

    for (int i=0; i<numberOfVariables; i++) {
      const double min = std::min(
          *(minLeft+i),
          *(minRight+i)
      );
      const double max = std::max(
          *(maxLeft+i),
          *(maxRight+i)
      );

      *(minLeft+i)  = min;
      *(minRight+i) = min;

      *(maxLeft+i)  = max;
      *(maxRight+i) = max;
    }
  } // else do nothing
}

void exahype::solvers::LimitingADERDGSolver::mergeNeighbours(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  // 1. Communicate the data between the solvers that is necessary for the solves
  SolverPatch& solverPatch1 = _solver->getCellDescription(cellDescriptionsIndex1,element1);
  SolverPatch& solverPatch2 = _solver->getCellDescription(cellDescriptionsIndex2,element2);
  mergeNeighboursBasedOnLimiterStatus(
      cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,
      solverPatch1.getLimiterStatus(),solverPatch2.getLimiterStatus(),
      pos1,pos2,
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);

  // 2. Merge the min and max of both cell description's solver's
  // solution value.
  // 2.1 Determine "left" and "right" solverPatch
  const int normalDirection = tarch::la::equalsReturnIndex(pos1, pos2);
  assertion(normalDirection >= 0 && normalDirection < DIMENSIONS);
  const int faceIndex1 = 2 * normalDirection +
      (pos2(normalDirection) > pos1(normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
  const int faceIndex2 = 2 * normalDirection +
      (pos1(normalDirection) > pos2(normalDirection) ? 1 : 0);   // !!! Be aware of the ">" !!!

  int cellDescriptionsIndexLeft  = cellDescriptionsIndex1;
  int elementLeft                = element1;
  int faceIndexLeft              = faceIndex1;

  int cellDescriptionsIndexRight = cellDescriptionsIndex2;
  int elementRight               = element2;
  int faceIndexRight             = faceIndex2;

  if (pos1(normalDirection) > pos2(normalDirection)) {
    cellDescriptionsIndexLeft  = cellDescriptionsIndex2;
    elementLeft                = element2;
    faceIndexLeft              = faceIndex2;

    cellDescriptionsIndexRight = cellDescriptionsIndex1;
    elementRight               = element1;
    faceIndexRight             = faceIndex1;
  }

  SolverPatch& solverPatchLeft  = _solver->getCellDescription(cellDescriptionsIndexLeft,elementLeft);
  SolverPatch& solverPatchRight = _solver->getCellDescription(cellDescriptionsIndexRight,elementRight);

  // 2.2. Merge min/max of both solver patches
  mergeSolutionMinMaxOnFace(solverPatchLeft,solverPatchRight,faceIndexLeft,faceIndexRight);
}

void exahype::solvers::LimitingADERDGSolver::mergeNeighboursBasedOnLimiterStatus(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const SolverPatch::LimiterStatus&         limiterStatus1,
    const SolverPatch::LimiterStatus&         limiterStatus2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  // 2. Merge the boundary data
  // Merge solver solution or limiter solution values in
  // non-overlapping parts of solver and limiter:
  SolverPatch& solverPatch1 = _solver->getCellDescription(cellDescriptionsIndex1,element1);
  SolverPatch& solverPatch2 = _solver->getCellDescription(cellDescriptionsIndex2,element2);
  int limiterElement1 = tryGetLimiterElement(cellDescriptionsIndex1,solverPatch1.getSolverNumber());
  int limiterElement2 = tryGetLimiterElement(cellDescriptionsIndex2,solverPatch2.getSolverNumber());

  switch (limiterStatus1) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      switch (limiterStatus2) {
        case SolverPatch::LimiterStatus::Ok:
        case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
          _solver->mergeNeighbours(cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
                                   tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          break;
        default:
          break;
      }
      break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      switch (limiterStatus2) {
        case SolverPatch::LimiterStatus::Troubled:
        case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
          _limiter->mergeNeighbours(cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
                                    tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          break;
        default:
          break;
      }
      break;
   default:
     break;
  }
  // Merge limiter solution values in overlapping part
  // of solver and limiter:
  switch (limiterStatus1) {
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      switch (limiterStatus2) {
        case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
        case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
          _limiter->mergeNeighbours(cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
                                    tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMax(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  // 1.1. Determine "left" and "right" solverPatch
  const int normalDirection = tarch::la::equalsReturnIndex(pos1, pos2);
  assertion(normalDirection >= 0 && normalDirection < DIMENSIONS);
  const int faceIndex1 = 2 * normalDirection +
      (pos2(normalDirection) > pos1(normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
  const int faceIndex2 = 2 * normalDirection +
      (pos1(normalDirection) > pos2(normalDirection) ? 1 : 0);   // !!! Be aware of the ">" !!!

  int cellDescriptionsIndexLeft  = cellDescriptionsIndex1;
  int elementLeft                = element1;
  int faceIndexLeft              = faceIndex1;

  int cellDescriptionsIndexRight = cellDescriptionsIndex2;
  int elementRight               = element2;
  int faceIndexRight             = faceIndex2;

  if (pos1(normalDirection) > pos2(normalDirection)) {
    cellDescriptionsIndexLeft  = cellDescriptionsIndex2;
    elementLeft                = element2;
    faceIndexLeft              = faceIndex2;

    cellDescriptionsIndexRight = cellDescriptionsIndex1;
    elementRight               = element1;
    faceIndexRight             = faceIndex1;
  }

  SolverPatch& solverPatchLeft  = _solver->getCellDescription(cellDescriptionsIndexLeft,elementLeft);
  SolverPatch& solverPatchRight = _solver->getCellDescription(cellDescriptionsIndexRight,elementRight);

  // 1.2. Merge min/max of both solver patches
  mergeSolutionMinMaxOnFace(solverPatchLeft,solverPatchRight,faceIndexLeft,faceIndexRight);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryData(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  mergeWithBoundaryDataBasedOnLimiterStatus(
      cellDescriptionsIndex,element,solverPatch.getLimiterStatus(),posCell,posBoundary,tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryDataBasedOnLimiterStatus(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const SolverPatch::LimiterStatus&         limiterStatus,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      _solver->mergeWithBoundaryData(cellDescriptionsIndex,element,posCell,posBoundary,
                                     tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
      break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      _limiter->mergeWithBoundaryData(cellDescriptionsIndex,limiterElement,posCell,posBoundary,
                                     tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
      break;
    default:
      break;
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch&  SolverPatch,
  int               faceIndex,
  double* min, double* max) const {
  if (SolverPatch.getType() == SolverPatch::Cell ||
      SolverPatch.getType() == SolverPatch::Ancestor ||
      SolverPatch.getType() == SolverPatch::Descendant
      ) {
    assertion( exahype::solvers::RegisteredSolvers[ SolverPatch.getSolverNumber() ]->getType()==exahype::solvers::Solver::Type::ADERDG );
    const int numberOfVariables = static_cast<exahype::solvers::ADERDGSolver*>(
        exahype::solvers::RegisteredSolvers[ SolverPatch.getSolverNumber() ])->getNumberOfVariables();

    for (int i=0; i<numberOfVariables; i++) {
      DataHeap::getInstance().getData( SolverPatch.getSolutionMin()  )[i+faceIndex*numberOfVariables]  =
        std::min( DataHeap::getInstance().getData( SolverPatch.getSolutionMin()  )[i+faceIndex*numberOfVariables], min[i] );
      DataHeap::getInstance().getData( SolverPatch.getSolutionMax()  )[i+faceIndex*numberOfVariables]  =
        std::max( DataHeap::getInstance().getData( SolverPatch.getSolutionMax()  )[i+faceIndex*numberOfVariables], max[i] );
    }
  }
}
