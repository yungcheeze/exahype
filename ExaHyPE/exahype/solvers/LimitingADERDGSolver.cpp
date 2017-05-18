/*
 *
 *  Created on: 26 Oct 2016
 *      Author: dominic
 */

#include "LimitingADERDGSolver.h"

#include "kernels/limiter/generic/Limiter.h"

#include "exahype/VertexOperations.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"

namespace exahype {
namespace solvers {

} /* namespace solvers */
} /* namespace exahype */

tarch::logging::Log exahype::solvers::LimitingADERDGSolver::_log("exahype::solvers::LimitingADERDGSolver");

bool exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalOrGlobalRecomputation() {
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG &&
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
        !=exahype::solvers::LimiterDomainChange::Regular
    ) {
      return true;
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalRecomputation(){
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG &&
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
        ==exahype::solvers::LimiterDomainChange::Irregular
    ) {
      return true;
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::oneSolverRequestedGlobalRecomputation(){
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG &&
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
        ==exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate
    ) {
      return true;
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::isValidCellDescriptionIndex(
    const int cellDescriptionsIndex) const  {
  return _solver->isValidCellDescriptionIndex(cellDescriptionsIndex);
}

exahype::solvers::Solver::SubcellPosition exahype::solvers::LimitingADERDGSolver::computeSubcellPositionOfCellOrAncestor(
      const int cellDescriptionsIndex,
      const int element) {
  return _solver->computeSubcellPositionOfCellOrAncestor(cellDescriptionsIndex,element);
}

exahype::solvers::LimitingADERDGSolver::LimitingADERDGSolver(
    const std::string& identifier,
    std::unique_ptr<exahype::solvers::ADERDGSolver> solver,
    std::unique_ptr<exahype::solvers::FiniteVolumesSolver> limiter,
    const double DMPRelaxationParameter,
    const double DMPDifferenceScaling,
    const int iterationsToCureTroubledCell)
    :
    exahype::solvers::Solver(identifier, Solver::Type::LimitingADERDG, solver->getNumberOfVariables(),
        solver->getNumberOfParameters(), solver->getNodesPerCoordinateAxis(), solver->getMaximumMeshSize(),
          solver->getMaximumAdaptiveMeshDepth(),
          solver->getTimeStepping()),
          _solver(std::move(solver)),
          _limiter(std::move(limiter)),
          _limiterDomainChange(LimiterDomainChange::Regular),
          _nextLimiterDomainChange(LimiterDomainChange::Regular),
          _DMPMaximumRelaxationParameter(DMPRelaxationParameter),
          _DMPDifferenceScaling(DMPDifferenceScaling),
          _iterationsToCureTroubledCell(iterationsToCureTroubledCell)
{
  assertion(_solver->getNumberOfParameters() == 0);
  assertion(_solver->getTimeStepping()==_limiter->getTimeStepping());
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
}

void exahype::solvers::LimitingADERDGSolver::initSolver(
    const double timeStamp,
    const tarch::la::Vector<DIMENSIONS,double>& domainOffset,
    const tarch::la::Vector<DIMENSIONS,double>& domainSize) {
  _domainOffset=domainOffset;
  _domainSize=domainSize;
  _coarsestMeshLevel =
      exahype::solvers::Solver::computeMeshLevel(_maximumMeshSize,domainSize[0]);

  _limiterDomainChange=LimiterDomainChange::Regular;

  _solver->initSolver(timeStamp, domainOffset, domainSize);
}

bool exahype::solvers::LimitingADERDGSolver::isCommunicating(exahype::records::State::AlgorithmSection& section) const {
  return
      section==exahype::records::State::AlgorithmSection::TimeStepping
      ||
      (section==exahype::records::State::AlgorithmSection::APosterioriRefinement && // assert: getMeshUpdateRequest())
      getLimiterDomainChange()==LimiterDomainChange::IrregularRequiringMeshUpdate)
      ||
      (section==exahype::records::State::AlgorithmSection::APrioriRefinement &&
       getMeshUpdateRequest() &&
       getLimiterDomainChange()==LimiterDomainChange::Regular)
      ||
      (section==exahype::records::State::AlgorithmSection::LocalRecomputation && // assert: !getMeshUpdateRequest()
      getLimiterDomainChange()==LimiterDomainChange::Regular);
}

bool exahype::solvers::LimitingADERDGSolver::isComputing(exahype::records::State::AlgorithmSection& section) const {
  return
      section==exahype::records::State::AlgorithmSection::TimeStepping
      ||
      (section==exahype::records::State::AlgorithmSection::APosterioriRefinement && // assert: getMeshUpdateRequest())
      getLimiterDomainChange()==LimiterDomainChange::IrregularRequiringMeshUpdate)
      ||
      (section==exahype::records::State::AlgorithmSection::APrioriRefinement &&
       getMeshUpdateRequest() &&
       getLimiterDomainChange()==LimiterDomainChange::Regular)
      ||
      (section==exahype::records::State::AlgorithmSection::LocalRecomputation && // assert: !getMeshUpdateRequest()
      getLimiterDomainChange()==LimiterDomainChange::Regular);
}

void exahype::solvers::LimitingADERDGSolver::synchroniseTimeStepping(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->synchroniseTimeStepping(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::startNewTimeStep() {
  _solver->startNewTimeStep();

  _minCellSize     = _nextMinCellSize;
  _maxCellSize     = _nextMaxCellSize;
  _nextMinCellSize = std::numeric_limits<double>::max();
  _nextMaxCellSize = -std::numeric_limits<double>::max(); // "-", min

  setNextMeshUpdateRequest();

  logDebug("startNewTimeStep()","_limiterDomainHasChanged="<<_limiterDomainChange<<",nextLimiterDomainChange="<<_nextLimiterDomainChange);

  _limiterDomainChange     = _nextLimiterDomainChange;
  _nextLimiterDomainChange = false;
}

void exahype::solvers::LimitingADERDGSolver::zeroTimeStepSizes() {
  _solver->zeroTimeStepSizes();
}

exahype::solvers::LimiterDomainChange exahype::solvers::LimitingADERDGSolver::getLimiterDomainChange() const {
  return _limiterDomainChange;
}

exahype::solvers::LimiterDomainChange exahype::solvers::LimitingADERDGSolver::getNextLimiterDomainChange() const {
  return _nextLimiterDomainChange;
}

void exahype::solvers::LimitingADERDGSolver::updateNextLimiterDomainChange(
    exahype::solvers::LimiterDomainChange limiterDomainChange) {
  _nextLimiterDomainChange =
      std::max( _nextLimiterDomainChange, limiterDomainChange );
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep() {
  _solver->rollbackToPreviousTimeStep();
}

void exahype::solvers::LimitingADERDGSolver::reinitialiseTimeStepData() {
  _solver->reinitialiseTimeStepData();
}

void exahype::solvers::LimitingADERDGSolver::reconstructStandardTimeSteppingDataAfterRollback() {
  _solver->reconstructStandardTimeSteppingDataAfterRollback();
}

void exahype::solvers::LimitingADERDGSolver::updateNextMinCellSize(double minCellSize) {
  _solver->updateNextMinCellSize(minCellSize);
  _nextMinCellSize = _solver->getNextMinCellSize();
}

void exahype::solvers::LimitingADERDGSolver::updateNextMaxCellSize(double maxCellSize) {
  _solver->updateNextMaxCellSize(maxCellSize);
  _nextMaxCellSize = _solver->getNextMaxCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getNextMinCellSize() const {
  return _solver->getNextMinCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getNextMaxCellSize() const {
  return _solver->getNextMaxCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getMinCellSize() const {
  return _solver->getMinCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getMaxCellSize() const {
  return _solver->getMaxCellSize();
}

///////////////////////////////////
// MODIFY CELL DESCRIPTION
///////////////////////////////////
bool exahype::solvers::LimitingADERDGSolver::markForRefinementBasedOnLimiterStatus(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber) {
  const int solverElement =
      _solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);

  if (solverElement!=exahype::solvers::Solver::NotFound) {
    SolverPatch& solverPatch =
        _solver->getCellDescription(fineGridCell.getCellDescriptionsIndex(),solverElement);
    bool refineFineGridCell =
        evaluateLimiterStatusBasedRefinementCriterion(
            fineGridCell.getCellDescriptionsIndex(),solverElement);
    if (refineFineGridCell) {
      solverPatch.setRefinementEvent(SolverPatch::RefinementEvent::RefiningRequested);
      return true;
    }

    if (solverPatch.getLimiterStatus()<
        computeMinimumLimiterStatusForRefinement(solverPatch.getLevel())) {
      return _solver->markForRefinement(
              fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
              coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
              fineGridPositionOfCell,
              initialGrid,
              solverNumber);
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::markForRefinement(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber)  {
  bool refineFineGridCell =
      markForRefinementBasedOnLimiterStatus(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
          coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
          fineGridPositionOfCell,
          initialGrid,
          solverNumber);

  const int fineGridSolverElement =
      _solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (fineGridSolverElement!=exahype::solvers::Solver::NotFound) {
    updateLimiterStatus(fineGridCell.getCellDescriptionsIndex(),fineGridSolverElement);
    deallocateLimiterPatchOnHelperCell(fineGridCell.getCellDescriptionsIndex(),fineGridSolverElement);
    ensureRequiredLimiterPatchIsAllocated(fineGridCell.getCellDescriptionsIndex(),fineGridSolverElement);

    vetoErasingChildrenRequestBasedOnLimiterStatus(
        fineGridCell.getCellDescriptionsIndex(),fineGridSolverElement,
        coarseGridCell.getCellDescriptionsIndex());
  }

  return refineFineGridCell;
}

bool exahype::solvers::LimitingADERDGSolver::updateStateInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber)  {
  bool refineFineGridCell =
      _solver->updateStateInEnterCell(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
          coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
          fineGridPositionOfCell,initialGrid,solverNumber);

  return refineFineGridCell;
}

void exahype::solvers::LimitingADERDGSolver::vetoErasingChildrenRequestBasedOnLimiterStatus(
    const int fineGridCellDescriptionsIndex,
    const int fineGridSolverElement,
    const int coarseGridCellDescriptionsIndex) const {
  SolverPatch& fineGridSolverPatch = _solver->getCellDescription(
          fineGridCellDescriptionsIndex,fineGridSolverElement);
  assertion1(coarseGridCellDescriptionsIndex==fineGridSolverPatch.getParentIndex(),fineGridSolverPatch.toString());

  if (
      fineGridSolverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::ErasingChildrenRequested
      ||
      fineGridSolverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::ChangeChildrenToDescendantsRequested
  ) {
    if (fineGridSolverPatch.getLimiterStatus()>=
        computeMinimumLimiterStatusForRefinement(fineGridSolverPatch.getLevel())-1) {
      fineGridSolverPatch.setRefinementEvent(SolverPatch::RefinementEvent::None);
    }
  }
}

bool exahype::solvers::LimitingADERDGSolver::updateStateInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber)  {
  const int solverElement =
      _solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  const int parentSolverElement =
      _solver->tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
  if (solverElement!=exahype::solvers::Solver::NotFound
      &&
      parentSolverElement!=exahype::solvers::Solver::NotFound
  ) {
    restrictLimiterStatus(
        fineGridCell.getCellDescriptionsIndex(),solverElement,
        coarseGridCell.getCellDescriptionsIndex(),parentSolverElement);
  }

  return _solver->updateStateInLeaveCell(
      fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
      coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
      fineGridPositionOfCell,solverNumber);
}

bool exahype::solvers::LimitingADERDGSolver::attainedStableState(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const int solverNumber) const {
  return _solver->attainedStableState(
      fineGridCell,fineGridVertices,fineGridVerticesEnumerator,solverNumber);
}

void exahype::solvers::LimitingADERDGSolver::finaliseStateUpdates(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) {
  _solver->finaliseStateUpdates(
      fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
      coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
      fineGridPositionOfCell,solverNumber);
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////
int exahype::solvers::LimitingADERDGSolver::computeMinimumLimiterStatusForRefinement(
    int level) const {
  const int levelDelta = getMaximumAdaptiveMeshLevel() - level;
  return std::min(static_cast<int>(SolverPatch::LimiterStatus::NeighbourOfTroubled2),
               static_cast<int>(SolverPatch::LimiterStatus::Ok) + levelDelta);
}

bool exahype::solvers::LimitingADERDGSolver::evaluateLimiterStatusBasedRefinementCriterion(
    const int cellDescriptionsIndex,const int solverElement) const {
  if (solverElement!=exahype::solvers::Solver::NotFound) {
    SolverPatch& solverPatch =
        _solver->getCellDescription(cellDescriptionsIndex,solverElement);
    if (
        solverPatch.getType()==SolverPatch::Type::Cell
        &&
        solverPatch.getLevel() < getMaximumAdaptiveMeshLevel()
        &&
        (solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None ||
         solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::DeaugmentingChildrenRequested ||
         solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::AugmentingRequested)
    ) {
      return solverPatch.getLimiterStatus() >
                computeMinimumLimiterStatusForRefinement(solverPatch.getLevel());
    }
  }
  return false;
}

bool
exahype::solvers::LimitingADERDGSolver::evaluateRefinementCriterionAfterSolutionUpdate(
      const int cellDescriptionsIndex,const int element) {
  // First evaluate the limiter status based refinement criterion
  bool refinementRequested =
      evaluateLimiterStatusBasedRefinementCriterion(cellDescriptionsIndex,element);

  // If no refinement was requested, evaluate the user's refinement criterion
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  if (!refinementRequested
      &&
      solverPatch.getLimiterStatus()<
      computeMinimumLimiterStatusForRefinement(solverPatch.getLevel())
  ) {
    refinementRequested =
        _solver->evaluateRefinementCriterionAfterSolutionUpdate(
            cellDescriptionsIndex,element);
  }

  return false;
}


double exahype::solvers::LimitingADERDGSolver::startNewTimeStep(
    const int cellDescriptionsIndex,
    const int element,
    double*   tempEigenvalues)  {
  double admissibleTimeStepSize =
      _solver->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
  return admissibleTimeStepSize;
}

void exahype::solvers::LimitingADERDGSolver::zeroTimeStepSizes(const int cellDescriptionsIndex, const int solverElement) {
  _solver->zeroTimeStepSizes(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep(
    const int cellDescriptionsIndex,
    const int solverElement) {
  _solver->rollbackToPreviousTimeStep(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::reconstructStandardTimeSteppingDataAfterRollback(
      const int cellDescriptionsIndex,
      const int element) const {
  _solver->reconstructStandardTimeSteppingDataAfterRollback(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::setInitialConditions(
    const int cellDescriptionsIndex,
    const int solverElement,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  _solver->setInitialConditions(
      cellDescriptionsIndex,solverElement,
      fineGridVertices,fineGridVerticesEnumerator);

  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    SolverPatch& solverPatch =
            _solver->getCellDescription(cellDescriptionsIndex,solverElement);
    LimiterPatch& limiterPatch =
            _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
    limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
    limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());

    _limiter->setInitialConditions(
        cellDescriptionsIndex,limiterElement,
        fineGridVertices,fineGridVerticesEnumerator);
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

  // 0. Erase old cells; now it's safe (TODO(Dominic): Add to docu)
  ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(cellDescriptionsIndex,element);
  deallocateLimiterPatchOnHelperCell(cellDescriptionsIndex,element);

  // 1. Write back the limiter status to the previous limiter status field
  solverPatch.setPreviousLimiterStatus(solverPatch.getLimiterStatus());

  // 2. Update the solution in the cells
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
      const int limiterElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());

      switch (solverPatch.getLimiterStatus()) {
      case SolverPatch::LimiterStatus::Ok:
        _solver->updateSolution(
            cellDescriptionsIndex,element,
            tempStateSizedVectors,tempUnknowns,
            fineGridVertices,fineGridVerticesEnumerator);
        break;
      case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled4: {
        _solver->updateSolution(
            cellDescriptionsIndex,element,
            tempStateSizedVectors,tempUnknowns,
            fineGridVertices,fineGridVerticesEnumerator);

        assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());
        LimiterPatch& limiterPatch =
            _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

        // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
        double* solverSolution = DataHeap::getInstance().getData(
            solverPatch.getSolution()).data();
        _limiter->swapSolutionAndPreviousSolution(cellDescriptionsIndex,limiterElement);
        double* limiterSolution = DataHeap::getInstance().getData(
            limiterPatch.getSolution()).data();
        kernels::limiter::generic::c::projectOnFVLimiterSpace(
            solverSolution,_solver->getNumberOfVariables(),
            _solver->getNodesPerCoordinateAxis(),
            _limiter->getGhostLayerWidth(),
            limiterSolution);
      } break;
      case SolverPatch::LimiterStatus::Troubled:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled2: {
        if (solverPatch.getType()==SolverPatch::Type::Cell) {
          assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());

          LimiterPatch& limiterPatch =
              _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

          limiterPatch.setPreviousTimeStepSize(solverPatch.getPreviousCorrectorTimeStepSize());
          limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
          limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());

          _limiter->updateSolution(
              cellDescriptionsIndex,limiterElement,
              tempStateSizedVectors,
              tempUnknowns,
              fineGridVertices,fineGridVerticesEnumerator);

          double* solverSolution = DataHeap::getInstance().getData(
              solverPatch.getSolution()).data();
          double* limiterSolution = DataHeap::getInstance().getData(
              limiterPatch.getSolution()).data();

          // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
          kernels::limiter::generic::c::projectOnDGSpace(
              limiterSolution,_solver->getNumberOfVariables(),
              _solver->getNodesPerCoordinateAxis(),
              _limiter->getGhostLayerWidth(),
              solverSolution);
        }
      } break;
      }
    } else {
      _solver->updateSolution(
          cellDescriptionsIndex,element,
          tempStateSizedVectors,tempUnknowns,
          fineGridVertices,fineGridVerticesEnumerator);
    }
  }
}

//void printSolutionMinOrMax(const double* minOrMax,const int numberOfVariables,const char* identifier) {
//  std::cout << identifier << "=" ;
//
//  for (int i = 0; i < DIMENSIONS_TIMES_TWO*numberOfVariables; ++i) {
//    std::cout << minOrMax[i] << ",";
//  }
//  std::cout << std::endl;
//}
//
//void printNormalFluxes(const double* flux,const int numberOfFaceUnknowns,const char* identifier) {
//  std::cout << identifier << "=" ;
//
//  for (int d=0; d<DIMENSIONS_TIMES_TWO;d++) {
//    for (int i = 0; i < numberOfFaceUnknowns; ++i) {
//      std::cout << flux[i+d*numberOfFaceUnknowns] << ",";
//    }
//    std::cout << "|";
//  }
//  std::cout << std::endl;
//}

exahype::solvers::LimiterDomainChange
exahype::solvers::LimitingADERDGSolver::updateLimiterStatusAndMinAndMaxAfterSolutionUpdate(
    const int cellDescriptionsIndex,
    const int solverElement) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);

  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    bool solutionIsValid =
        evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(solverPatch)
        && evaluatePhysicalAdmissibilityCriterion(solverPatch); // after min and max was found

    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
      switch (solverPatch.getLimiterStatus()) {
      case SolverPatch::LimiterStatus::Troubled:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled2: {
        const int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
        assertion2(limiterElement!=exahype::solvers::Solver::NotFound,limiterElement,cellDescriptionsIndex);
        LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
        determineLimiterMinAndMax(solverPatch,limiterPatch);
      }
      break;
      case SolverPatch::LimiterStatus::Ok:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
        // Keep the previously computed min and max values
        break;
      }
    }

    LimiterDomainChange limiterDomainChange =
        determineLimiterStatusAfterSolutionUpdate(solverPatch,!solutionIsValid);
    allocateLimiterPatchAfterSolutionUpdate(cellDescriptionsIndex,solverElement);
    return limiterDomainChange;
  } else {
    return updateLimiterStatus(cellDescriptionsIndex,solverElement);
  }
}

void exahype::solvers::LimitingADERDGSolver::allocateLimiterPatchAfterSolutionUpdate(
    const int cellDescriptionsIndex,const int solverElement) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);

  bool limiterPatchAllocated =
      ensureRequiredLimiterPatchIsAllocated(cellDescriptionsIndex,solverPatch.getSolverNumber());

  if (limiterPatchAllocated) {
    switch (solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4: {
      const int limiterElement =
          tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
      assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());
      LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
    } break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2: {
      // TODO(Dominic): This should actually never be entered unless for the case Troubled
    } break;
    }
  }
}


exahype::solvers::LimiterDomainChange
exahype::solvers::LimitingADERDGSolver::updateLimiterStatusAndMinAndMaxAfterSetInitialConditions(
    const int cellDescriptionsIndex,
    const int solverElement) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);

  if (tarch::la::equals(solverPatch.getCorrectorTimeStepSize(),0.0) &&
      solverPatch.getType()==SolverPatch::Type::Cell) {
    determineSolverMinAndMax(solverPatch);
    return determineLimiterStatusAfterSolutionUpdate(
            solverPatch,
            !evaluatePhysicalAdmissibilityCriterion(solverPatch)); // only evaluate PAD here
  }

  return false;
}

exahype::solvers::LimiterDomainChange
exahype::solvers::LimitingADERDGSolver::determineLimiterStatusAfterSolutionUpdate(
    SolverPatch& solverPatch,const bool isTroubled) const {
  assertion1(solverPatch.getType()==SolverPatch::Type::Cell,solverPatch.toString());

  LimiterDomainChange limiterDomainChange = LimiterDomainChange::Regular;

  if (isTroubled) {
    solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
    switch (solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Troubled:
      // do nothing
      break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
      solverPatch.setLimiterStatus(SolverPatch::LimiterStatus::Troubled);
      ADERDGSolver::writeLimiterStatusOnBoundary(solverPatch);
      break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
    case SolverPatch::LimiterStatus::Ok:
      limiterDomainChange = LimiterDomainChange::IrregularRequiringMeshUpdate;
      if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
        limiterDomainChange = LimiterDomainChange::Irregular;
      }
      solverPatch.setLimiterStatus(SolverPatch::LimiterStatus::Troubled);
      ADERDGSolver::writeLimiterStatusOnBoundary(solverPatch);
      break;
    }
  } else {
    switch (solverPatch.getPreviousLimiterStatus()) {
    case SolverPatch::LimiterStatus::Troubled:
      solverPatch.setLimiterStatus(SolverPatch::LimiterStatus::Troubled);
      solverPatch.setIterationsToCureTroubledCell(
          solverPatch.getIterationsToCureTroubledCell()-1);
      if (solverPatch.getIterationsToCureTroubledCell()==0) {
        solverPatch.setLimiterStatus(SolverPatch::LimiterStatus::NeighbourOfTroubled1);
        solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell); // TODO(Dominic): Probably not necessary
      }
      ADERDGSolver::writeLimiterStatusOnBoundary(solverPatch);
      break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
    case SolverPatch::LimiterStatus::Ok:
      // do nothing
      break;
    }
  }
  solverPatch.setLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));
  ADERDGSolver::writeLimiterStatusOnBoundary(solverPatch);

  return limiterDomainChange;
}

bool exahype::solvers::LimitingADERDGSolver::evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(SolverPatch& solverPatch) {
  double* solution = DataHeap::getInstance().getData(
      solverPatch.getSolution()).data();

  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    double* observablesMin = DataHeap::getInstance().getData(
        solverPatch.getSolutionMin()).data();
    double* observablesMax = DataHeap::getInstance().getData(
        solverPatch.getSolutionMax()).data();

    // 1. Check if the DMP is satisfied and search for the min and max
    // Write the new min and max to the storage reserved for face 0
    bool dmpIsSatisfied = kernels::limiter::generic::c::discreteMaximumPrincipleAndMinAndMaxSearch(
          solution,_solver.get(),
          _DMPMaximumRelaxationParameter, _DMPDifferenceScaling,
          observablesMin,observablesMax);

    // 2. Copy the result on the other faces as well
    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
      std::copy_n(
          observablesMin,numberOfObservables, // past-the-end element
          observablesMin+i*numberOfObservables);
      std::copy_n(
          observablesMax,numberOfObservables, // past-the-end element
          observablesMax+i*numberOfObservables);
    }

    return dmpIsSatisfied;
  } else {
    return true;
  }
}

bool exahype::solvers::LimitingADERDGSolver::evaluatePhysicalAdmissibilityCriterion(SolverPatch& solverPatch) {
  double* observablesMin = nullptr;
  double* observablesMax = nullptr;

  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables > 0) {
    observablesMin = DataHeap::getInstance().getData(
        solverPatch.getSolutionMin()).data();
    observablesMax = DataHeap::getInstance().getData(
        solverPatch.getSolutionMax()).data();
  }

  const double* const solution = DataHeap::getInstance().getData(
        solverPatch.getSolution()).data();

  return _solver->isPhysicallyAdmissible(
      solution,
      observablesMin,observablesMax,numberOfObservables,
      solverPatch.getOffset()+0.5*solverPatch.getSize(),solverPatch.getSize(),
      solverPatch.getCorrectorTimeStamp(),solverPatch.getCorrectorTimeStepSize());
}

void exahype::solvers::LimitingADERDGSolver::determineMinAndMax(
    const int cellDescriptionsIndex,
    const int element) {
  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    switch(solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
      determineSolverMinAndMax(solverPatch);
      break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
      assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());
      exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      determineLimiterMinAndMax(solverPatch,limiterPatch);
      break;
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::determineSolverMinAndMax(SolverPatch& solverPatch) {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    const double* const solution = DataHeap::getInstance().getData(
        solverPatch.getSolution()).data();

    double* solutionMin = DataHeap::getInstance().getData(
        solverPatch.getSolutionMin()).data();
    double* solutionMax = DataHeap::getInstance().getData(
        solverPatch.getSolutionMax()).data();

    // Write the result to the face with index "0"
    kernels::limiter::generic::c::findCellLocalMinAndMax(
        solution,_solver.get(),
        solutionMin,solutionMax);

    // Copy the result on the other faces as well
    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
      std::copy_n(
          solutionMin,numberOfObservables, // past-the-end element
          solutionMin+i*numberOfObservables);
      std::copy_n(
          solutionMax,numberOfObservables, // past-the-end element
          solutionMax+i*numberOfObservables);
    }

    for (int i=0; i<DIMENSIONS_TIMES_TWO*numberOfObservables; ++i) {
      assertion(*(solutionMin+i)<std::numeric_limits<double>::max());
      assertion(*(solutionMax+i)>-std::numeric_limits<double>::max());
    } // Dead code elimination will get rid of this loop
  }
}

void exahype::solvers::LimitingADERDGSolver::determineLimiterMinAndMax(SolverPatch& solverPatch,LimiterPatch& limiterPatch) {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    double* limiterSolution = DataHeap::getInstance().getData(
        limiterPatch.getSolution()).data();

    double* observablesMin = DataHeap::getInstance().getData(
        solverPatch.getSolutionMin()).data();
    double* observablesMax = DataHeap::getInstance().getData(
        solverPatch.getSolutionMax()).data();

    // Write the result to the face with index "0"
    kernels::limiter::generic::c::findCellLocalLimiterMinAndMax(
        limiterSolution,_solver.get(),
        _limiter->getGhostLayerWidth(),observablesMin,observablesMax);

    // Copy the result on the other faces as well
    const int numberOfObservables = _solver->getDMPObservables();
    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
      std::copy_n(
          observablesMin,numberOfObservables, // past-the-end element
          observablesMin+i*numberOfObservables);
      std::copy_n(
          observablesMax,numberOfObservables, // past-the-end element
          observablesMax+i*numberOfObservables);
    }

    for (int i=0; i<DIMENSIONS_TIMES_TWO*numberOfObservables; ++i) {
      assertion(*(observablesMin+i)<std::numeric_limits<double>::max());
      assertion(*(observablesMax+i)>-std::numeric_limits<double>::max());
    } // Dead code elimination will get rid of this loop
  }
}

int exahype::solvers::LimitingADERDGSolver::convertToCoarserMeshLevelLimiterStatus(
    const int limiterStatusAsInt) const {
  return std::max(0, std::min(3, limiterStatusAsInt-2));
}

void exahype::solvers::LimitingADERDGSolver::mergeWithLimiterStatus(
    SolverPatch& solverPatch,
    const int direction,
    const int otherLimiterStatusAsInt) const {
  int limiterStatusAsInt = solverPatch.getLimiterStatus();
  limiterStatusAsInt =
      ( limiterStatusAsInt>=static_cast<int>(SolverPatch::LimiterStatus::Troubled) )
      ?
      limiterStatusAsInt
      :
      std::max (
          0,
          std::max( limiterStatusAsInt, otherLimiterStatusAsInt ) - 1
      );

  solverPatch.setFacewiseLimiterStatus( direction, limiterStatusAsInt );
}

void exahype::solvers::LimitingADERDGSolver::deallocateLimiterPatch(
        const int cellDescriptionsIndex,
        const int solverElement) const {
  assertion(solverElement!=exahype::solvers::Solver::NotFound);
  const int limiterElement =
        tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  assertion(limiterElement!=exahype::solvers::Solver::NotFound);
  LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
  limiterPatch.setType(LimiterPatch::Type::Erased);
  _limiter->ensureNoUnnecessaryMemoryIsAllocated(limiterPatch);

  tarch::multicore::Lock lock(_heapSemaphore);
  LimiterHeap::getInstance().getData(cellDescriptionsIndex).erase(
      LimiterHeap::getInstance().getData(cellDescriptionsIndex).begin()+limiterElement);
  lock.free();
}

void exahype::solvers::LimitingADERDGSolver::deallocateLimiterPatchOnHelperCell(
    const int cellDescriptionsIndex,
    const int solverElement) const {
  assertion(solverElement!=exahype::solvers::Solver::NotFound);
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);
  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);

  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    switch(solverPatch.getType()) {
    case SolverPatch::Type::EmptyAncestor:
    case SolverPatch::Type::Ancestor:
    case SolverPatch::Type::Descendant:
    case SolverPatch::Type::EmptyDescendant:
    case SolverPatch::Type::Erased: {
      deallocateLimiterPatch(cellDescriptionsIndex,solverElement);
    } break;
    case SolverPatch::Type::Cell: {
      // do nothing
    } break;
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(
    const int cellDescriptionsIndex,
    const int solverElement) const {
  assertion(solverElement!=exahype::solvers::Solver::NotFound);
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);
  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);

  if (limiterElement!=exahype::solvers::Solver::NotFound
      && solverPatch.getLimiterStatus()==SolverPatch::LimiterStatus::Ok) {
    deallocateLimiterPatch(cellDescriptionsIndex,solverElement);
  }
}

int exahype::solvers::LimitingADERDGSolver::allocateLimiterPatch(
        const int cellDescriptionsIndex,
        const int solverElement) const {
  assertion(solverElement!=exahype::solvers::Solver::NotFound);
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);

  #if defined(Asserts)
  const int previouslimiterElement =
          tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  #endif
  assertion(previouslimiterElement==exahype::solvers::Solver::NotFound);

  tarch::multicore::Lock lock(_heapSemaphore);
  exahype::solvers::FiniteVolumesSolver::addNewCellDescription(
      cellDescriptionsIndex,
      solverPatch.getSolverNumber(),
      LimiterPatch::Type::Cell,
      LimiterPatch::RefinementEvent::None,
      solverPatch.getLevel(),
      solverPatch.getParentIndex(),
      solverPatch.getSize(),
      solverPatch.getOffset());
  lock.free();

  assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getPreviousSolution()),solverPatch.toString());
  assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolution()),solverPatch.toString());

  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  assertion(limiterElement!=exahype::solvers::Solver::NotFound);
  LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
  _limiter->ensureNecessaryMemoryIsAllocated(limiterPatch);

  // Copy more geometry information from the solver patch
  limiterPatch.setIsInside(solverPatch.getIsInside());

  return limiterElement;
}

bool exahype::solvers::LimitingADERDGSolver::ensureRequiredLimiterPatchIsAllocated(
        const int cellDescriptionsIndex,
        const int solverElement) const {
  assertion(solverElement!=exahype::solvers::Solver::NotFound);
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);
  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()
      && limiterElement==exahype::solvers::Solver::NotFound
      && solverPatch.getType()==SolverPatch::Type::Cell
  ) {
    switch (solverPatch.getLimiterStatus()) {
      case SolverPatch::LimiterStatus::Troubled:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled4: {
        assertion1(solverPatch.getPreviousLimiterStatus()==SolverPatch::LimiterStatus::Ok,solverPatch.toString());
        allocateLimiterPatch(cellDescriptionsIndex,solverElement);
      } return true;
    }
  }
  return false;
}


void exahype::solvers::LimitingADERDGSolver::projectDGSolutionOnFVSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
  double*       limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution()).data();

  // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
  kernels::limiter::generic::c::projectOnFVLimiterSpace(
      solverSolution,_solver->getNumberOfVariables(),
      _solver->getNodesPerCoordinateAxis(),
      _limiter->getGhostLayerWidth(),
      limiterSolution);
}

void exahype::solvers::LimitingADERDGSolver::reinitialiseSolvers(
    const int cellDescriptionsIndex,
    const int solverElement,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);

  // 0. Update the limiter status (do not overwrite the previous limiter status)
  solverPatch.setLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));

  // 1. Allocate or deallocate a limiter patch
  deallocateLimiterPatchOnHelperCell(cellDescriptionsIndex,solverElement);
  ensureRequiredLimiterPatchIsAllocated(cellDescriptionsIndex,solverElement);
  ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(cellDescriptionsIndex,solverElement);

  // 2. Rollback with limiter or solver solution depending on limiter status
  if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getType()==SolverPatch::Type::Cell) {
    switch (solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4: {
      switch (solverPatch.getPreviousLimiterStatus()) {
      case SolverPatch::LimiterStatus::Troubled:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled2: {
        const int limiterElement =
            tryGetLimiterElementFromSolverElement(fineGridCell.getCellDescriptionsIndex(),solverElement);
        assertion(limiterElement!=exahype::solvers::Solver::NotFound);
        _limiter->rollbackSolution(
            fineGridCell.getCellDescriptionsIndex(),limiterElement,
            fineGridVertices,fineGridVerticesEnumerator);
      } break;
      case SolverPatch::LimiterStatus::Ok:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled4: {
        _solver->rollbackSolution(
            fineGridCell.getCellDescriptionsIndex(),solverElement,
            fineGridVertices,fineGridVerticesEnumerator);

        if (solverPatch.getPreviousLimiterStatus()==SolverPatch::LimiterStatus::Ok) {
          const int limiterElement =
              tryGetLimiterElementFromSolverElement(fineGridCell.getCellDescriptionsIndex(),solverElement);
          assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());
          assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getPreviousSolution()),solverPatch.toString());
          assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolution()),solverPatch.toString());

          LimiterPatch& limiterPatch = _limiter->getCellDescription(fineGridCell.getCellDescriptionsIndex(),limiterElement);
          projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
        } else {
          const int limiterElement =
              tryGetLimiterElementFromSolverElement(fineGridCell.getCellDescriptionsIndex(),solverElement);
          _limiter->rollbackSolution(cellDescriptionsIndex,limiterElement,fineGridVertices,fineGridVerticesEnumerator);
        }
      } break;
      }
    } break;
    case SolverPatch::LimiterStatus::Ok: {
      if (getLimiterDomainChange()==LimiterDomainChange::IrregularRequiringMeshUpdate) {
        _solver->rollbackSolution(
            fineGridCell.getCellDescriptionsIndex(),solverElement,
            fineGridVertices,fineGridVerticesEnumerator);
      }
      #if defined(Asserts)
      const int limiterElement =
          tryGetLimiterElementFromSolverElement(fineGridCell.getCellDescriptionsIndex(),solverElement);
      #endif
      assertion(limiterElement==exahype::solvers::Solver::NotFound);
    } break;
    }
  }

  // 3. Reset the iterationsToCure on all troubled cells to maximum value if cell is troubled
  if (solverPatch.getLimiterStatus()==SolverPatch::LimiterStatus::Troubled) {
    solverPatch.setIterationsToCureTroubledCell(1+_iterationsToCureTroubledCell);
  }
}

void exahype::solvers::LimitingADERDGSolver::projectFVSolutionOnDGSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution()).data();
  double*       solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

  // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
  kernels::limiter::generic::c::projectOnDGSpace(
      limiterSolution,_solver->getNumberOfVariables(),
      _solver->getNodesPerCoordinateAxis(),
      _limiter->getGhostLayerWidth(),
      solverSolution);
}

void exahype::solvers::LimitingADERDGSolver::recomputeSolution(
        const int cellDescriptionsIndex, const int solverElement,
        exahype::solvers::SolutionUpdateTemporaryVariables& solutionUpdateTemporaryVariables,
        exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);

  if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getType()==SolverPatch::Type::Cell) {
    switch (solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2: {
      const int limiterElement =
          tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      // Copy time step data from the solver patch
      limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
      limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
      // 1. Evolve solution to desired  time step again
      _limiter->updateSolution(cellDescriptionsIndex,limiterElement,
          solutionUpdateTemporaryVariables._tempStateSizedVectors[solverPatch.getSolverNumber()],
          solutionUpdateTemporaryVariables._tempUnknowns[solverPatch.getSolverNumber()],
          fineGridVertices,fineGridVerticesEnumerator);
      // 2. Project FV solution on ADER-DG space
      projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
    } break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4: { // TODO(Dominic): Add to docu: Here, were just went back one time step to supply the NT neighbours with old limiter unknowns.
      switch (solverPatch.getPreviousLimiterStatus()) {
      case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
      case SolverPatch::LimiterStatus::Ok: {
        const int limiterElement =
            tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
        assertion(limiterElement!=exahype::solvers::Solver::NotFound);
        _solver->swapSolutionAndPreviousSolution(cellDescriptionsIndex,solverElement);
        LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
        projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
      } break;
      case SolverPatch::LimiterStatus::Troubled: // TODO(Dominic): Update docu
      case SolverPatch::LimiterStatus::NeighbourOfTroubled1:   // TODO(Dominic): Update docu
      case SolverPatch::LimiterStatus::NeighbourOfTroubled2: { // TODO(Dominic): Update docu
        const int limiterElement =
                    tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
        assertion(limiterElement!=exahype::solvers::Solver::NotFound);
        _limiter->swapSolutionAndPreviousSolution(cellDescriptionsIndex,limiterElement);
        LimiterPatch& limiterPatch    = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
        projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
      } break;
      }
    } break;
    case SolverPatch::LimiterStatus::Ok: {
      #if defined(Asserts)
      const int limiterElement =
          tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
      #endif
      assertion(limiterElement==exahype::solvers::Solver::NotFound);
      // do nothing
    }  break;
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::recomputePredictor(
    const int cellDescriptionsIndex,
    const int element,
    exahype::solvers::PredictionTemporaryVariables& predictionTemporaryVariables,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()
      && solverPatch.getType()==SolverPatch::Type::Cell
      && solverPatch.getPredictorTimeStepSize() > 0) {
    switch (solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Troubled:
      // do nothing
      break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
      _solver->performPredictionAndVolumeIntegral(
          solverPatch,
          predictionTemporaryVariables._tempSpaceTimeUnknowns    [solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempSpaceTimeFluxUnknowns[solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempUnknowns             [solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempFluxUnknowns         [solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempStateSizedVectors    [solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempPointForceSources    [solverPatch.getSolverNumber()]);
      break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
    case SolverPatch::LimiterStatus::Ok:
      switch (solverPatch.getPreviousLimiterStatus()) {
        case SolverPatch::LimiterStatus::Troubled:
        _solver->performPredictionAndVolumeIntegral(
            solverPatch,
            predictionTemporaryVariables._tempSpaceTimeUnknowns    [solverPatch.getSolverNumber()],
            predictionTemporaryVariables._tempSpaceTimeFluxUnknowns[solverPatch.getSolverNumber()],
            predictionTemporaryVariables._tempUnknowns             [solverPatch.getSolverNumber()],
            predictionTemporaryVariables._tempFluxUnknowns         [solverPatch.getSolverNumber()],
            predictionTemporaryVariables._tempStateSizedVectors    [solverPatch.getSolverNumber()],
            predictionTemporaryVariables._tempPointForceSources    [solverPatch.getSolverNumber()]);
        break;
      default:
        break;
      }
      break;
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::prepareNextNeighbourMerging(
    const int cellDescriptionsIndex,const int solverElement,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  _solver->prepareNextNeighbourMerging(
      cellDescriptionsIndex,solverElement,
      fineGridVertices,fineGridVerticesEnumerator);

  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->prepareNextNeighbourMerging(
        cellDescriptionsIndex,limiterElement,
        fineGridVertices,fineGridVerticesEnumerator);
  }
}

exahype::solvers::LimiterDomainChange exahype::solvers::LimitingADERDGSolver::updateLimiterStatus(
    const int cellDescriptionsIndex,const int solverElement) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);
  solverPatch.setLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));
  ADERDGSolver::writeLimiterStatusOnBoundary(solverPatch);

  if (
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()
      &&
      solverPatch.getLimiterStatus()>SolverPatch::LimiterStatus::Ok
      &&
      (solverPatch.getType()==SolverPatch::EmptyDescendant ||
       solverPatch.getType()==SolverPatch::Descendant)) {
    return LimiterDomainChange::IrregularRequiringMeshUpdate;
  }

  return LimiterDomainChange::Regular;
}

void exahype::solvers::LimitingADERDGSolver::preProcess(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->preProcess(cellDescriptionsIndex,element);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  _limiter->preProcess(cellDescriptionsIndex,limiterElement);
}

void exahype::solvers::LimitingADERDGSolver::postProcess(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->postProcess(cellDescriptionsIndex,element);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  _limiter->postProcess(cellDescriptionsIndex,limiterElement);
}

void exahype::solvers::LimitingADERDGSolver::prolongateDataAndPrepareDataRestriction(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->prolongateDataAndPrepareDataRestriction(cellDescriptionsIndex,element);

  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  if (
      (solverPatch.getType()==SolverPatch::Type::EmptyAncestor ||
      solverPatch.getType()==SolverPatch::Type::Ancestor)
      &&
      solverPatch.getPreviousLimiterStatus()==SolverPatch::LimiterStatus::Troubled) {
    solverPatch.setLimiterStatus(SolverPatch::LimiterStatus::NeighbourOfTroubled1);
    solverPatch.setFacewiseLimiterStatus(SolverPatch::LimiterStatus::NeighbourOfTroubled1);
  }
}

void exahype::solvers::LimitingADERDGSolver::restrictLimiterStatus(
    const int cellDescriptionsIndex,
    const int element,
    const int parentCellDescriptionsIndex,
    const int parentElement) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  if (solverPatch.getLimiterStatus()==SolverPatch::LimiterStatus::Troubled) {
    SolverPatch& parentSolverPatch =
        _solver->getCellDescription(parentCellDescriptionsIndex,parentElement);
    if(parentSolverPatch.getType()==SolverPatch::Type::Ancestor ||
       parentSolverPatch.getType()==SolverPatch::Type::EmptyAncestor) {
      parentSolverPatch.setLimiterStatus(SolverPatch::LimiterStatus::Troubled);
      ADERDGSolver::writeLimiterStatusOnBoundary(parentSolverPatch);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::restrictToNextParent(
        const int fineGridCellDescriptionsIndex,
        const int fineGridElement,
        const int coarseGridCellDescriptionsIndex,
        const int coarseGridElement) {
  restrictLimiterStatus(
      fineGridCellDescriptionsIndex,
      fineGridElement,coarseGridCellDescriptionsIndex,coarseGridElement);
}

void exahype::solvers::LimitingADERDGSolver::restrictToTopMostParent(
    const int cellDescriptionsIndex,
    const int solverElement,
    const int parentCellDescriptionsIndex,
    const int parentSolverElement,
    const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  _solver->restrictToTopMostParent(cellDescriptionsIndex,solverElement,parentCellDescriptionsIndex,parentSolverElement,subcellIndex);
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
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  // 1. Solve the riemann problems
  mergeNeighboursBasedOnLimiterStatus(
      cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,
      pos1,pos2,
      false, /* isRecomputation */
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);

  // 2. Merge the min and max of both cell description's solver's
  // solution value.
  mergeSolutionMinMaxOnFace(
      cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);

  // 3. Merge the limiter status
  mergeNeighboursLimiterStatus(
      cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
}

// TODO(Dominic): Remove limiterstatus1 and limiterStatus2 argument.
// They depend on the isRecomputation value
void exahype::solvers::LimitingADERDGSolver::mergeNeighboursBasedOnLimiterStatus(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    const bool                                isRecomputation,
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) const {
  SolverPatch& solverPatch1 = _solver->getCellDescription(cellDescriptionsIndex1,element1);
  SolverPatch& solverPatch2 = _solver->getCellDescription(cellDescriptionsIndex2,element2);
  const int limiterElement1 = tryGetLimiterElement(cellDescriptionsIndex1,solverPatch1.getSolverNumber());
  const int limiterElement2 = tryGetLimiterElement(cellDescriptionsIndex2,solverPatch2.getSolverNumber());

  // We only limit on the finest mesh level
  if (solverPatch1.getLevel()==getMaximumAdaptiveMeshLevel()) {
    // 1. Merge solver solution or limiter solution values in
    // non-overlapping parts of solver and limiter domain:
    switch (solverPatch1.getLimiterStatus()) {
      case SolverPatch::LimiterStatus::Ok:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
        switch (solverPatch2.getLimiterStatus()) {
          case SolverPatch::LimiterStatus::Ok:
          case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
          case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
            if (!isRecomputation) {
              _solver->mergeNeighbours(
                  cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
                  tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
            }
            break;
          default:
            break;
        }
        break;
          case SolverPatch::LimiterStatus::Troubled:
          case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
          case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
            switch (solverPatch2.getLimiterStatus()) {
              case SolverPatch::LimiterStatus::Troubled:
              case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
              case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
                _limiter->mergeNeighbours(
                    cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
                    tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
                break;
              default:
                break;
            }
            break;
    }
    // 2. Merge limiter solution values in overlapping part
    // of solver and limiter domain:
    switch (solverPatch1.getLimiterStatus()) {
      case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
        switch (solverPatch2.getLimiterStatus()) {
          case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
          case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
            assertion2(limiterElement1!=exahype::solvers::Solver::NotFound,solverPatch1.toString(),solverPatch2.toString());
            assertion2(limiterElement2!=exahype::solvers::Solver::NotFound,solverPatch2.toString(),solverPatch1.toString());
            _limiter->mergeNeighbours(
                cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
                tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.

            if (!isRecomputation) {
              _solver->mergeNeighbours(
                  cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
                  tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
            }
            break;
          default:
            break;
        }
        break;
          case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
          case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
            switch (solverPatch2.getLimiterStatus()) {
              case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
              case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
                assertion2(limiterElement1!=exahype::solvers::Solver::NotFound,solverPatch1.toString(),solverPatch2.toString());
                assertion2(limiterElement2!=exahype::solvers::Solver::NotFound,solverPatch2.toString(),solverPatch1.toString());
                _limiter->mergeNeighbours(
                    cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
                    tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
                if (!isRecomputation) {
                  _solver->mergeNeighbours(
                      cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
                      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
                }
                break;
              default:
                break;
            }
            break;
          default:
            break;
    }
  // On the other levels, we work with the ADER-DG solver only
  } else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
    if (!isRecomputation) {
      _solver->mergeNeighbours(
          cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
          tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  SolverPatch& solverPatch1 = _solver->getCellDescription(cellDescriptionsIndex1,element1);
  SolverPatch& solverPatch2 = _solver->getCellDescription(cellDescriptionsIndex2,element2);

  // 1.2. Merge min/max of both solver patches
  const int direction    = tarch::la::equalsReturnIndex(pos1,pos2);
  const int orientation1 = (1 + pos2(direction) - pos1(direction))/2;
  const int orientation2 = 1-orientation1;

  mergeSolutionMinMaxOnFace(
      solverPatch1,solverPatch2,
      2*direction+orientation1,2*direction+orientation2);
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch& solverPatch1,
  SolverPatch& solverPatch2,
  const int faceIndex1,
  const int faceIndex2
) const {
  if (
      _solver->getDMPObservables() > 0
      &&
      (solverPatch1.getType()==SolverPatch::Cell ||
      solverPatch2.getType()==SolverPatch::Cell)
  ) {
    assertion( solverPatch1.getSolverNumber() == solverPatch2.getSolverNumber() );
    const int numberOfObservables = _solver->getDMPObservables();
    double* min1 = DataHeap::getInstance().getData( solverPatch1.getSolutionMin()  ).data() + faceIndex1 * numberOfObservables;
    double* min2 = DataHeap::getInstance().getData( solverPatch2.getSolutionMin()  ).data() + faceIndex2 * numberOfObservables;
    double* max1 = DataHeap::getInstance().getData( solverPatch1.getSolutionMax()  ).data() + faceIndex1 * numberOfObservables;
    double* max2 = DataHeap::getInstance().getData( solverPatch2.getSolutionMax()  ).data() + faceIndex2 * numberOfObservables;

    for (int i=0; i<numberOfObservables; i++) {
      const double min = std::min(
          *(min1+i),
          *(min2+i)
      );
      const double max = std::max(
          *(max1+i),
          *(max2+i)
      );

      *(min1+i)  = min;
      *(min2+i) = min;

      *(max1+i)  = max;
      *(max2+i) = max;
    }
  } // else do nothing
}

void exahype::solvers::LimitingADERDGSolver::mergeNeighboursLimiterStatus(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) const {
  SolverPatch& solverPatch1 = _solver->getCellDescription(cellDescriptionsIndex1,element1);
  SolverPatch& solverPatch2 = _solver->getCellDescription(cellDescriptionsIndex2,element2);

  const int direction    = tarch::la::equalsReturnIndex(pos1,pos2);
  const int orientation1 = (1 + pos2(direction) - pos1(direction))/2;
  const int orientation2 = 1-orientation1;

  const int limiterStatus1 = solverPatch1.getLimiterStatus();
  const int limiterStatus2 = solverPatch2.getLimiterStatus();

  mergeWithLimiterStatus(solverPatch1,2*direction+orientation1,limiterStatus2);
  mergeWithLimiterStatus(solverPatch2,2*direction+orientation2,limiterStatus1);
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
      cellDescriptionsIndex,element,
      solverPatch.getLimiterStatus(),
      posCell,posBoundary,
      false, /* isRecomputation */
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryDataBasedOnLimiterStatus(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const int                                 limiterStatus, //TODO(Dominic): Still necessary?
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      const bool                                isRecomputation,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());

      switch (limiterStatus) {
      case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
      case SolverPatch::LimiterStatus::Ok:
        assertion(limiterStatus==SolverPatch::LimiterStatus::Ok || limiterElement!=exahype::solvers::Solver::NotFound);
        if (!isRecomputation) {
          _solver->mergeWithBoundaryData(cellDescriptionsIndex,element,posCell,posBoundary,
              tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
        }
        break;
      case SolverPatch::LimiterStatus::Troubled:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
        assertion(limiterElement!=exahype::solvers::Solver::NotFound);
        _limiter->mergeWithBoundaryData(cellDescriptionsIndex,limiterElement,posCell,posBoundary,
            tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
        break;
      default:
        break;
      }
    } else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
      if (!isRecomputation) {
        _solver->mergeWithBoundaryData(cellDescriptionsIndex,element,posCell,posBoundary,
            tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
      }
    }
  }
}

#ifdef Parallel
const int exahype::solvers::LimitingADERDGSolver::DataMessagesPerNeighbourCommunication    = 1;
const int exahype::solvers::LimitingADERDGSolver::DataMessagesPerForkOrJoinCommunication   = 0;
const int exahype::solvers::LimitingADERDGSolver::DataMessagesPerMasterWorkerCommunication = 0;

///////////////////////////////////
// NEIGHBOUR - Mesh refinement
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourMetadata(
      const exahype::MetadataHeap::HeapEntries& neighbourMetadata,
      const tarch::la::Vector<DIMENSIONS, int>& src,
      const tarch::la::Vector<DIMENSIONS, int>& dest,
      const int cellDescriptionsIndex,
      const int element) {
  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithNeighbourMetadata(neighbourMetadata,src,dest,cellDescriptionsIndex,limiterElement);
  } // There is no drop method for metadata necessary
  _solver->mergeWithNeighbourMetadata(neighbourMetadata,src,dest,cellDescriptionsIndex,element);

  // Refine according to the merged limiter status
  const SolverPatch::LimiterStatus neighbourLimiterStatus =
      neighbourMetadata[exahype::MetadataLimiterStatus].getU();
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  // TODO(Dominic): Update
  if (
      neighbourLimiterStatus==SolverPatch::LimiterStatus::Troubled
      &&
      solverPatch.getLevel()<_coarsestMeshLevel+_maximumAdaptiveMeshDepth
      &&
      solverPatch.getType()==SolverPatch::Type::Cell
      &&
      (solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None ||
      solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::DeaugmentingChildrenRequested ||
      solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::AugmentingRequested)
  ) {
    solverPatch.setRefinementEvent(SolverPatch::RefinementEvent::RefiningRequested);
  }

  // TODO(Dominic)
  if (tarch::la::equals(src,dest)==DIMENSIONS-1) {
    const int direction = peano::utils::dLinearised(src-dest+1,3);
    mergeWithLimiterStatus(solverPatch,direction,neighbourLimiterStatus);

    // old code keep for reference in case we go back to face index
    // const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
    // assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
    // const int faceIndex = 2 * normalOfExchangedFace +
    //        (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!
    //                                                                            // |src|dest| : 1; |dest|src| : 0
  }
}

///////////////////////////////////
// NEIGHBOUR - Time marching
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  sendMinAndMaxToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);

  sendDataToNeighbourBasedOnLimiterStatus(
      toRank,cellDescriptionsIndex,element,src,dest,
      false,/* isRecomputation */
      x,level);
}

void exahype::solvers::LimitingADERDGSolver::sendMinAndMaxToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!
                                                                          // |src|dest| : 1; |dest|src| : 0

  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMin()));
  assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMax()));

  // We append all the max values to the min values.
  // And then append the limiter status as double
  const int numberOfObservables = _solver->getDMPObservables();
  std::vector<double> minAndMaxToSend(2*numberOfObservables);
  for (int i=0; i<numberOfObservables; i++) {
    minAndMaxToSend[i]                     = DataHeap::getInstance().getData( solverPatch.getSolutionMin() )[faceIndex*numberOfObservables+i];
    minAndMaxToSend[i+numberOfObservables] = DataHeap::getInstance().getData( solverPatch.getSolutionMax() )[faceIndex*numberOfObservables+i];
  }

  DataHeap::getInstance().sendData(
      minAndMaxToSend, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}


void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbourBasedOnLimiterStatus(
        const int                                    toRank,
        const int                                    cellDescriptionsIndex,
        const int                                    element,
        const tarch::la::Vector<DIMENSIONS, int>&    src,
        const tarch::la::Vector<DIMENSIONS, int>&    dest,
        const bool                                   isRecomputation,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus limiterStatus = solverPatch.getLimiterStatus();

  logDebug("sendDataToNeighbourBasedOnLimiterStatus(...)", "send data for solver " << _identifier << " from rank " <<
               toRank << " at vertex x=" << x << ", level=" << level <<
               ", src=" << src << ", dest=" << dest <<", limiterStatus="<<limiterStatus);

  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok: {
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
      _limiter->sendEmptyDataToNeighbour(toRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
      } break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4: {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
      } break;
    case SolverPatch::LimiterStatus::Troubled: {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
      } break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:{
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
      } break;
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerNeighbourCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);

  _solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
  _limiter->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const exahype::MetadataHeap::HeapEntries&    neighbourMetadata,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    double**                                     tempFaceUnknowns,
    double**                                     tempStateSizedVectors,
    double**                                     tempStateSizedSquareMatrices,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  mergeWithNeighbourDataBasedOnLimiterStatus(
      fromRank,neighbourMetadata,cellDescriptionsIndex,element,src,dest,
      false,/*isRecomputation*/
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,x,level);

  mergeWithNeighbourMinAndMax(fromRank,cellDescriptionsIndex,element,src,dest,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourDataBasedOnLimiterStatus(
    const int                                    fromRank,
    const exahype::MetadataHeap::HeapEntries&    neighbourMetadata,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const bool                                   isRecomputation,
    double**                                     tempFaceUnknowns,
    double**                                     tempStateSizedVectors,
    double**                                     tempStateSizedSquareMatrices,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus limiterStatus = solverPatch.getLimiterStatus();

  logDebug("mergeWithNeighbourDataBasedOnLimiterStatus(...)", "receive data for solver " << _identifier << " from rank " <<
          fromRank << " at vertex x=" << x << ", level=" << level <<
          ", src=" << src << ", dest=" << dest << ",limiterStatus=" << SolverPatch::toString(limiterStatus));

  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4: {
      _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
      if (!isRecomputation) {
        _solver->mergeWithNeighbourData(
                 fromRank,neighbourMetadata,cellDescriptionsIndex,element,
                 src,dest,tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,x,level);
      } else {
        _solver->dropNeighbourData(fromRank,src,dest,x,level);
      }
      }break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:{
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->mergeWithNeighbourData(
          fromRank,neighbourMetadata,cellDescriptionsIndex,limiterElement,
          src,dest,tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,x,level);
      _solver->dropNeighbourData(fromRank,src,dest,x,level);
      } break;
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourMinAndMax(
    const int                                    fromRank,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) > dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

  const int numberOfObservables = _solver->getDMPObservables();
  const int receivedMinMaxIndex = DataHeap::getInstance().createData(2*numberOfObservables, 2*numberOfObservables);
  assertion(DataHeap::getInstance().getData(receivedMinMaxIndex).size()==static_cast<unsigned int>(2*numberOfObservables));
  double* receivedMinAndMax = DataHeap::getInstance().getData(receivedMinMaxIndex).data();

  DataHeap::getInstance().receiveData(receivedMinAndMax, 2*numberOfObservables, fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
  mergeSolutionMinMaxOnFace(solverPatch,faceIndex,receivedMinAndMax,receivedMinAndMax+numberOfObservables);

  DataHeap::getInstance().deleteData(receivedMinMaxIndex,true);
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch&  SolverPatch,
  const int           faceIndex,
  const double* const min, const double* const max) const {
  if (SolverPatch.getType() == SolverPatch::Cell ||
      SolverPatch.getType() == SolverPatch::Ancestor ||
      SolverPatch.getType() == SolverPatch::Descendant
      ) {
    double* solutionMin = DataHeap::getInstance().getData( SolverPatch.getSolutionMin()  ).data();
    double* solutionMax = DataHeap::getInstance().getData( SolverPatch.getSolutionMax()  ).data();

    const int numberOfObservables = _solver->getDMPObservables();
    for (int i=0; i<numberOfObservables; i++) {
      solutionMin[i+faceIndex*numberOfObservables]  = std::min( solutionMin[i+faceIndex*numberOfObservables], min[i] );
      solutionMax[i+faceIndex*numberOfObservables]  = std::max( solutionMax[i+faceIndex*numberOfObservables], max[i] );
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  for(int receives=0; receives<DataMessagesPerNeighbourCommunication; ++receives)
      DataHeap::getInstance().receiveData(
          fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);

  _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
  _solver->dropNeighbourData(fromRank,src,dest,x,level);
}

///////////////////////////////////////////////////////////
// NEIGHBOUR - Solution recomputation
///////////////////////////////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendEmptySolverAndLimiterDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
  _limiter->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourSolverAndLimiterData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  logDebug("dropNeighbourSolverAndLimiterData(...)", "drop data for solver " << _identifier << " from rank " <<
            fromRank << " at vertex x=" << x << ", level=" << level <<
            ", src=" << src << ", dest=" << dest);

  _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
  _solver->dropNeighbourData(fromRank,src,dest,x,level);
}

///////////////////////////////////////
// FORK OR JOIN
///////////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendDataToWorkerOrMasterDueToForkOrJoin(
      toRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->sendDataToWorkerOrMasterDueToForkOrJoin(
          toRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
          toRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
        toRank,x,level);
  _limiter->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
        toRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->mergeWithWorkerOrMasterDataDueToForkOrJoin(
      fromRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithWorkerOrMasterDataDueToForkOrJoin(
        fromRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->dropWorkerOrMasterDataDueToForkOrJoin(
        fromRank,x,level);
  } // !!! Receive order must be the same in master<->worker comm.
}

void exahype::solvers::LimitingADERDGSolver::dropWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->dropWorkerOrMasterDataDueToForkOrJoin(
          fromRank,x,level);
  _limiter->dropWorkerOrMasterDataDueToForkOrJoin(
          fromRank,x,level);
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->sendDataToMaster(masterRank,x,level);
//  _limiter->sendDataToMaster(masterRank,x,level);

  // Send the information to master if limiter status has changed or not
  std::vector<double> dataToSend(0,1);
  dataToSend.push_back(_limiterDomainChange ? 1.0 : -1.0); // TODO(Dominic): ugly
  assertion1(dataToSend.size()==1,dataToSend.size());
  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToMaster(...)","Sending data to master:" <<
             " data[0]=" << dataToSend[0]);
  }
  DataHeap::getInstance().sendData(
      dataToSend.data(), dataToSend.size(),
      masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithWorkerData(workerRank,x,level); // !!! Receive order must be the same in master<->worker comm.
//  _limiter->mergeWithWorkerData(workerRank,x,level); // TODO(Dominic): Revision

  // Receive the information if limiter status has changed
  std::vector<double> receivedData(1); // !!! Creates and fills the vector
  DataHeap::getInstance().receiveData(
      receivedData.data(),receivedData.size(),workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
  assertion(tarch::la::equals(receivedData[0],1.0) ||
            tarch::la::equals(receivedData[0],-1.0)); // TODO(Dominic): ugly

  bool workerLimiterDomainHasChanged = tarch::la::equals(receivedData[0],1.0) ? true : false; // TODO(Dominic): Fix this
  updateNextLimiterDomainChange(workerLimiterDomainHasChanged); // !!! It is important that we merge with the "next" field here

  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","Received data from worker:" <<
            " data[0]=" << receivedData[0]);
    logDebug("mergeWithWorkerData(...)","nextLimiterDomainChange=" << _nextLimiterDomainChange);
  }
}

bool exahype::solvers::LimitingADERDGSolver::hasToSendDataToMaster(
      const int cellDescriptionsIndex,
      const int element) {
  #if defined(Asserts) || defined(Debug)
  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    assertion(_limiter->hasToSendDataToMaster(cellDescriptionsIndex,limiterElement));
  }
  #endif

  return _solver->hasToSendDataToMaster(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::sendDataToMaster(
    const int                                     masterRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendDataToMaster(
      masterRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->sendDataToMaster(
        masterRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->sendEmptyDataToMaster(masterRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToMaster(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendEmptyDataToMaster(masterRank,x,level);
  _limiter->sendEmptyDataToMaster(masterRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const exahype::MetadataHeap::HeapEntries&    workerMetadata,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithWorkerData(
      workerRank,workerMetadata,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithWorkerData( // !!! Receive order must be the same in master<->worker comm.
        workerRank,workerMetadata,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->dropWorkerData(workerRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::dropWorkerData(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->dropWorkerData(workerRank,x,level);// !!! Receive order must be the same in master<->worker comm.
  _limiter->dropWorkerData(workerRank,x,level);
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToWorker(
    const                                        int workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->sendDataToWorker(workerRank,x,level);
//  _limiter->sendDataToWorker(workerRank,x,level); // TODO(Dominic): Revision

  // TODO(Dominic): Add information that limiter status has been
  // changed for this solver.
  // Send the information to master if limiter status has changed or not
  std::vector<double> dataToSend(0,1);
  dataToSend.push_back(_limiterDomainChange ? 1.0 : -1.0);
  assertion1(dataToSend.size()==1,dataToSend.size());
  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToWorker(...)","Sending data to worker:" <<
            " data[0]=" << dataToSend[0]);
  }
  DataHeap::getInstance().sendData(
      dataToSend.data(), dataToSend.size(),
      workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithMasterData(masterRank,x,level); // !!! Receive order must be the same in master<->worker comm.
//  _limiter->mergeWithMasterData(masterRank,x,level); // TODO(Dominic): Revision

  // Receive the information if limiter status has changed
  std::vector<double> receivedData(1); // !!! Creates and fills the vector
  DataHeap::getInstance().receiveData(
      receivedData.data(),receivedData.size(),masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
  assertion(tarch::la::equals(receivedData[0],1.0) ||
            tarch::la::equals(receivedData[0],-1.0)); // TODO(Dominic): ugly

  bool masterLimiterDomainHasChanged = tarch::la::equals(receivedData[0],1.0) ? true : false;
  _limiterDomainChange = masterLimiterDomainHasChanged; // !!! It is important that we merge with the "next" field here

  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithMasterData(...)","Received data from worker:" <<
            " data[0]=" << receivedData[0]);
    logDebug("mergeWithMasterData(...)","_limiterDomainHasChanged=" << _limiterDomainChange);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendDataToWorker(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendDataToWorker(
      workerRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->sendDataToWorker(
        workerRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->sendEmptyDataToWorker(workerRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToWorker(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendEmptyDataToWorker(workerRank,x,level);
  _limiter->sendEmptyDataToWorker(workerRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterData(
    const int                                     masterRank,
    const exahype::MetadataHeap::HeapEntries&     masterMetadata,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->mergeWithMasterData(
      masterRank,masterMetadata,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithMasterData( // !!! Receive order must be the same in master<->worker comm.
        masterRank,masterMetadata,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->dropMasterData(masterRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::dropMasterData(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->dropMasterData(masterRank,x,level);
  _limiter->dropMasterData(masterRank,x,level);
}
#endif

std::string exahype::solvers::LimitingADERDGSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::LimitingADERDGSolver::toString (std::ostream& out) const {
  out << getIdentifier() << "{_ADERDG: ";
  out << _solver->toString() << "}\n";
  out << getIdentifier() << "{_FV: ";
  out << _limiter->toString() << "}";
}
