/*
 * LimitedADERDGSolver.h
 *
 *  Created on: 26 Oct 2016
 *      Author: dominic
 */

#ifndef LIMITEDADERDGSOLVER_H_
#define LIMITEDADERDGSOLVER_H_

#include "exahype/solvers/Solver.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "tarch/multicore/BooleanSemaphore.h"

// *.cpp
#include "tarch/multicore/Lock.h"

#include "kernels/limiter/generic/c/Limiter.h"

namespace exahype {
namespace solvers {

class LimitedADERDGSolver : public Solver {};

} /* namespace solvers */
} /* namespace exahype */

class exahype::solvers::LimitedADERDGSolver : public exahype::solvers::Solver {
private:
  /**
   * A semaphore for serialising the adding and removing of limiter patches
   * for the limiter heap.
   */
  tarch::multicore::BooleanSemaphore _semaphoreForLimiterHeapAccess;

  /**
   * A flag indicating that the limiter domain has changed.
   * This might be the case if either a cell has been
   * newly marked as troubled or healed.
   */
  bool _limiterDomainChanged;

  typedef exahype::records::ADERDGCellDescription SolverPatch;
  typedef peano::heap::PlainHeap<SolverPatch> SolverHeap;

  typedef exahype::records::FiniteVolumesCellDescription LimiterPatch;
  typedef peano::heap::PlainHeap<LimiterPatch> LimiterHeap;


  /**
   * Log device.
   */
  static tarch::logging::Log _log;

  /**
   * The ADERDG solver.
   */
  exahype::solvers::ADERDGSolver* _solver;

  /**
   * The finite volumes solver used for the a posteriori subcell limiting.
   */
  exahype::solvers::FiniteVolumesSolver* _limiter;


  void determineLimiterStatusOfCellDescription(SolverPatch& solverPatch) {
    SolverPatch::LimiterStatus limiterStatus = SolverPatch::LimiterStatus::Ok;
    for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
      // 1. Change the newly communicated statuses to "old" ones.
      switch (limiterStatus) {
        case SolverPatch::LimiterStatus::NewlyNeighbourIsTroubledCell:
          solverPatch.setLimiterStatus(i,SolverPatch::LimiterStatus::NeighbourIsTroubledCell);
          break;
        case SolverPatch::LimiterStatus::NewlyNeighbourIsNeighbourOfTroubledCell:
          solverPatch.setLimiterStatus(i,SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell);
          break;
        default:
          break;
      }

      // 2. Determine new limiter status.
      switch (limiterStatus) {
        case SolverPatch::LimiterStatus::Ok:
          switch (solverPatch.getLimiterStatus(i)) {
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
          switch (solverPatch.getLimiterStatus(i)) {
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

    // 3. Finally, set the limiter status on all faces to the determined value.
    for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
      solverPatch.setLimiterStatus(i,limiterStatus);
    }
  }
public:
  /*
   * A time stamp minimised over all the ADERDG and FV solver
   * patches.
   */
  double getMinTimeStamp() const override;

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  double getMinTimeStepSize() const override;

  void updateNextTimeStepSize( double value ) override;

  void initInitialTimeStamp(double value) override;

  void synchroniseTimeStepping(
          const int cellDescriptionsIndex,
          const int element) override;

  void startNewTimeStep() override;

  void reinitTimeStepData() override;

  double getNextMinTimeStepSize() const override;

  bool isValidCellDescriptionIndex(
      const int cellDescriptionsIndex) const override {
    return _solver->isValidCellDescriptionIndex(cellDescriptionsIndex);
  }

  int tryGetElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const override {
    return _solver->tryGetElement(cellDescriptionsIndex,solverNumber);
  }

  ///////////////////////////////////
  // MODIFY CELL DESCRIPTION
  ///////////////////////////////////
  bool enterCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) override {
    return _solver->enterCell(
        fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
        coarseGridVertices,coarseGridVerticesEnumerator,coarseGridCell,
        fineGridPositionOfCell,solverNumber);
  }

  bool leaveCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) override {
    return _solver->leaveCell(
            fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
            coarseGridVertices,coarseGridVerticesEnumerator,coarseGridCell,
            fineGridPositionOfCell,solverNumber);
  }

  ///////////////////////////////////
  // CELL-LOCAL
  //////////////////////////////////
  double startNewTimeStep(
      const int cellDescriptionsIndex,
      const int element,
      double*   tempEigenvalues) override {
    exahype::solvers::ADERDGSolver::CellDescription& solverPatch =
        _solver->getCellDescription(cellDescriptionsIndex,element);
    exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch;
    double admissibleTimeStepSize = std::numeric_limits<double>::max();

    switch (solverPatch.getLimiterStatus(0)) {
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

  void setInitialConditions(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override {
    _solver->setInitialConditions(
        cellDescriptionsIndex,element,
        fineGridVertices,fineGridVerticesEnumerator);
  }

  /**
   * Checks if updated solution
   * of the ADER-DG solver is Troubled,
   * or Ok.
   */
  bool solutionIsTroubled(SolverPatch& solverPatch) {
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
   *
   */
  bool determineLimiterStatus(
      const int cellDescriptionsIndex,
      const int element) {
    SolverPatch& solverPatch =
        _solver->getCellDescription(cellDescriptionsIndex,element);
    switch (solverPatch.getLimiterStatus(0)) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      if (solutionIsTroubled(solverPatch)) {
        for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
          solverPatch.setLimiterStatus(i,SolverPatch::LimiterStatus::Troubled);
        }
        return true;
      }
      break;
    case SolverPatch::LimiterStatus::Troubled:
      if (!solutionIsTroubled(solverPatch)) {
        for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
          solverPatch.setLimiterStatus(i,SolverPatch::LimiterStatus::Ok);
        }
        return true;
      }
      break;
    default:
      break;
    }

    return false;
  }

  /**
   * This method assumes the ADERDG solver's cell-local limiter status has
   * already been determined.
   */
  void updateSolution(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override {
    SolverPatch& solverPatch =
        _solver->getCellDescription(cellDescriptionsIndex,element);
    exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch;
    int finiteVolumesElement;

    double* solution    = nullptr;
    double* update      = nullptr;
    double* solutionMin = nullptr;
    double* solutionMax = nullptr;

    double* limiterSolution = nullptr;

    // 1. Update the solution in the cells
    switch (solverPatch.getLimiterStatus(0)) {
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

      solution = DataHeap::getInstance().getData(
          solverPatch.getSolution()).data();
      limiterSolution = DataHeap::getInstance().getData(
          limiterPatch.getSolution()).data();

      // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      kernels::limiter::generic::c::projectOnFVLimiterSpace(
          solution,_solver->getNumberOfVariables(),
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
          _solver->getNodesPerCoordinateAxis(),solution);
      break;
    }
  }

  void preProcess(
      const int cellDescriptionsIndex,
      const int element) override;

  void postProcess(
      const int cellDescriptionsIndex,
      const int element) override;

  void prolongateDataAndPrepareDataRestriction(
      const int cellDescriptionsIndex,
      const int element) override;

  void restrictData(
      const int cellDescriptionsIndex,
      const int element,
      const int parentCellDescriptionsIndex,
      const int parentElement,
      const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) override;


  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  void mergeNeighbours(
      const int                                 cellDescriptionsIndex1,
      const int                                 element1,
      const int                                 cellDescriptionsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2,
      double**                                  tempFaceUnknownsArrays,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) override;

  void mergeWithBoundaryData(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      double**                                  tempFaceUnknownsArrays,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) override;


#ifdef Parallel
  static void sendCellDescriptions(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const peano::heap::MessageType&               messageType,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level);

  static void sendEmptyCellDescriptions(
      const int                                     toRank,
      const peano::heap::MessageType&               messageType,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level);

  /**
   * Receives cell descriptions from rank \p fromRank
   * and resets the data heap indices to -1.
   *
   * If a received cell description has the same
   * solver number as a cell description in the
   * array at address \p cellDescriptionsIndex,
   * we merge the metadata (time stamps, time step size)
   * of both cell descriptions.
   *
   * If no cell description in the array at address
   * \p cellDescriptionsIndex can be found with the
   * same solver number than a received cell description,
   * we push the received cell description to
   * the back of the array at address \p cellDescriptions
   * Index.
   *
   * This operation is intended to be used in combination
   * with the solver method mergeWithWorkerOrMasterDataDueToForkOrJoin(...).
   * Here, we would merge first the cell descriptions sent by the master and worker
   * and then merge the data that is sent out right after.
   */
  static void mergeCellDescriptionsWithRemoteData(
      const int                                     fromRank,
      const int                                     cellDescriptionsIndex,
      const peano::heap::MessageType&               messageType,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level);

  /**
   * Drop cell descriptions received from \p fromRank.
   */
  static void dropCellDescriptions(
      const int                                     fromRank,
      const peano::heap::MessageType&               messageType,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level);

  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  void mergeWithNeighbourMetadata(
        const int neighbourTypeAsInt,
        const int cellDescriptionsIndex,
        const int element) override {
    assertionMsg(false,"Please implement!");
  }

  void sendDataToNeighbour(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     elementIndex,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false,"Please implement!");
  }

  void sendEmptyDataToNeighbour(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void mergeWithNeighbourData(
      const int                                    fromRank,
      const int                                    neighbourTypeAsInt,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      double**                                     tempFaceUnknownsArrays,
      double**                                     tempStateSizedVectors,
      double**                                     tempStateSizedSquareMatrices,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override {
    assertionMsg(false,"Please implement!");
  }


  void dropNeighbourData(
      const int                                     fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void sendDataToWorkerOrMasterDueToForkOrJoin(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void mergeWithWorkerOrMasterDataDueToForkOrJoin(
      const int                                     fromRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void dropWorkerOrMasterDataDueToForkOrJoin(
      const int                                     fromRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  ///////////////////////////////////
  // WORKER->MASTER
  ///////////////////////////////////
  bool hasToSendDataToMaster(
        const int cellDescriptionsIndex,
        const int element) override;

  void sendDataToMaster(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void mergeWithWorkerData(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void sendDataToMaster(
      const int                                     masterRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void sendEmptyDataToMaster(
      const int                                     masterRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void mergeWithWorkerData(
      const int                                    workerRank,
      const int                                    workerTypeAsInt,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void dropWorkerData(
      const int                                     workerRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  ///////////////////////////////////
  // MASTER->WORKER
  ///////////////////////////////////
  void sendDataToWorker(
      const                                        int workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void mergeWithMasterData(
      const                                        int masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void sendDataToWorker(
      const int                                     workerRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void sendEmptyDataToWorker(
      const int                                     workerRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void mergeWithMasterData(
      const int                                     masterRank,
      const int                                     masterTypeAsInt,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void dropMasterData(
      const int                                     masterRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;
#endif

  std::string toString() const override;

  void toString (std::ostream& out) const override;
};


#endif /* LIMITEDADERDGSOLVER_H_ */
