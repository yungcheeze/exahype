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

// *.cpp
#include "kernels/limiter/generic/c/Limiter.h"

namespace exahype {
namespace solvers {

class LimitedADERDGSolver : public Solver {};

} /* namespace solvers */
} /* namespace exahype */

class exahype::solvers::LimitedADERDGSolver : public exahype::solvers::Solver {
private:
  /**
   * Log device.
   */
  static tarch::logging::Log _log;

  /**
   * The ADERDG solver.
   */
  exahype::solvers::ADERDGSolver* _ADERDGSolver;

  /**
   * The finite volumes solver used for the a posteriori subcell limiting.
   */
  exahype::solvers::FiniteVolumesSolver* _finiteVolumesSolver;

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
    return _ADERDGSolver->isValidCellDescriptionIndex(cellDescriptionsIndex);
  }

  int tryGetElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const override {
    return _ADERDGSolver->tryGetElement(cellDescriptionsIndex,solverNumber);
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
    return _ADERDGSolver->enterCell(
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
    return _ADERDGSolver->leaveCell(
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
    exahype::solvers::ADERDGSolver::CellDescription& ADERDGCellDescription =
        _ADERDGSolver->getCellDescription(cellDescriptionsIndex,element);
    exahype::solvers::FiniteVolumesSolver::CellDescription& finiteVolumesCellDescription;
    double admissibleTimeStepSize = std::numeric_limits<double>::max();

    switch (ADERDGCellDescription.getLimiterStatus(0)) {
      case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Ok:
        admissibleTimeStepSize =
            _ADERDGSolver->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
        break;
      case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
        admissibleTimeStepSize =
                    _ADERDGSolver->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
        int finiteVolumesElement =
            _finiteVolumesSolver->tryGetElement(cellDescriptionsIndex,ADERDGCellDescription.getSolverNumber());
        assertion(finiteVolumesElement!=exahype::solvers::Solver::NotFound);
        finiteVolumesCellDescription = _finiteVolumesSolver->getCellDescription(cellDescriptionsIndex,element);
        finiteVolumesCellDescription.setTimeStamp(ADERDGCellDescription.getCorrectorTimeStamp());
        finiteVolumesCellDescription.setTimeStepSize(ADERDGCellDescription.getCorrectorTimeStepSize());
        break;
      case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Troubled:
      case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsTroubledCell:
        admissibleTimeStepSize =
            _finiteVolumesSolver->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
        int finiteVolumesElement =
            _finiteVolumesSolver->tryGetElement(cellDescriptionsIndex,ADERDGCellDescription.getSolverNumber());
        assertion(finiteVolumesElement!=exahype::solvers::Solver::NotFound);
        finiteVolumesCellDescription = _finiteVolumesSolver->getCellDescription(cellDescriptionsIndex,element);
        ADERDGCellDescription.setCorrectorTimeStamp(finiteVolumesCellDescription.getTimeStamp());
        ADERDGCellDescription.setCorrectorTimeStepSize(finiteVolumesCellDescription.getTimeStepSize());
        ADERDGCellDescription.setPredictorTimeStamp(finiteVolumesCellDescription.getTimeStamp());
        ADERDGCellDescription.setPredictorTimeStepSize(finiteVolumesCellDescription.getTimeStepSize());
        break;
    }

    return admissibleTimeStepSize;
  }

  void setInitialConditions(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override {
    _ADERDGSolver->setInitialConditions(
        cellDescriptionsIndex,element,
        fineGridVertices,fineGridVerticesEnumerator);
  }

  void updateSolution(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override {
    exahype::solvers::ADERDGSolver::CellDescription& ADERDGCellDescription =
        _ADERDGSolver->getCellDescription(cellDescriptionsIndex,element);
    exahype::solvers::FiniteVolumesSolver::CellDescription& finiteVolumesCellDescription;

    // 1. Perform solution update of valid (master) solver.
    switch (ADERDGCellDescription.getLimiterStatus(0)) {
      case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Ok:
            _ADERDGSolver->updateSolution(
                cellDescriptionsIndex,element,
                fineGridVertices,fineGridVerticesEnumerator);
        break;
      case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
        _ADERDGSolver->updateSolution(
            cellDescriptionsIndex,element,
            fineGridVertices,fineGridVerticesEnumerator);
        int finiteVolumesElement =
            _finiteVolumesSolver->tryGetElement(cellDescriptionsIndex,ADERDGCellDescription.getSolverNumber());
        assertion(finiteVolumesElement!=exahype::solvers::Solver::NotFound);
        finiteVolumesCellDescription = _finiteVolumesSolver->getCellDescription(cellDescriptionsIndex,element);

        double* ADERDGSolution        = DataHeap::getInstance().getData(
            ADERDGCellDescription.getSolution()).data();
        double* finiteVolumesSolution = DataHeap::getInstance().getData(
            finiteVolumesCellDescription.getSolution()).data();

        // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
        kernels::limiter::generic::c::projectOnFVLimiterSpace(ADERDGSolution,_ADERDGSolver->getNumberOfVariables(),
            _ADERDGSolver->getNodesPerCoordinateAxis(),finiteVolumesSolution);
        break;
      case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Troubled:
      case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsTroubledCell:
        _finiteVolumesSolver->updateSolution(
                        cellDescriptionsIndex,element,
                        fineGridVertices,fineGridVerticesEnumerator);
        int finiteVolumesElement =
            _finiteVolumesSolver->tryGetElement(cellDescriptionsIndex,ADERDGCellDescription.getSolverNumber());
        assertion(finiteVolumesElement!=exahype::solvers::Solver::NotFound);
        finiteVolumesCellDescription = _finiteVolumesSolver->getCellDescription(cellDescriptionsIndex,element);

        // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
        kernels::limiter::generic::c::projectOnDGSpace(finiteVolumesSolution,_ADERDGSolver->getNumberOfVariables(),
                    _ADERDGSolver->getNodesPerCoordinateAxis(),ADERDGSolution);
        break;
    }

    // 1. Compute limiter status after update.
    switch (ADERDGCellDescription.getLimiterStatus(0)) {

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
