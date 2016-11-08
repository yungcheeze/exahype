/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include "exahype/solvers/FiniteVolumesSolver.h"

#include <string>
#include <limits>

#include "exahype/Cell.h"
#include "exahype/Vertex.h"
#include "exahype/VertexOperations.h"

#include "peano/utils/Loop.h"

#include "tarch/multicore/Lock.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"


namespace {
constexpr const char* tags[]{"solutionUpdate", "stableTimeStepSize"};
}  // namespace

tarch::logging::Log exahype::solvers::FiniteVolumesSolver::_log( "exahype::solvers::FiniteVolumesSolver");

exahype::solvers::FiniteVolumesSolver::FiniteVolumesSolver(
    const std::string& identifier, int numberOfVariables,
    int numberOfParameters, int nodesPerCoordinateAxis, int ghostLayerWidth,
    double maximumMeshSize,
    exahype::solvers::Solver::TimeStepping timeStepping,
    std::unique_ptr<profilers::Profiler> profiler)
    : Solver(identifier, exahype::solvers::Solver::Type::FiniteVolumes,
             numberOfVariables, numberOfParameters, nodesPerCoordinateAxis,
             maximumMeshSize, timeStepping, std::move(profiler)),
      _unknownsPerPatch((numberOfVariables + numberOfParameters) *
                       power(nodesPerCoordinateAxis, DIMENSIONS + 0)),
      _ghostLayerWidth(ghostLayerWidth),
      _ghostValuesPerPatch((numberOfVariables + numberOfParameters) *
                       power(nodesPerCoordinateAxis+2*ghostLayerWidth, DIMENSIONS + 0) - _unknownsPerPatch),
      _unknownsPerPatchFace(
          (numberOfVariables + numberOfParameters)*power(nodesPerCoordinateAxis, DIMENSIONS - 1)),
      _unknownsPerPatchBoundary(
          DIMENSIONS_TIMES_TWO *_unknownsPerPatchFace),
      _minTimeStamp(std::numeric_limits<double>::max()),
      _minTimeStepSize(std::numeric_limits<double>::max()),
      _nextMinTimeStepSize(std::numeric_limits<double>::max()) {
  assertion3(_unknownsPerPatch > 0, numberOfVariables, numberOfParameters,
             nodesPerCoordinateAxis);

  // register tags with profiler
  for (const char* tag : tags) {
    _profiler->registerTag(tag);
  }
}

int exahype::solvers::FiniteVolumesSolver::getUnknownsPerPatch() const {
  return _unknownsPerPatch;
}

int exahype::solvers::FiniteVolumesSolver::getGhostLayerWidth() const {
  return _ghostLayerWidth;
}

int exahype::solvers::FiniteVolumesSolver::getGhostValuesPerPatch() const {
  return _ghostValuesPerPatch;
}

int exahype::solvers::FiniteVolumesSolver::getUnknownsPerFace() const {
  return _unknownsPerPatchFace;
}

double exahype::solvers::FiniteVolumesSolver::getMinTimeStamp() const {
  return _minTimeStamp;
}

double exahype::solvers::FiniteVolumesSolver::getMinTimeStepSize() const {
  return _minTimeStepSize;
}

void exahype::solvers::FiniteVolumesSolver::updateNextTimeStepSize(
    double value) {
  _nextMinTimeStepSize = std::min(_nextMinTimeStepSize, value);
}

void exahype::solvers::FiniteVolumesSolver::initInitialTimeStamp(double value) {
  _minTimeStamp = value;
}

void exahype::solvers::FiniteVolumesSolver::synchroniseTimeStepping(
    CellDescription& cellDescription) {
  switch (_timeStepping) {
    case TimeStepping::Global:
      cellDescription.setTimeStamp(_minTimeStamp);
      cellDescription.setTimeStepSize(_minTimeStepSize);
      break;
    case TimeStepping::GlobalFixed:
      cellDescription.setTimeStamp(_minTimeStamp);
      cellDescription.setTimeStepSize(_minTimeStepSize);
      break;
  }
}

void exahype::solvers::FiniteVolumesSolver::synchroniseTimeStepping(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  synchroniseTimeStepping(cellDescription);
}

void exahype::solvers::FiniteVolumesSolver::startNewTimeStep() {
  switch (_timeStepping) {
    case TimeStepping::Global:
      _minTimeStamp        += _minTimeStepSize;
      _minTimeStepSize      = _nextMinTimeStepSize;
      _nextMinTimeStepSize  = std::numeric_limits<double>::max();
      break;
    case TimeStepping::GlobalFixed:
      _minTimeStamp        += _minTimeStepSize;
      _minTimeStepSize      = _nextMinTimeStepSize;
      break;
  }
}

void exahype::solvers::FiniteVolumesSolver::reinitTimeStepData() {
  switch (_timeStepping) {
    case TimeStepping::Global:
      // do nothing
      break;
    case TimeStepping::GlobalFixed:
      // do nothing
      break;
  }
}

double exahype::solvers::FiniteVolumesSolver::getNextMinTimeStepSize() const {
  return _nextMinTimeStepSize;
}

int exahype::solvers::FiniteVolumesSolver::tryGetElement(
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  if (Heap::getInstance().isValidIndex(cellDescriptionsIndex)) {
    int element=0;
    for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
      if (p.getSolverNumber()==solverNumber) {
        return element;
      }
      ++element;
    }
  }
  return NotFound;
}

exahype::solvers::Solver::SubcellPosition exahype::solvers::FiniteVolumesSolver::computeSubcellPositionOfCellOrAncestor(
        const int cellDescriptionsIndex,
        const int element) {
  // TODO(Dominic): Comment code in as soon as we have all the required fields
  // on the cell description.
//  CellDescription& cellDescription =
//      getCellDescription(cellDescriptionsIndex,element);
//
//  return
//      exahype::amr::computeSubcellPositionOfCellOrAncestor
//      <CellDescription,Heap>(cellDescription);
  SubcellPosition empty;
  return empty;
}

///////////////////////////////////
// MODIFY CELL DESCRIPTION
///////////////////////////////////
bool exahype::solvers::FiniteVolumesSolver::enterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber) {
  // Fine grid cell based uniform mesh refinement.
  int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (fineGridCellElement==exahype::solvers::Solver::NotFound &&
      tarch::la::allSmallerEquals(fineGridVerticesEnumerator.getCellSize(),getMaximumMeshSize()) &&
      tarch::la::allGreater(coarseGridVerticesEnumerator.getCellSize(),getMaximumMeshSize())) {
    addNewCell(fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
               multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex,
               solverNumber);
    // Fine grid cell based adaptive mesh refinement operations.
  }

  return false;
}

void exahype::solvers::FiniteVolumesSolver::addNewCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const int coarseGridCellDescriptionsIndex,
    const int solverNumber) {
  fineGridCell.addNewCellDescription(
              solverNumber,
              CellDescription::Cell,
//              CellDescription::None,
              fineGridVerticesEnumerator.getLevel(),
              coarseGridCellDescriptionsIndex,
              fineGridVerticesEnumerator.getCellSize(),
              fineGridVerticesEnumerator.getVertexPosition());
  int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  CellDescription& fineGridCellDescription =
      getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridCellElement);
  ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
  exahype::Cell::determineInsideAndOutsideFaces(
            fineGridCellDescription,
            fineGridVertices,
            fineGridVerticesEnumerator);
}

void exahype::solvers::FiniteVolumesSolver::ensureNoUnnecessaryMemoryIsAllocated(CellDescription& cellDescription) {
  if (DataHeap::getInstance().isValidIndex(cellDescription.getSolution())) {
    switch (cellDescription.getType()) {
      case CellDescription::Erased:
      case CellDescription::EmptyAncestor:
      case CellDescription::EmptyDescendant:
      case CellDescription::Ancestor:
      case CellDescription::Descendant:
        {
        waitUntilAllBackgroundTasksHaveTerminated();
        tarch::multicore::Lock lock(_heapSemaphore);

        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()));
        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getOldSolution()));

        DataHeap::getInstance().deleteData(cellDescription.getSolution());
        DataHeap::getInstance().deleteData(cellDescription.getOldSolution());

        cellDescription.setSolution(-1);
        cellDescription.setOldSolution(-1);
        }
        break;
      default:
        break;
    }
  }

//  if (DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor())) {
//    switch (cellDescription.getType()) {
//      case CellDescription::Erased:
//      case CellDescription::EmptyAncestor:
//      case CellDescription::EmptyDescendant:
//        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()));
//        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionMin()));
//        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionMax()));
//
//        DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictor());
//        DataHeap::getInstance().deleteData(cellDescription.getFluctuation());
//        DataHeap::getInstance().deleteData(cellDescription.getSolutionMin());
//        DataHeap::getInstance().deleteData(cellDescription.getSolutionMax());
//
//        cellDescription.setExtrapolatedPredictor(-1);
//        cellDescription.setFluctuation(-1);
//        break;
//      default:
//        break;
//    }
//  }
}

void exahype::solvers::FiniteVolumesSolver::ensureNecessaryMemoryIsAllocated(CellDescription& cellDescription) {
  switch (cellDescription.getType()) {
    case CellDescription::Cell:
      if (!DataHeap::getInstance().isValidIndex(cellDescription.getSolution())) {
        assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getSolution()));
        // Allocate volume DoF for limiter
        const int size = _unknownsPerPatch+_ghostValuesPerPatch;

        waitUntilAllBackgroundTasksHaveTerminated();
        tarch::multicore::Lock lock(_heapSemaphore);

        cellDescription.setSolution(DataHeap::getInstance().createData(size, size, DataHeap::Allocation::DoNotUseAnyRecycledEntry));
        cellDescription.setOldSolution(DataHeap::getInstance().createData(size, size, DataHeap::Allocation::DoNotUseAnyRecycledEntry));
      }
      break;
    default:
      break;
  }

//  switch (cellDescription.getType()) {
//    case CellDescription::Cell:
//    case CellDescription::Ancestor:
//    case CellDescription::Descendant:
//      if (!DataHeap::getInstance().isValidIndex(
//          cellDescription.getExtrapolatedPredictor())) {
//        assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()));
//
//        // Allocate face DoF
//        const int unknownsPerCellBoundary = getUnknownsPerCellBoundary();
//        cellDescription.setExtrapolatedPredictor(DataHeap::getInstance().createData(
//            unknownsPerCellBoundary, unknownsPerCellBoundary));
//        cellDescription.setFluctuation(DataHeap::getInstance().createData(
//            unknownsPerCellBoundary, unknownsPerCellBoundary));
//
//        // Allocate volume DoF for limiter (we need for every of the 2*DIMENSIONS faces an array of min values
//        // and array of max values of the neighbour at this face).
//        const int unknownsPerCell = _unknownsPerPatch;
//        cellDescription.setSolutionMin(DataHeap::getInstance().createData(
//            unknownsPerCell * 2 * DIMENSIONS, unknownsPerCell * 2 * DIMENSIONS));
//        cellDescription.setSolutionMax(DataHeap::getInstance().createData(
//            unknownsPerCell * 2 * DIMENSIONS, unknownsPerCell * 2 * DIMENSIONS));
//
//        // !!!
//        // TODO(Dominic): Make sure this everywhere initialised correctly.
//        // !!!
//        for (int i=0; i<unknownsPerCell * 2 * DIMENSIONS; i++) {
//          DataHeap::getInstance().getData( cellDescription.getSolutionMax() )[i] = std::numeric_limits<double>::max();
//          DataHeap::getInstance().getData( cellDescription.getSolutionMin() )[i] = -std::numeric_limits<double>::max(); // "-", min
//        }
//      }
//      break;
//    default:
//      break;
//  }
}

bool exahype::solvers::FiniteVolumesSolver::leaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber) {
//  assertionMsg(false,"Not implemented."); // TODO(Dominic): Implement.
  return false;
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////
double exahype::solvers::FiniteVolumesSolver::startNewTimeStep(
    const int cellDescriptionsIndex,
    const int element,
    double*   tempEigenvalues) {
  CellDescription& p = getCellDescription(cellDescriptionsIndex,element);

  if (p.getType()==exahype::records::FiniteVolumesCellDescription::Cell) {
    //         assertion1(p.getRefinementEvent()==exahype::records::FiniteVolumesCellDescription::None,p.toString()); // todo
    double* solution = exahype::DataHeap::getInstance().getData(p.getSolution()).data();

    double admissibleTimeStepSize = stableTimeStepSize(
        solution, tempEigenvalues, p.getSize());

    assertion(!std::isnan(admissibleTimeStepSize));

    p.setTimeStamp(p.getTimeStamp()+p.getTimeStepSize());
    p.setTimeStepSize(admissibleTimeStepSize);

    return admissibleTimeStepSize;
  }

  return std::numeric_limits<double>::max();
}

void exahype::solvers::FiniteVolumesSolver::setInitialConditions(
    const int cellDescriptionsIndex,
    const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  // reset helper variables
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  exahype::Cell::resetNeighbourMergeHelperVariables(
      cellDescription,fineGridVertices,fineGridVerticesEnumerator); // TODO(Dominic): Add flags.

  if (cellDescription.getType()==CellDescription::Cell
//      && cellDescription.getRefinementEvent()==CellDescription::None
      ) {
    double* solution = exahype::DataHeap::getInstance().getData(cellDescription.getSolution()).data();

    if (hasToAdjustSolution(
        cellDescription.getOffset()+0.5*cellDescription.getSize(),
        cellDescription.getSize(),
        cellDescription.getTimeStamp(),
        cellDescription.getTimeStepSize())) {
      solutionAdjustment(
          solution,
          cellDescription.getOffset()+0.5*cellDescription.getSize(),
          cellDescription.getSize(),
          cellDescription.getTimeStamp(),
          cellDescription.getTimeStepSize());
    }

    for (int i=0; i<_unknownsPerPatch+_ghostValuesPerPatch; i++) {
      assertion3(std::isfinite(solution[i]),cellDescription.toString(),"setInitialConditions(...)",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
  }
}

void exahype::solvers::FiniteVolumesSolver::updateSolution(
    const int cellDescriptionsIndex,
    const int element,
    double** tempStateSizedArrays,
    double** tempUnknowns,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  // reset helper variables
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  exahype::Cell::resetNeighbourMergeHelperVariables(
      cellDescription,fineGridVertices,fineGridVerticesEnumerator);

  double* solution    = DataHeap::getInstance().getData(cellDescription.getOldSolution()).data();
  double* newSolution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
  std::copy(newSolution,newSolution+_unknownsPerPatch+_ghostValuesPerPatch,solution); // Copy (current solution) in old solution field.

  dfor(i,_nodesPerCoordinateAxis+_ghostLayerWidth) {
    if (tarch::la::allSmaller(i,_nodesPerCoordinateAxis+_ghostLayerWidth)
    && tarch::la::allGreater(i,_ghostLayerWidth-1)) {
      for (int unknown=0; unknown < _numberOfVariables; unknown++) {
        int iScalar = peano::utils::dLinearisedWithoutLookup(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth)*_numberOfVariables+unknown;
        assertion3(std::isfinite(newSolution[iScalar]),cellDescription.toString(),newSolution[iScalar],i.toString());
      }
    }
  }

  double admissibleTimeStepSize=0;
  solutionUpdate(
      newSolution,solution,tempStateSizedArrays,tempUnknowns,
      cellDescription.getSize(),cellDescription.getTimeStepSize(),admissibleTimeStepSize);

  if (admissibleTimeStepSize * 1.001 < cellDescription.getTimeStepSize()) { //TODO JMG 1.001 factor to prevent same dt computation to throw logerror
    logWarning("updateSolution(...)","Finite volumes solver time step size harmed CFL condition. dt="<<cellDescription.getTimeStepSize()<<", dt_adm=" << admissibleTimeStepSize);
  }

  if (hasToAdjustSolution(
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getTimeStamp()+cellDescription.getTimeStepSize(),cellDescription.getTimeStepSize())) {
    solutionAdjustment(
        newSolution,
        cellDescription.getOffset()+0.5*cellDescription.getSize(),
        cellDescription.getSize(),
        cellDescription.getTimeStamp()+cellDescription.getTimeStepSize(),
        cellDescription.getTimeStepSize());
  }

  dfor(i,_nodesPerCoordinateAxis+_ghostLayerWidth) {
    if (tarch::la::allSmaller(i,_nodesPerCoordinateAxis+_ghostLayerWidth)
    && tarch::la::allGreater(i,_ghostLayerWidth-1)) {
      for (int unknown=0; unknown < _numberOfVariables; unknown++) {
        int iScalar = peano::utils::dLinearisedWithoutLookup(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth)*_numberOfVariables+unknown;
        assertion3(std::isfinite(newSolution[iScalar]),cellDescription.toString(),newSolution[iScalar],i.toString());
      }
    }
  }
}


void exahype::solvers::FiniteVolumesSolver::rollbackSolution(
    const int cellDescriptionsIndex,
    const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  // Simply swap the array pointers.
  const int oldSolution = cellDescription.getOldSolution();
  cellDescription.setOldSolution(cellDescription.getSolution());
  cellDescription.setSolution(oldSolution);
}

void exahype::solvers::FiniteVolumesSolver::preProcess(
        const int cellDescriptionsIndex,
        const int element) {
  // TODO(Dominic)
//  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::postProcess(
        const int cellDescriptionsIndex,
        const int element) {
  // TODO(Dominic)
//  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::prolongateDataAndPrepareDataRestriction(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Cell)
    return;

  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::restrictData(
    const int cellDescriptionsIndex,
    const int element,
    const int parentCellDescriptionsIndex,
    const int parentElement,
    const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  assertionMsg(false,"Please implement!");
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::FiniteVolumesSolver::mergeNeighbours(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    double**                                  tempFaceUnknownsArrays,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  CellDescription& cellDescription1 = getCellDescription(cellDescriptionsIndex1,element1);
  CellDescription& cellDescription2 = getCellDescription(cellDescriptionsIndex2,element2);

  if (cellDescription1.getType()==CellDescription::Cell ||
      cellDescription2.getType()==CellDescription::Cell) {
    synchroniseTimeStepping(cellDescription1);
    synchroniseTimeStepping(cellDescription2);

    double* solution1 = DataHeap::getInstance().getData(cellDescription1.getSolution()).data();
    double* solution2 = DataHeap::getInstance().getData(cellDescription2.getSolution()).data();

    ghostLayerFilling(solution1,solution2,pos1-pos2);
    ghostLayerFilling(solution2,solution1,pos2-pos1);
  }

  return;

  assertionMsg(false,"Not implemented.");
}

void exahype::solvers::FiniteVolumesSolver::mergeWithBoundaryData(
    const int                                 cellDescriptionsIndex,
    const int                                 element,
    const tarch::la::Vector<DIMENSIONS, int>& posCell,
    const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
    double**                                  tempFaceUnknownsArrays,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Cell) {
    synchroniseTimeStepping(cellDescription);

    double* luh       = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
    double* luhbndIn  = tempFaceUnknownsArrays[0];
    double* luhbndOut = tempFaceUnknownsArrays[1];
    assertion2(tarch::la::countEqualEntries(posCell,posBoundary)==DIMENSIONS-1,posCell.toString(),posBoundary.toString());

    boundaryLayerExtraction(luhbndIn,luh,posBoundary-posCell);

    const int normalNonZero = tarch::la::equalsReturnIndex(posCell, posBoundary);
    assertion(normalNonZero >= 0 && normalNonZero < DIMENSIONS);
    const int faceIndex = 2 * normalNonZero +
        (posCell(normalNonZero) < posBoundary(normalNonZero) ? 1 : 0);

    boundaryConditions(
        luhbndOut,luhbndIn,
        cellDescription.getOffset()+0.5*cellDescription.getSize(),
        cellDescription.getSize(),
        cellDescription.getTimeStamp(),
        cellDescription.getTimeStepSize(),
        faceIndex,
        normalNonZero);

    ghostLayerFillingAtBoundary(luh,luhbndOut,posBoundary-posCell);
  }
}


#ifdef Parallel
const int exahype::solvers::FiniteVolumesSolver::DataMessagesPerNeighbourCommunication    = 1;
const int exahype::solvers::FiniteVolumesSolver::DataMessagesPerForkOrJoinCommunication   = 1;
const int exahype::solvers::FiniteVolumesSolver::DataMessagesPerMasterWorkerCommunication = 1;

void exahype::solvers::FiniteVolumesSolver::sendCellDescriptions(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::sendEmptyCellDescriptions(
    const int                                     toRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::HeapEntries emptyMessage(0);
  Heap::getInstance().sendData(emptyMessage,
      toRank,x,level,messageType);
}

void exahype::solvers::FiniteVolumesSolver::mergeCellDescriptionsWithRemoteData(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  waitUntilAllBackgroundTasksHaveTerminated();
  tarch::multicore::Lock lock(_heapSemaphore);


  assertionMsg(false,"Please implement!");
}

/**
 * Drop cell descriptions received from \p fromRank.
 */
void exahype::solvers::FiniteVolumesSolver::dropCellDescriptions(
    const int                                     fromRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::getInstance().receiveData(fromRank,x,level,messageType);
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////

void exahype::solvers::FiniteVolumesSolver::sendEmptyDataToNeighbour(
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
}

void exahype::solvers::FiniteVolumesSolver::dropNeighbourData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  for(int receives=0; receives<DataMessagesPerNeighbourCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
}

///////////////////////////////////
// FORK OR JOIN
///////////////////////////////////

void exahype::solvers::FiniteVolumesSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertionMsg(false,"Please implement!");
}


void exahype::solvers::FiniteVolumesSolver::sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerForkOrJoinCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, toRank, x, level,
        peano::heap::MessageType::ForkOrJoinCommunication);
}


void exahype::solvers::FiniteVolumesSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::dropWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  for(int receives=0; receives<DataMessagesPerForkOrJoinCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        fromRank, x, level,
        peano::heap::MessageType::ForkOrJoinCommunication);
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////
void exahype::solvers::FiniteVolumesSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  assertionMsg(false,"Please implement!");
}


void exahype::solvers::FiniteVolumesSolver::sendDataToMaster(
    const int                                     masterRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::sendEmptyDataToMaster(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerMasterWorkerCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::FiniteVolumesSolver::mergeWithWorkerData(
    const int                                     workerRank,
    const int                                     workerTypeAsInt,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::dropWorkerData(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  for(int receives=0; receives<DataMessagesPerMasterWorkerCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
void exahype::solvers::FiniteVolumesSolver::sendDataToWorker(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  assertionMsg(false,"Please implement!");

  // todo send time step size
}

bool exahype::solvers::FiniteVolumesSolver::hasToSendDataToMaster(
    const int cellDescriptionsIndex,
    const int element) {
  assertionMsg(false,"Please implement!");
  return false;
}

void exahype::solvers::FiniteVolumesSolver::sendDataToWorker(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::sendEmptyDataToWorker(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerMasterWorkerCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::FiniteVolumesSolver::mergeWithMasterData(
    const int                                     masterRank,
    const int                                     masterTypeAsInt,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::dropMasterData(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
        const int                                     level) {
  for(int receives=0; receives<DataMessagesPerMasterWorkerCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
}
#endif

std::string exahype::solvers::FiniteVolumesSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::FiniteVolumesSolver::toString (std::ostream& out) const {
  out << "(";
  out << "_identifier:" << _identifier;
  out << ",";
  out << "_type:" << exahype::solvers::Solver::toString(_type);
  out << ",";
  out << "_numberOfVariables:" << _numberOfVariables;
  out << ",";
  out << "_numberOfParameters:" << _numberOfParameters;
  out << ",";
  out << "_nodesPerCoordinateAxis:" << _nodesPerCoordinateAxis;
  out << ",";
  out << "_maximumMeshSize:" << _maximumMeshSize;
  out << ",";
  out << "_timeStepping:" << exahype::solvers::Solver::toString(_timeStepping); // only solver attributes
  out << ",";
  out << "_unknownsPerPatchFace:" << _unknownsPerPatchFace;
  out << ",";
  out << "_unknownsPerPatchBoundary:" << _unknownsPerPatchBoundary;
  out << ",";
  out << "_unknownsPerPatch:" << _unknownsPerPatch;
  out << ",";
  out << "_minTimeStamp:" << _minTimeStamp;
  out << ",";
  out << "_minTimeStepSize:" << _minTimeStepSize;
  out << ",";
  out << "_nextMinTimeStepSize:" << _nextMinTimeStepSize;
  out <<  ")";
}
