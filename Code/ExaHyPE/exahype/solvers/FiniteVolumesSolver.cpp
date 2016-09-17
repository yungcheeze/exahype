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

#include <limits>

#include "exahype/Cell.h"
#include "exahype/Vertex.h"
#include "exahype/VertexOperations.h"

namespace {
constexpr const char* tags[]{"solutionUpdate", "stableTimeStepSize"};
}  // namespace

tarch::logging::Log exahype::solvers::FiniteVolumesSolver::_log( "exahype::solvers::FiniteVolumesSolver");

exahype::solvers::FiniteVolumesSolver::FiniteVolumesSolver(
    const std::string& identifier, int numberOfVariables,
    int numberOfParameters, int nodesPerCoordinateAxis, double maximumMeshSize,
    exahype::solvers::Solver::TimeStepping timeStepping,
    std::unique_ptr<profilers::Profiler> profiler)
    : Solver(identifier, exahype::solvers::Solver::Type::FiniteVolumes,
             numberOfVariables, numberOfParameters, nodesPerCoordinateAxis,
             maximumMeshSize, timeStepping, std::move(profiler)),
      _unknownsPerCell((numberOfVariables + numberOfParameters) *
                       power(nodesPerCoordinateAxis, DIMENSIONS + 0)),
      _unknownsPerFace(
          (numberOfVariables)*power(nodesPerCoordinateAxis, DIMENSIONS - 1)),
      _unknownsPerCellBoundary(
          DIMENSIONS_TIMES_TWO *
          (numberOfVariables)*power(nodesPerCoordinateAxis, DIMENSIONS - 1)),
      _minTimeStamp(std::numeric_limits<double>::max()),
      _minTimeStepSize(std::numeric_limits<double>::max()),
      _nextMinTimeStepSize(std::numeric_limits<double>::max()) {
  assertion3(_unknownsPerCell > 0, numberOfVariables, numberOfParameters,
             nodesPerCoordinateAxis);

  // register tags with profiler
  for (const char* tag : tags) {
    _profiler->registerTag(tag);
  }
}

int exahype::solvers::FiniteVolumesSolver::getUnknownsPerCell() const {
  return _unknownsPerCell;
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
    const int cellDescriptionsIndex,
    const int element) {
  assertionMsg(false,"Please implement!");
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
  assertionMsg(false,"Not implemented.");
  return false;
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
  assertionMsg(false,"Not implemented.");
  return false;
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////
double exahype::solvers::FiniteVolumesSolver::startNewTimeStep(
    const int cellDescriptionsIndex,
    const int element) {
  // do nothing (for now). See discussion within body of function mappings/SolutionUpdate::enterCell(...).
  return std::numeric_limits<double>::max();
}

void exahype::solvers::FiniteVolumesSolver::setInitialConditions(
    const int cellDescriptionsIndex,
    const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  // reset helper variables
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
//  exahype::Cell::resetNeighbourMergeHelperVariables(
//      cellDescription,fineGridVertices,fineGridVerticesEnumerator); // TODO(Dominic): Add flags.

  if (cellDescription.getType()==CellDescription::Cell
//      && cellDescription.getRefinementEvent()==CellDescription::None
      ) {
    double* luh = exahype::DataHeap::getInstance().getData(cellDescription.getSolution()).data();

    if (hasToAdjustSolution(
        cellDescription.getOffset()+0.5*cellDescription.getSize(),
        cellDescription.getSize(),
        cellDescription.getTimeStamp())) {
      solutionAdjustment(
          luh,
          cellDescription.getOffset()+0.5*cellDescription.getSize(),
          cellDescription.getSize(),
          cellDescription.getTimeStamp(), cellDescription.getTimeStepSize());
    }

    for (int i=0; i<getUnknownsPerCell(); i++) {
      assertion3(std::isfinite(luh[i]),cellDescription.toString(),"setInitialConditions(...)",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
  }
}

void exahype::solvers::FiniteVolumesSolver::updateSolution(
    const int cellDescriptionsIndex,
    const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  // reset helper variables
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
//  exahype::Cell::resetNeighbourMergeHelperVariables(
//      cellDescription,fineGridVertices,fineGridVerticesEnumerator); // TODO(Dominic): Add helper variables.

  // TODO(Dominic): update solution (This should be separated in mergeNeighbours and updateSolution
  const tarch::la::Vector<THREE_POWER_D, int> neighbourCellDescriptionsIndices =
      multiscalelinkedcell::getIndicesAroundCell(
          exahype::VertexOperations::readCellDescriptionsIndex(fineGridVerticesEnumerator, fineGridVertices));

  // todo MPI
  // todo Boundary
  #ifdef SharedTBB
  assertionMsg(false,"Not implemented yet!");
  #endif
  assertion1(multiscalelinkedcell::HangingVertexBookkeeper::allAdjacencyInformationIsAvailable(
      VertexOperations::readCellDescriptionsIndex(fineGridVerticesEnumerator, fineGridVertices)),
             fineGridVerticesEnumerator.toString());

  double* finiteVolumeSolutions[THREE_POWER_D];
  for (int nScalar=0; nScalar<THREE_POWER_D; ++nScalar) {
    if (Heap::getInstance().isValidIndex(neighbourCellDescriptionsIndices[nScalar])) {
      CellDescription& pNeighbour =
          Heap::getInstance().getData(neighbourCellDescriptionsIndices[nScalar])[cellDescription.getSolverNumber()]; // todo assumes same number of patches per cell
      finiteVolumeSolutions[nScalar] = DataHeap::getInstance().getData(pNeighbour.getSolution()).data();
    } else {
      finiteVolumeSolutions[nScalar] = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
    }
  }

  double* finiteVolumeSolution  = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
  assertion(!std::isnan(finiteVolumeSolution[0]));

  double admissibleTimeStepSize=0;
  solutionUpdate(finiteVolumeSolutions,cellDescription.getSize(),cellDescription.getTimeStepSize(),admissibleTimeStepSize);

  if (admissibleTimeStepSize < cellDescription.getTimeStepSize()) {
    logWarning("updateSolution(...)","Finite volumes solver time step size harmed CFL condition. dt="<<cellDescription.getTimeStepSize()<<", dt_adm=" << admissibleTimeStepSize);
  }

  if (hasToAdjustSolution(
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getTimeStamp())) {
    solutionAdjustment(
        finiteVolumeSolution,
        cellDescription.getOffset()+0.5*cellDescription.getSize(),
        cellDescription.getSize(),
        cellDescription.getTimeStamp(), cellDescription.getTimeStepSize());
  }

  for (int i=0; i<getUnknownsPerCell(); i++) {
    assertion3(std::isfinite(finiteVolumeSolution[i]),cellDescription.toString(),"finiteVolumeSolution[i]",i);
  } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
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
      const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  assertionMsg(false,"Not implemented.");
}

void exahype::solvers::FiniteVolumesSolver::mergeWithBoundaryData(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary) {
  assertionMsg(false,"Not implemented.");
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

void exahype::solvers::FiniteVolumesSolver::sendToRank(int rank, int tag) {
  assertionMsg(false, "not implemented yet");
}

void exahype::solvers::FiniteVolumesSolver::receiveFromMasterRank(int rank, int tag) {
  assertionMsg(false, "not implemented yet");
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
  out << "_unknownsPerFace:" << _unknownsPerFace;
  out << ",";
  out << "_unknownsPerCellBoundary:" << _unknownsPerCellBoundary;
  out << ",";
  out << "_unknownsPerCell:" << _unknownsPerCell;
  out << ",";
  out << "_minTimeStamp:" << _minTimeStamp;
  out << ",";
  out << "_minTimeStepSize:" << _minTimeStepSize;
  out << ",";
  out << "_nextMinTimeStepSize:" << _nextMinTimeStepSize;
  out <<  ")";
}
