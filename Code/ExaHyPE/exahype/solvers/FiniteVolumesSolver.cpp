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

namespace {
constexpr const char* tags[]{"solutionUpdate", "stableTimeStepSize"};
}  // namespace

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

#ifdef Parallel
const int exahype::solvers::FiniteVolumesSolver::DataMessagesPerNeighbour = 1;

void exahype::solvers::FiniteVolumesSolver::sendToRank(int rank, int tag) {
  assertionMsg(false, "not implemented yet");
}

void exahype::solvers::FiniteVolumesSolver::receiveFromMasterRank(int rank, int tag) {
  assertionMsg(false, "not implemented yet");
}

void exahype::solvers::FiniteVolumesSolver::sendEmptyDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    int                                           level) {
  std::vector<double> emptyMessage(0,0);
  for(int sends=0; sends<DataMessagesPerNeighbour; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
}

void exahype::solvers::FiniteVolumesSolver::dropNeighbourData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    int                                           level) {
  for(int receives=0; receives<DataMessagesPerNeighbour; ++receives)
    DataHeap::getInstance().receiveData(
        fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
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
