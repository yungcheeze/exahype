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
 
#include "exahype/Cell.h"

#include "exahype/solvers/ADERDGSolver.h"

#include "tarch/la/VectorVectorOperations.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

namespace {
constexpr const char* tags[]{"solutionUpdate",
                             "volumeIntegral",
                             "surfaceIntegral",
                             "riemannSolver",
                             "spaceTimePredictor",
                             "stableTimeStepSize",
                             "solutionAdjustment",
                             "faceUnknownsProlongation",
                             "faceUnknownsRestriction",
                             "volumeUnknownsProlongation",
                             "volumeUnknownsRestriction",
                             "boundaryConditions"};
}  // namespace

exahype::solvers::ADERDGSolver::ADERDGSolver(
    const std::string& identifier, int numberOfVariables,
    int numberOfParameters, int nodesPerCoordinateAxis, double maximumMeshSize,
    exahype::solvers::Solver::TimeStepping timeStepping,
    std::unique_ptr<profilers::Profiler> profiler)
    : Solver(identifier, Solver::Type::ADER_DG, numberOfVariables,
             numberOfParameters, nodesPerCoordinateAxis, maximumMeshSize,
             timeStepping, std::move(profiler)),
      _unknownsPerFace(numberOfVariables *
                       power(nodesPerCoordinateAxis, DIMENSIONS - 1)),
      _unknownsPerCellBoundary(DIMENSIONS_TIMES_TWO * _unknownsPerFace),
      _unknownsPerCell(numberOfVariables *
                       power(nodesPerCoordinateAxis, DIMENSIONS + 0)),
      _fluxUnknownsPerCell(_unknownsPerCell *
                           (DIMENSIONS + 1)),  // +1 for sources
      _spaceTimeUnknownsPerCell(numberOfVariables *
                                power(nodesPerCoordinateAxis, DIMENSIONS + 1)),
      _spaceTimeFluxUnknownsPerCell(_spaceTimeUnknownsPerCell *
                                    (DIMENSIONS + 1)),  // +1 for sources
      _dataPerCell(numberOfVariables *
                   power(nodesPerCoordinateAxis, DIMENSIONS + 0)),
      _minCorrectorTimeStamp(std::numeric_limits<double>::max()),
      _minPredictorTimeStamp(std::numeric_limits<double>::max()),
      _minCorrectorTimeStepSize(std::numeric_limits<double>::max()),
      _minPredictorTimeStepSize(std::numeric_limits<double>::max()),
      _minNextPredictorTimeStepSize(std::numeric_limits<double>::max()) {
  assertion(numberOfParameters == 0);
  // register tags with profiler
  for (const char* tag : tags) {
    _profiler->registerTag(tag);
  }
}

std::string exahype::solvers::ADERDGSolver::getIdentifier() const {
  return _identifier;
}

exahype::solvers::ADERDGSolver::Type exahype::solvers::ADERDGSolver::getType() const {
  return _type;
}

int exahype::solvers::ADERDGSolver::getNumberOfVariables() const {
  return _numberOfVariables;
}

int exahype::solvers::ADERDGSolver::getNumberOfParameters() const {
  return _numberOfParameters;
}

int exahype::solvers::ADERDGSolver::getNodesPerCoordinateAxis() const {
  return _nodesPerCoordinateAxis;
}

double exahype::solvers::ADERDGSolver::getMaximumMeshSize() const {
  return _maximumMeshSize;
}

int exahype::solvers::ADERDGSolver::getUnknownsPerFace() const {
  return _unknownsPerFace;
}

int exahype::solvers::ADERDGSolver::getUnknownsPerCellBoundary() const {
  return _unknownsPerCellBoundary;
}

int exahype::solvers::ADERDGSolver::getUnknownsPerCell() const {
  return _unknownsPerCell;
}

int exahype::solvers::ADERDGSolver::getFluxUnknownsPerCell() const {
  return _fluxUnknownsPerCell;
}

int exahype::solvers::ADERDGSolver::getSpaceTimeUnknownsPerCell() const {
  return _spaceTimeUnknownsPerCell;
}

int exahype::solvers::ADERDGSolver::getSpaceTimeFluxUnknownsPerCell() const {
  return _spaceTimeFluxUnknownsPerCell;
}

int exahype::solvers::ADERDGSolver::getDataPerCell() const {
  return _dataPerCell;
}

void exahype::solvers::ADERDGSolver::synchroniseTimeStepping(
    exahype::records::ADERDGCellDescription& p) const {
  switch (_timeStepping) {
    case TimeStepping::Global:
      p.setCorrectorTimeStamp(_minCorrectorTimeStamp);
      p.setCorrectorTimeStepSize(_minCorrectorTimeStepSize);
      p.setPredictorTimeStamp(_minPredictorTimeStamp);
      p.setPredictorTimeStepSize(_minPredictorTimeStepSize);
      break;
    case TimeStepping::GlobalFixed:
      p.setCorrectorTimeStamp(_minCorrectorTimeStamp);
      p.setCorrectorTimeStepSize(_minCorrectorTimeStepSize);
      p.setPredictorTimeStamp(_minPredictorTimeStamp);
      p.setPredictorTimeStepSize(_minPredictorTimeStepSize);
      break;
  }
}

void exahype::solvers::ADERDGSolver::synchroniseTimeStepping(
      const int cellDescriptionsIndex,
      const int element) {
  synchroniseTimeStepping(Heap::getInstance().getData(cellDescriptionsIndex)[element]);
}


void exahype::solvers::ADERDGSolver::startNewTimeStep() {
  switch (_timeStepping) {
    case TimeStepping::Global:
      _minCorrectorTimeStamp    = _minPredictorTimeStamp;
      _minCorrectorTimeStepSize = _minPredictorTimeStepSize;

      _minPredictorTimeStepSize = _minNextPredictorTimeStepSize;
      _minPredictorTimeStamp    = _minPredictorTimeStamp + _minNextPredictorTimeStepSize;

      _minNextPredictorTimeStepSize = std::numeric_limits<double>::max();
      break;
    case TimeStepping::GlobalFixed:
      _minCorrectorTimeStamp    = _minPredictorTimeStamp;
      _minCorrectorTimeStepSize = _minPredictorTimeStepSize;

      _minPredictorTimeStepSize = _minNextPredictorTimeStepSize;
      _minPredictorTimeStamp    = _minPredictorTimeStamp + _minNextPredictorTimeStepSize;
      break;
  }
}

void exahype::solvers::ADERDGSolver::updateMinNextPredictorTimeStepSize(
    const double& minNextPredictorTimeStepSize) {
  switch (_timeStepping) {
    case TimeStepping::Global:
      _minNextPredictorTimeStepSize =
          std::min(_minNextPredictorTimeStepSize, minNextPredictorTimeStepSize);
      break;
    case TimeStepping::GlobalFixed:
      _minNextPredictorTimeStepSize =
          _minNextPredictorTimeStepSize == std::numeric_limits<double>::max()
              ? std::min(_minNextPredictorTimeStepSize,
                         minNextPredictorTimeStepSize)
              : _minNextPredictorTimeStepSize;
      break;
  }
}

double exahype::solvers::ADERDGSolver::getMinNextPredictorTimeStepSize() const {
  return _minNextPredictorTimeStepSize;
}

void exahype::solvers::ADERDGSolver::setMinCorrectorTimeStamp(
    double minCorrectorTimeStamp) {
  _minCorrectorTimeStamp = minCorrectorTimeStamp;
}

double exahype::solvers::ADERDGSolver::getMinCorrectorTimeStamp() const {
  return _minCorrectorTimeStamp;
}

void exahype::solvers::ADERDGSolver::setMinPredictorTimeStamp(
    double minPredictorTimeStamp) {
  _minPredictorTimeStamp = minPredictorTimeStamp;
}

double exahype::solvers::ADERDGSolver::getMinPredictorTimeStamp() const {
  return _minPredictorTimeStamp;
}

double exahype::solvers::ADERDGSolver::getMinCorrectorTimeStepSize() const {
  return _minCorrectorTimeStepSize;
}

void exahype::solvers::ADERDGSolver::setMinPredictorTimeStepSize(
    double minPredictorTimeStepSize) {
  _minPredictorTimeStepSize = minPredictorTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getMinPredictorTimeStepSize() const {
  return _minPredictorTimeStepSize;
}

int exahype::solvers::ADERDGSolver::tryGetElement(
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

void exahype::solvers::ADERDGSolver::mergeNeighbours(
    const int                                     cellDescriptionsIndex1,
    const int                                     element1,
    const int                                     cellDescriptionsIndex2,
    const int                                     element2,
    const tarch::la::Vector<DIMENSIONS, int>&     pos1,
    const tarch::la::Vector<DIMENSIONS, int>&     pos2) {
  if (tarch::la::countEqualEntries(pos1,pos2)!=1) {
    return; // We only consider faces; no corners.
  }
  // !!! In Riemann solve we consider "left" face of "right" cell and
  // "right" face of "left" cell. !!!
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

  auto& pLeft  = Heap::getInstance().getData(cellDescriptionsIndexLeft) [elementLeft];
  auto& pRight = Heap::getInstance().getData(cellDescriptionsIndexRight)[elementRight];

  solveRiemannProblemAtInterface(pLeft,pRight,faceIndexLeft,faceIndexRight);

  mergeSolutionMinMaxOnFace(pLeft,pRight,faceIndexLeft,faceIndexRight);
}

void exahype::solvers::ADERDGSolver::solveRiemannProblemAtInterface(
    exahype::records::ADERDGCellDescription& pLeft,
    exahype::records::ADERDGCellDescription& pRight,
    const int faceIndexLeft,
    const int faceIndexRight) {
  if (pLeft.getType()  == exahype::records::ADERDGCellDescription::Cell ||
      pRight.getType() == exahype::records::ADERDGCellDescription::Cell) {
    assertion1(pLeft.getType() == exahype::records::ADERDGCellDescription::Cell ||
        pLeft.getType() == exahype::records::ADERDGCellDescription::Ancestor ||
        pLeft.getType() == exahype::records::ADERDGCellDescription::Descendant,
        pLeft.toString());
    assertion1(pRight.getType() == exahype::records::ADERDGCellDescription::Cell ||
        pRight.getType() == exahype::records::ADERDGCellDescription::Ancestor ||
        pRight.getType() == exahype::records::ADERDGCellDescription::Descendant,
        pRight.toString());
    assertion1(pLeft.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,
        pLeft.toString());
    assertion1(pRight.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,
        pRight.toString());
    assertionEquals4(pLeft.getRiemannSolvePerformed(faceIndexLeft),
        pRight.getRiemannSolvePerformed(faceIndexRight),
        faceIndexLeft, faceIndexRight,
        pLeft.toString(),
        pRight.toString());
    exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*> (
        exahype::solvers::RegisteredSolvers[pLeft.getSolverNumber()] );

    const int numberOfFaceDof = solver->getUnknownsPerFace();

    double* QL = DataHeap::getInstance() .getData(pLeft.getExtrapolatedPredictor()).data() +
        (faceIndexLeft * numberOfFaceDof);
    double* QR = DataHeap::getInstance().getData(pRight.getExtrapolatedPredictor()).data() +
        (faceIndexRight * numberOfFaceDof);
    double* FL = DataHeap::getInstance().getData(pLeft.getFluctuation()).data() +
        (faceIndexLeft * numberOfFaceDof);
    double* FR = DataHeap::getInstance().getData(pRight.getFluctuation()).data() +
        (faceIndexRight * numberOfFaceDof);

    for(int ii=0; ii<numberOfFaceDof; ++ii) {
      assertion(std::isfinite(QL[ii]));
      assertion(std::isfinite(QR[ii]));
      assertion(std::isfinite(FL[ii]));
      assertion(std::isfinite(FR[ii]));
    }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

    // Synchronise time stepping.
    solver->synchroniseTimeStepping(pLeft);
    solver->synchroniseTimeStepping(pRight);

    // todo Time step must be interpolated in local time stepping case
    // both time step sizes are the same, so the min has no effect here.
    assertion1(faceIndexRight%2==0,faceIndexRight);
    const int normalDirection = (faceIndexRight - (faceIndexRight %2));

    solver->riemannSolver(
        FL, FR, QL, QR,
        std::min(pLeft.getCorrectorTimeStepSize(),
            pRight.getCorrectorTimeStepSize()),
            normalDirection);

    for(int i=0; i<numberOfFaceDof; ++i) {
      assertion(std::isfinite(FL[i]));
      assertion(std::isfinite(FR[i]));
    }  // Dead code elimination will get rid of this loop if Asserts flag is not set.
  }
}

void exahype::solvers::ADERDGSolver::mergeWithBoundaryData(
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, int>&     posCell,
      const tarch::la::Vector<DIMENSIONS, int>&     posBoundary) {
  if (tarch::la::countEqualEntries(posCell,posBoundary)!=1) {
    return; // We only consider faces; no corners.
  }

  // !!! Left face of right cell.
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(posCell, posBoundary);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (posCell(normalOfExchangedFace) < posBoundary(normalOfExchangedFace) ? 1 : 0);

  exahype::records::ADERDGCellDescription& p = Heap::getInstance().getData(cellDescriptionsIndex)[element];
  applyBoundaryConditions(p,faceIndex);
}

// Verified correct calling of this method for 9x9 grid on [0,1]x[0,1].
void exahype::solvers::ADERDGSolver::applyBoundaryConditions(
    exahype::records::ADERDGCellDescription& p,
    const int faceIndex) {
  assertion1(p.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,p.toString());
  const int numberOfFaceDof = getUnknownsPerFace();

  double* stateIn = DataHeap::getInstance().getData(p.getExtrapolatedPredictor()).data() +
      (faceIndex * numberOfFaceDof);
  double* fluxIn = DataHeap::getInstance().getData(p.getFluctuation()).data() +
      (faceIndex * numberOfFaceDof);

  const int normalDirection = (faceIndex - faceIndex % 2)/2;
  for(int ii=0; ii<numberOfFaceDof; ++ii) {
    assertion5(std::isfinite(stateIn[ii]), p.toString(),
        faceIndex, normalDirection, ii, stateIn[ii]);
    assertion5(std::isfinite(fluxIn[ii]), p.toString(),
        faceIndex, normalDirection, ii, fluxIn[ii]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  double* stateOut = new double[numberOfFaceDof];
  double* fluxOut  = new double[numberOfFaceDof];

  // Synchronise time stepping.
  synchroniseTimeStepping(p);

  boundaryConditions(fluxOut,stateOut,
      fluxIn,stateIn,
      p.getOffset() + 0.5*p.getSize(), // centre
      p.getSize(),
      p.getCorrectorTimeStamp(),
      p.getCorrectorTimeStepSize(),
      faceIndex,
      normalDirection);

  for(int ii=0; ii<numberOfFaceDof; ++ii) {
    assertion5(std::isfinite(stateOut[ii]), p.toString(), faceIndex, normalDirection, ii, stateOut[ii]);
    assertion5(std::isfinite(fluxOut[ii]), p.toString(), faceIndex, normalDirection, ii, fluxOut[ii]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  // @todo(Dominic): Add to docu why we need this. Left or right input
  if (faceIndex % 2 == 0) {
    riemannSolver(fluxOut, fluxIn, stateOut, stateIn,
        p.getCorrectorTimeStepSize(),
        normalDirection);
  } else {
    riemannSolver(fluxIn, fluxOut, stateIn, stateOut,
        p.getCorrectorTimeStepSize(),
        normalDirection);
  }

  for(int ii=0; ii<numberOfFaceDof; ++ii) {
    assertion5(std::isfinite(fluxIn[ii]), p.toString(),
               faceIndex, normalDirection, ii, fluxIn[ii]);
    assertion5(std::isfinite(fluxOut[ii]), p.toString(),
               faceIndex, normalDirection, ii, fluxOut[ii]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  delete[] stateOut;
  delete[] fluxOut;
}

void exahype::solvers::ADERDGSolver::mergeSolutionMinMaxOnFace(
  exahype::records::ADERDGCellDescription& pLeft,
  exahype::records::ADERDGCellDescription& pRight,
  const int faceIndexLeft,
  const int faceIndexRight
) const {
  if (
      (pLeft.getType() == exahype::records::ADERDGCellDescription::Cell
          && (pRight.getType() == exahype::records::ADERDGCellDescription::Cell ||
              pRight.getType() == exahype::records::ADERDGCellDescription::Ancestor ||
              pRight.getType() == exahype::records::ADERDGCellDescription::Descendant))
              ||
              (pRight.getType() == exahype::records::ADERDGCellDescription::Cell
                  && (pLeft.getType() == exahype::records::ADERDGCellDescription::Cell ||
                      pLeft.getType() == exahype::records::ADERDGCellDescription::Ancestor ||
                      pLeft.getType() == exahype::records::ADERDGCellDescription::Descendant))
  ) {
    assertion( pLeft.getSolverNumber() == pRight.getSolverNumber() );
    assertion( exahype::solvers::RegisteredSolvers[ pLeft.getSolverNumber() ]->getType()==exahype::solvers::Solver::Type::ADER_DG );
    const int numberOfVariables = static_cast<exahype::solvers::ADERDGSolver*>(
        exahype::solvers::RegisteredSolvers[ pLeft.getSolverNumber() ])->getNumberOfVariables();
    for (int i=0; i<numberOfVariables; i++) {
      double min = std::min(
          DataHeap::getInstance().getData( pLeft.getSolutionMin()  )[i+faceIndexLeft *numberOfVariables],
          DataHeap::getInstance().getData( pRight.getSolutionMin() )[i+faceIndexRight*numberOfVariables]
      );
      double max = std::min(
          DataHeap::getInstance().getData( pLeft.getSolutionMax()  )[i+faceIndexLeft *numberOfVariables],
          DataHeap::getInstance().getData( pRight.getSolutionMax() )[i+faceIndexRight*numberOfVariables]
      );

      DataHeap::getInstance().getData( pLeft.getSolutionMin()  )[i+faceIndexLeft *numberOfVariables] = min;
      DataHeap::getInstance().getData( pRight.getSolutionMin() )[i+faceIndexRight*numberOfVariables] = min;

      DataHeap::getInstance().getData( pLeft.getSolutionMax()  )[i+faceIndexLeft *numberOfVariables] = max;
      DataHeap::getInstance().getData( pRight.getSolutionMax() )[i+faceIndexRight*numberOfVariables] = max;
    }
  } // else do nothing
}

#ifdef Parallel
const int exahype::solvers::ADERDGSolver::DataMessagesPerNeighbourCommunication    = 3;
const int exahype::solvers::ADERDGSolver::DataMessagesPerForkOrJoinCommunication   = 1;
const int exahype::solvers::ADERDGSolver::DataMessagesPerMasterWorkerCommunication = 2;

void exahype::solvers::ADERDGSolver::sendCellDescriptions(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),
      cellDescriptionsIndex);

  for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
    switch(p.getType()) {
      case exahype::records::ADERDGCellDescription::Descendant:
      case exahype::records::ADERDGCellDescription::EmptyDescendant:
        // TODO(Dominic): Change to RemoteBoundaryDescendant
        p.setType(exahype::records::ADERDGCellDescription::Descendant);
        exahype::Cell::ensureNecessaryMemoryIsAllocated(p);
        break;
      case exahype::records::ADERDGCellDescription::Ancestor:
      case exahype::records::ADERDGCellDescription::EmptyAncestor:
        // TODO(Dominic): Change to RemoteBoundaryAncestor if top most parent
        // is ancestor or RemoteBoundaryAncestor
        p.setType(exahype::records::ADERDGCellDescription::Ancestor);
        exahype::Cell::ensureNecessaryMemoryIsAllocated(p);
        break;
      default:
        break;
    }
  }

  Heap::getInstance().sendData(cellDescriptionsIndex,
          toRank,x,level,messageType);
}

void exahype::solvers::ADERDGSolver::sendEmptyCellDescriptions(
    const int                                     toRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::HeapEntries emptyMessage(0);
  Heap::getInstance().sendData(emptyMessage,
      toRank,x,level,messageType);
}

void exahype::solvers::ADERDGSolver::mergeCellDescriptionsWithRemoteData(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),
      cellDescriptionsIndex);

  int receivedCellDescriptionsIndex =
      Heap::getInstance().createData(0,exahype::solvers::RegisteredSolvers.size());
  Heap::getInstance().receiveData(cellDescriptionsIndex,
      fromRank,x,level,messageType);
  resetDataHeapIndices(
      receivedCellDescriptionsIndex,
      multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex);

  Heap::getInstance().getData(cellDescriptionsIndex).reserve(
      std::max(Heap::getInstance().getData(cellDescriptionsIndex).size(),
               Heap::getInstance().getData(receivedCellDescriptionsIndex).size()));

  for (auto& pReceived : Heap::getInstance().getData(receivedCellDescriptionsIndex)) {
    bool found = false;
    for (auto& pLocal : Heap::getInstance().getData(cellDescriptionsIndex)) {
      if (pReceived.getSolverNumber()==pLocal.getSolverNumber()) {
        found = true;

        assertion(pReceived.getType()==pLocal.getType());
        if (pLocal.getType()==exahype::records::ADERDGCellDescription::Type::Cell ||
            pLocal.getType()==exahype::records::ADERDGCellDescription::Type::Ancestor ||
//            pLocal.getType()==exahype::records::ADERDGCellDescription::Type::RemoteBoundaryAncestor ||
//            pLocal.getType()==exahype::records::ADERDGCellDescription::Type::RemoteBoundaryDescendant ||
            pLocal.getType()==exahype::records::ADERDGCellDescription::Type::Descendant
            ) {
          assertionNumericalEquals2(pLocal.getCorrectorTimeStamp(),pReceived.getCorrectorTimeStamp(),
              pLocal.toString(),pReceived.toString());
          assertionNumericalEquals2(pLocal.getCorrectorTimeStepSize(),pReceived.getCorrectorTimeStepSize(),
              pLocal.toString(),pReceived.toString());
          assertionNumericalEquals2(pLocal.getPredictorTimeStamp(),pReceived.getPredictorTimeStamp(),
              pLocal.toString(),pReceived.toString());
          assertionNumericalEquals2(pLocal.getPredictorTimeStepSize(),pReceived.getPredictorTimeStepSize(),
              pLocal.toString(),pReceived.toString());
        }
      }
    }

    if (!found) {
      exahype::Cell::ensureNecessaryMemoryIsAllocated(pReceived);
      Heap::getInstance().getData(cellDescriptionsIndex).
          push_back(pReceived);
    }
  }
}

void exahype::solvers::ADERDGSolver::resetDataHeapIndices(
    const int cellDescriptionsIndex,
    const int parentIndex) {
  for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
    p.setParentIndex(parentIndex);

    // Default field data indices
    p.setSolution(-1);
    p.setUpdate(-1);
    p.setExtrapolatedPredictor(-1);
    p.setFluctuation(-1);

    // Limiter meta data (oscillations identificator)
    p.setSolutionMin(-1);
    p.setSolutionMax(-1);
  }
}

/**
 * Drop cell descriptions received from \p fromRank.
 */
void exahype::solvers::ADERDGSolver::dropCellDescriptions(
    const int                                     fromRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::getInstance().receiveData(fromRank,x,level,messageType);
}

//std::vector<double> exahype::solvers::ADERDGSolver::collectTimeStampsAndStepSizes() {
//  std::vector<double> timeStampsAndStepSizes(0,5);
//
//  timeStampsAndStepSizes.push_back(_minCorrectorTimeStamp);
//  timeStampsAndStepSizes.push_back(_minCorrectorTimeStepSize);
//  timeStampsAndStepSizes.push_back(_minPredictorTimeStepSize);
//  timeStampsAndStepSizes.push_back(_minPredictorTimeStamp);
//  timeStampsAndStepSizes.push_back(_minNextPredictorTimeStepSize);
//  return timeStampsAndStepSizes;
//}
//
//void exahype::solvers::ADERDGSolver::setTimeStampsAndStepSizes(
//    std::vector<double>& timeSteppingData) {
//  _minCorrectorTimeStamp        = timeSteppingData[0];
//  _minCorrectorTimeStepSize     = timeSteppingData[1];
//  _minPredictorTimeStepSize     = timeSteppingData[2];
//  _minPredictorTimeStamp        = timeSteppingData[3];
//  _minNextPredictorTimeStepSize = timeSteppingData[4];
//}

void exahype::solvers::ADERDGSolver::sendToRank(int rank, int tag) {
  MPI_Send(&_minCorrectorTimeStamp, 1, MPI_DOUBLE, rank, tag,
      tarch::parallel::Node::getInstance().getCommunicator()); // This is necessary since we might have performed a predictor rerun.
  MPI_Send(&_minCorrectorTimeStepSize, 1, MPI_DOUBLE, rank, tag,
      tarch::parallel::Node::getInstance().getCommunicator()); // This is necessary since we might have performed a predictor rerun.
  MPI_Send(&_minPredictorTimeStepSize, 1, MPI_DOUBLE, rank, tag,
      tarch::parallel::Node::getInstance().getCommunicator());
  MPI_Send(&_minPredictorTimeStamp, 1, MPI_DOUBLE, rank, tag,
      tarch::parallel::Node::getInstance().getCommunicator());
  MPI_Send(&_minNextPredictorTimeStepSize, 1, MPI_DOUBLE, rank, tag,
      tarch::parallel::Node::getInstance().getCommunicator());
}

void exahype::solvers::ADERDGSolver::receiveFromMasterRank(int rank, int tag) {
  MPI_Recv(&_minCorrectorTimeStamp, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator(),
           MPI_STATUS_IGNORE); // This is necessary since we might have performed a predictor rerun.
  MPI_Recv(&_minCorrectorTimeStepSize, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator(),
           MPI_STATUS_IGNORE); // This is necessary since we might have performed a predictor rerun.
  MPI_Recv(&_minPredictorTimeStepSize, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator(),
           MPI_STATUS_IGNORE);
  MPI_Recv(&_minPredictorTimeStamp, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator(),
           MPI_STATUS_IGNORE);
  MPI_Recv(&_minNextPredictorTimeStepSize, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator(),
           MPI_STATUS_IGNORE);
}

///////////////////////////////////
// FORK OR JOIN
///////////////////////////////////

void exahype::solvers::ADERDGSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  auto& p = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  double* solution = 0;
  switch(p.getType()) {
    case exahype::records::ADERDGCellDescription::Cell:
      solution        = DataHeap::getInstance().getData(p.getSolution()).data();

      logDebug("sendDataToWorkerOrMasterDueToForkOrJoin(...)","solution of solver " << p.getSolverNumber() << " sent to rank "<<toRank<<
               ", cell: "<< cellCentre << ", level: " << level);

      DataHeap::getInstance().sendData(
          solution, getUnknownsPerCell(), toRank, x, level,
          peano::heap::MessageType::ForkOrJoinCommunication);
      break;
    default:
      break;
  }
}


void exahype::solvers::ADERDGSolver::sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerForkOrJoinCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, toRank, x, level,
        peano::heap::MessageType::ForkOrJoinCommunication);
}


void exahype::solvers::ADERDGSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  auto& p = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex)[element];
  assertion4(tarch::la::equals(x,p.getOffset()+0.5*p.getSize()),x,p.getOffset()+0.5*p.getSize(),level,p.getLevel());
  assertion2(p.getLevel()==level,p.getLevel(),level);

  if (p.getType()==exahype::records::ADERDGCellDescription::Cell) {
    logDebug("mergeWithRemoteDataDueToForkOrJoin(...)","[solution] receive from rank "<<fromRank<<
             ", cell: "<< cellCentre << ", level: " << level);

    DataHeap::getInstance().receiveData(
        p.getSolution(),fromRank,x,level,
        peano::heap::MessageType::ForkOrJoinCommunication);
  }
}

void exahype::solvers::ADERDGSolver::dropWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  for(int receives=0; receives<DataMessagesPerForkOrJoinCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        fromRank, x, level,
        peano::heap::MessageType::ForkOrJoinCommunication);
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////

void exahype::solvers::ADERDGSolver::sendDataToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  if (tarch::la::countEqualEntries(src,dest)!=1) {
    return; // We only consider faces; no corners.
  }

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!

  auto& p = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  if (holdsFaceData(p.getType())) {
    assertion(DataHeap::getInstance().isValidIndex(p.getExtrapolatedPredictor()));
    assertion(DataHeap::getInstance().isValidIndex(p.getFluctuation()));
    assertion(DataHeap::getInstance().isValidIndex(p.getSolutionMin()));
    assertion(DataHeap::getInstance().isValidIndex(p.getSolutionMax()));

    const int numberOfFaceDof = getUnknownsPerFace();
    const double* lQhbnd = DataHeap::getInstance().getData(
        p.getExtrapolatedPredictor()).data() +
        (faceIndex * numberOfFaceDof);
    const double* lFhbnd = DataHeap::getInstance().getData(
        p.getFluctuation()).data() +
        (faceIndex * numberOfFaceDof);

    logDebug(
        "sendDataToNeighbour(...)",
        "send "<<exahype::Cell::DataExchangesPerADERDGSolver<<" arrays to rank " <<
        toRank << " for cell="<<p.getOffset()<< " and face=" << faceIndex << " from vertex x=" << x << ", level=" << level <<
        ", src type=" << multiscalelinkedcell::indexToString(srcCellDescriptionIndex) <<
        ", dest type=" << multiscalelinkedcell::indexToString(destCellDescriptionIndex) <<
        ", src=" << src << ", dest=" << dest <<
        ", adjacent ranks=" << vertex.getAdjacentRanks() <<
        ", counter=" << p.getFaceDataExchangeCounter(faceIndex)
    );

    // We append all the max values to the min values.
    std::vector<double> sentMinMax( 2*getNumberOfVariables() );
    for (int i=0; i<getNumberOfVariables(); i++) {
      sentMinMax[i]                                = DataHeap::getInstance().getData( p.getSolutionMin() )[faceIndex*getNumberOfVariables()+i];
      sentMinMax[i+getNumberOfVariables()] = DataHeap::getInstance().getData( p.getSolutionMax() )[faceIndex*getNumberOfVariables()+i];
    }

    // Send order: minMax,lQhbnd,lFhbnd
    // Receive order: lFhbnd,lQhbnd,minMax
    DataHeap::getInstance().sendData(
        sentMinMax, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
    DataHeap::getInstance().sendData(
        lQhbnd, numberOfFaceDof, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
    DataHeap::getInstance().sendData(
        lFhbnd, numberOfFaceDof, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
    // TODO(Dominic): If anarchic time stepping send the time step over too.
  } else {
    std::vector<double> emptyArray(0,0);

    for(int sends=0; sends<DataMessagesPerNeighbourCommunication; ++sends) {
      DataHeap::getInstance().sendData(
          emptyArray, toRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
    }
  }
}

void exahype::solvers::ADERDGSolver::sendEmptyDataToNeighbour(
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

void exahype::solvers::ADERDGSolver::mergeWithNeighbourData(
    const int                                     fromRank,
    const int                                     neighbourTypeAsInt,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  if (tarch::la::countEqualEntries(src,dest)!=1) {
    return; // We only consider faces; no corners.
  }

  auto& p = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  exahype::records::ADERDGCellDescription::Type neighbourType =
      static_cast<exahype::records::ADERDGCellDescription::Type>(neighbourTypeAsInt);

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) > dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

  if(holdsFaceData(neighbourType) && holdsFaceData(p.getType())){
    assertion4(!p.getRiemannSolvePerformed(faceIndex),
        faceIndex,cellDescriptionsIndex,p.getOffset().toString(),p.getLevel());
    assertion(DataHeap::getInstance().isValidIndex(p.getExtrapolatedPredictor()));
    assertion(DataHeap::getInstance().isValidIndex(p.getFluctuation()));

    logDebug(
        "mergeNeighbourData(...)", "receive ">>DataMessagesPerNeighbourCommunication>>" arrays from rank " <<
        fromRank << " for vertex x=" << x << ", level=" << level <<
        ", src type=" << multiscalelinkedcell::indexToString(srcCellDescriptionIndex) <<
        ", src=" << src << ", dest=" << dest <<
        ", counter=" << p.getFaceDataExchangeCounter(faceIndex)
    );

    const int numberOfFaceDof = getUnknownsPerFace();
    int receivedlQhbndIndex = DataHeap::getInstance().createData(0, numberOfFaceDof);
    int receivedlFhbndIndex = DataHeap::getInstance().createData(0, numberOfFaceDof);
    int receivedMinMax      = DataHeap::getInstance().createData(0, 2*getNumberOfVariables());

    assertion(DataHeap::getInstance().getData(receivedlQhbndIndex).empty());
    assertion(DataHeap::getInstance().getData(receivedlFhbndIndex).empty());
    assertion(DataHeap::getInstance().getData(receivedMinMax).empty());

    // Send order: minMax,lQhbnd,lFhbnd
    // Receive order: lFhbnd,lQhbnd,minMax
    DataHeap::getInstance().receiveData(receivedlFhbndIndex, fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
    DataHeap::getInstance().receiveData(receivedlQhbndIndex, fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
    DataHeap::getInstance().receiveData(receivedMinMax,  fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);

    logDebug(
        "receiveADERDGFaceData(...)", "[pre] solve Riemann problem with received data." <<
        " cellDescription=" << p.toString() <<
        ",faceIndexForCell=" << faceIndex <<
        ",normalOfExchangedFac=" << normalOfExchangedFace <<
        ",x=" << x.toString() << ", level=" << level <<
        ", counter=" << p.getFaceDataExchangeCounter(faceIndex)
    );

    solveRiemannProblemAtInterface(
        p,
        faceIndex,
        receivedlQhbndIndex,
        receivedlFhbndIndex);

    mergeSolutionMinMaxOnFace(
        p,
        faceIndex,
        DataHeap::getInstance().getData(receivedMinMax).data(),
        DataHeap::getInstance().getData(receivedMinMax).data() + getNumberOfVariables() );

    // TODO(Dominic): If anarchic time stepping, receive the time step too.

    DataHeap::getInstance().deleteData(receivedlQhbndIndex);
    DataHeap::getInstance().deleteData(receivedlFhbndIndex);
    DataHeap::getInstance().deleteData(receivedMinMax);
  } else  {
    logDebug(
        "receiveADERDGFaceData(...)", "drop three arrays from rank " <<
        fromRank << " for vertex x=" << x << ", level=" << level <<
        ", src type=" << multiscalelinkedcell::indexToString(srcCellDescriptionIndex) <<
        ", src=" << src << ", dest=" << dest <<
        ", counter=" << p.getFaceDataExchangeCounter(faceIndex)
    );

    dropNeighbourData(fromRank,src,dest,x,level);
  }
}

void exahype::solvers::ADERDGSolver::solveRiemannProblemAtInterface(
    records::ADERDGCellDescription& cellDescription,
    const int faceIndex,
    const int indexOfQValues,
    const int indexOfFValues) {
  cellDescription.setRiemannSolvePerformed(faceIndex, true);

  const int numberOfFaceDof = getUnknownsPerFace();

  logDebug("solveRiemannProblemAtInterface(...)",
      "cell-description=" << cellDescription.toString());

  double* QL = 0;
  double* QR = 0;
  double* FL = 0;
  double* FR = 0;

  assertionEquals(DataHeap::getInstance().getData(indexOfQValues).size(),
      static_cast<unsigned int>(numberOfFaceDof));
  assertionEquals(DataHeap::getInstance().getData(indexOfFValues).size(),
      static_cast<unsigned int>(numberOfFaceDof));

  // @todo Doku im Header warum wir das hier brauchen,
  if (faceIndex % 2 == 0) {
    QL = DataHeap::getInstance().getData(indexOfQValues).data();
    QR = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
        (faceIndex * numberOfFaceDof);
    FL = DataHeap::getInstance().getData(indexOfFValues).data();
    FR = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
        (faceIndex * numberOfFaceDof);
  } else {
    QR = DataHeap::getInstance().getData(indexOfQValues).data();
    QL = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
        (faceIndex * numberOfFaceDof);
    FR = DataHeap::getInstance().getData(indexOfFValues).data();
    FL = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
        (faceIndex * numberOfFaceDof);
  }

  // Synchronise time stepping.
  synchroniseTimeStepping(cellDescription);

  const int normalDirection = (faceIndex - faceIndex%2)/2; // faceIndex=2*normalNonZero+f, f=0,1
  riemannSolver(FL, FR, QL, QR,
      cellDescription.getCorrectorTimeStepSize(),
      normalDirection);

  for (int ii = 0; ii<numberOfFaceDof; ii++) {
    assertion8(std::isfinite(QR[ii]), cellDescription.toString(),
        faceIndex, normalDirection, indexOfQValues, indexOfFValues,
        ii, QR[ii], QL[ii]);
    assertion8(std::isfinite(QL[ii]), cellDescription.toString(),
        faceIndex, normalDirection, indexOfQValues, indexOfFValues,
        ii, QR[ii], QL[ii]);
    assertion8(std::isfinite(FR[ii]), cellDescription.toString(),
        faceIndex, normalDirection, indexOfQValues, indexOfFValues,
        ii, QR[ii], QL[ii]);
    assertion8(std::isfinite(FL[ii]), cellDescription.toString(),
        faceIndex, normalDirection, indexOfQValues, indexOfFValues,
        ii, QR[ii], QL[ii]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.
}

void exahype::solvers::ADERDGSolver::mergeSolutionMinMaxOnFace(
  exahype::records::ADERDGCellDescription&  cellDescription,
  int                                       faceIndex,
  double* min, double* max) const {
  if (cellDescription.getType() == exahype::records::ADERDGCellDescription::Cell ||
      cellDescription.getType() == exahype::records::ADERDGCellDescription::Ancestor ||
      cellDescription.getType() == exahype::records::ADERDGCellDescription::Descendant
      ) {
    assertion( exahype::solvers::RegisteredSolvers[ cellDescription.getSolverNumber() ]->getType()==exahype::solvers::Solver::Type::ADER_DG );
    const int numberOfVariables = static_cast<exahype::solvers::ADERDGSolver*>(
        exahype::solvers::RegisteredSolvers[ cellDescription.getSolverNumber() ])->getNumberOfVariables();

    for (int i=0; i<numberOfVariables; i++) {
      DataHeap::getInstance().getData( cellDescription.getSolutionMin()  )[i+faceIndex*numberOfVariables]  = min[i];
      DataHeap::getInstance().getData( cellDescription.getSolutionMax()  )[i+faceIndex*numberOfVariables]  = max[i];
    }
  }
}

void exahype::solvers::ADERDGSolver::dropNeighbourData(
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
// WORKER->MASTER
///////////////////////////////////

void exahype::solvers::ADERDGSolver::sendDataToMaster(
    const int                                     masterRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  auto& p = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  if (p.getType()==exahype::records::ADERDGCellDescription::Ancestor
      // || exahype::records::ADERDGCellDescription::RemoteBoundaryAncestor: // TODO(Dominic)
  ) {
    double* extrapolatedPredictor = DataHeap::getInstance().getData(p.getExtrapolatedPredictor()).data();
    double* fluctuations          = DataHeap::getInstance().getData(p.getFluctuation()).data();

    logDebug("sendDataToWorkerOrMasterDueToForkOrJoin(...)","solution of solver " << p.getSolverNumber() << " sent to rank "<<masterRank<<
             ", cell: "<< cellCentre << ", level: " << level);

    // !!! Be aware of inverted receive order !!!
    // Send order:    extrapolatedPredictor, fluctuations
    // Receive order: fluctuations, extrapolatedPredictor
    DataHeap::getInstance().sendData(
        extrapolatedPredictor, getUnknownsPerCellBoundary(), masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
    DataHeap::getInstance().sendData(
        fluctuations, getUnknownsPerCellBoundary(), masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);

  }
}

void exahype::solvers::ADERDGSolver::sendEmptyDataToMaster(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerMasterWorkerCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::ADERDGSolver::mergeWithWorkerData(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  auto& p = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  // TODO(Dominic): RemoteBoundaryAncestor should be only valid type at
  // master-worker boundary. Change both ancestor and empty ancestor
  // to this type if necessary.
  if (p.getType()==exahype::records::ADERDGCellDescription::Ancestor
      // || exahype::records::ADERDGCellDescription::RemoteBoundaryAncestor:
  ) {
    logDebug("sendDataToWorkerOrMasterDueToForkOrJoin(...)","solution of solver " << p.getSolverNumber() << " sent to rank "<<masterRank<<
             ", cell: "<< cellCentre << ", level: " << level);

    // !!! Be aware of inverted receive order !!!
    // Send order:    extrapolatedPredictor, fluctuations
    // Receive order: fluctuations, extrapolatedPredictor
    DataHeap::getInstance().receiveData(
        p.getFluctuation(), workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
    DataHeap::getInstance().receiveData(
        p.getExtrapolatedPredictor(), workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);

  }
}

void exahype::solvers::ADERDGSolver::dropWorkerData(
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

void exahype::solvers::ADERDGSolver::sendDataToWorker(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  auto& p = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  if (p.getType()==exahype::records::ADERDGCellDescription::Descendant
      // || exahype::records::ADERDGCellDescription::RemoteBoundaryDescendant: // TODO(Dominic)
  ) {
    double* extrapolatedPredictor = DataHeap::getInstance().getData(p.getExtrapolatedPredictor()).data();
    double* fluctuations          = DataHeap::getInstance().getData(p.getFluctuation()).data();

    logDebug("sendDataToWorkerOrMasterDueToForkOrJoin(...)","solution of solver " << p.getSolverNumber() << " sent to rank "<<masterRank<<
             ", cell: "<< cellCentre << ", level: " << level);

    // !!! Be aware of inverted receive order !!!
    // Send order:    extrapolatedPredictor, fluctuations
    // Receive order: fluctuations, extrapolatedPredictor
    DataHeap::getInstance().sendData(
        extrapolatedPredictor, getUnknownsPerCellBoundary(), workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
    DataHeap::getInstance().sendData(
        fluctuations, getUnknownsPerCellBoundary(), workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);

  }
}

void exahype::solvers::ADERDGSolver::sendEmptyDataToWorker(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerMasterWorkerCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::ADERDGSolver::mergeWithMasterData(
    const int                                     masterRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  auto& p = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  // TODO(Dominic): RemoteBoundaryAncestor should be only valid type at
  // master-worker boundary. Change both ancestor and empty ancestor
  // to this type if necessary.
  if (p.getType()==exahype::records::ADERDGCellDescription::Descendant
      // || exahype::records::ADERDGCellDescription::RemoteBoundaryDescendant:
  ) {
    logDebug("sendDataToWorkerOrMasterDueToForkOrJoin(...)","solution of solver " << p.getSolverNumber() << " sent to rank "<<masterRank<<
             ", cell: "<< cellCentre << ", level: " << level);

    // !!! Be aware of inverted receive order !!!
    // Send order:    extrapolatedPredictor, fluctuations
    // Receive order: fluctuations, extrapolatedPredictor
    DataHeap::getInstance().receiveData(
        p.getFluctuation(), masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
    DataHeap::getInstance().receiveData(
        p.getExtrapolatedPredictor(), masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);

  }
}

void exahype::solvers::ADERDGSolver::dropMasterData(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
        const int                                     level) {
  for(int receives=0; receives<DataMessagesPerMasterWorkerCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
}
#endif

std::string exahype::solvers::ADERDGSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::ADERDGSolver::toString (std::ostream& out) const {
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
  out << "_fluxUnknownsPerCell:" << _fluxUnknownsPerCell;
  out << ",";
  out << "_spaceTimeUnknownsPerCell:" << _spaceTimeUnknownsPerCell;
  out << ",";
  out << "_spaceTimeFluxUnknownsPerCell:" << _spaceTimeFluxUnknownsPerCell;
  out << ",";
  out << "_minCorrectorTimeStamp:" << _minCorrectorTimeStamp;
  out << ",";
  out << "_minPredictorTimeStamp:" << _minPredictorTimeStamp;
  out << ",";
  out << "_minCorrectorTimeStepSize:" << _minCorrectorTimeStepSize;
  out << ",";
  out << "_minPredictorTimeStepSize:" << _minPredictorTimeStepSize;
  out << ",";
  out << "_minNextPredictorTimeStepSize:" << _minNextPredictorTimeStepSize;
  out <<  ")";
}
