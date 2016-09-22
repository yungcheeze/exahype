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
 
#include "exahype/Vertex.h"
#include "peano/utils/Loop.h"
#include "peano/grid/Checkpoint.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/State.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"


exahype::Vertex::Vertex() : Base() {
  _vertexData._persistentRecords._CellDescriptionsIndex =
      multiscalelinkedcell::HangingVertexBookkeeper::getInstance()
          .createVertexLinkMapForNewVertex();
}

exahype::Vertex::Vertex(const Base::DoNotCallStandardConstructor& value)
    : Base(value) {
  // do nothing
}

exahype::Vertex::Vertex(const Base::PersistentVertex& argument)
    : Base(argument) {
  // do nothing
}

tarch::la::Vector<TWO_POWER_D, int>&
exahype::Vertex::getCellDescriptionsIndex() {
  return _vertexData._persistentRecords._CellDescriptionsIndex;
}

bool exahype::Vertex::hasToMergeNeighbours(
      const tarch::la::Vector<DIMENSIONS,int>& pos1,
      const tarch::la::Vector<DIMENSIONS,int>& pos2) const {
  if (tarch::la::countEqualEntries(pos1,pos2)==1) {
    const int pos1Scalar = peano::utils::dLinearisedWithoutLookup(pos1,2);
    const int pos2Scalar = peano::utils::dLinearisedWithoutLookup(pos2,2);
    const int cellDescriptionsIndex1 =
        _vertexData._persistentRecords._CellDescriptionsIndex[pos1Scalar];
    const int cellDescriptionsIndex2 =
        _vertexData._persistentRecords._CellDescriptionsIndex[pos2Scalar];

    if (exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex1) &&
        exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex2)) {
      assertion1(pos1Scalar!=pos2Scalar,pos1Scalar);
      assertion1(cellDescriptionsIndex1!=cellDescriptionsIndex2,cellDescriptionsIndex1);
      assertion1(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex1),
          cellDescriptionsIndex1);
      assertion1(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex2),
          cellDescriptionsIndex2);

      const int normalOfExchangedFace = tarch::la::equalsReturnIndex(pos1, pos2);
      assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
      const int faceIndex1 = 2 * normalOfExchangedFace +
          (pos2(normalOfExchangedFace) > pos1(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!
      const int faceIndex2 = 2 * normalOfExchangedFace +
          (pos1(normalOfExchangedFace) > pos2(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

      // Here, we check all cell descriptions.
      for (auto& p1 : exahype::solvers::
          ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex1)) {
        if (p1.getRiemannSolvePerformed(faceIndex1)) {
          return false;
        }
      }
      // TODO(Dominic): Finite Volumes.

      for (auto& p2 : exahype::solvers::
          ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex2)) {
        if (p2.getRiemannSolvePerformed(faceIndex2)) {
          return false;
        }
      }
      // TODO(Dominic): Finite Volumes.

      return true;
    }
  }

  return false;
}

bool exahype::Vertex::hasToMergeWithBoundaryData(
      const tarch::la::Vector<DIMENSIONS,int>& pos1,
      const tarch::la::Vector<DIMENSIONS,int>& pos2) const {
  if (tarch::la::countEqualEntries(pos1,pos2)==1) {
    const int pos1Scalar = peano::utils::dLinearisedWithoutLookup(pos1,2);
    const int pos2Scalar = peano::utils::dLinearisedWithoutLookup(pos2,2);
    const int cellDescriptionsIndex1 =
        _vertexData._persistentRecords._CellDescriptionsIndex[pos1Scalar];
    const int cellDescriptionsIndex2 =
        _vertexData._persistentRecords._CellDescriptionsIndex[pos2Scalar];

    const bool validIndexNextToBoundaryOrRemoteIndex =
        (exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex1)
            &&
            (cellDescriptionsIndex2==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex
            #if !defined(PeriodBC)
            || cellDescriptionsIndex2==multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex
            #endif
            ))
            ||
            (exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex2)
                &&
                (cellDescriptionsIndex1==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex
                #if !defined(PeriodBC)
                || cellDescriptionsIndex1==multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex
                #endif
                ));

    if (validIndexNextToBoundaryOrRemoteIndex) {
      const int normalOfExchangedFace = tarch::la::equalsReturnIndex(pos1, pos2);
      assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
      const int faceIndex1 = 2 * normalOfExchangedFace +
          (pos2(normalOfExchangedFace) > pos1(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!
      const int faceIndex2 = 2 * normalOfExchangedFace +
          (pos1(normalOfExchangedFace) > pos2(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

      if (exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex1)) {
        for (auto& p1 : exahype::solvers::
            ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex1)) {
          if (!p1.getIsInside(faceIndex1)
              && !p1.getRiemannSolvePerformed(faceIndex1)) {
            return true;
          }
        }
        // TODO(Dominic): Finite volumes
      } else if (exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex2)) {
        for (auto& p2 : exahype::solvers::
            ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex2)) {
          if (!p2.getIsInside(faceIndex2)
              && !p2.getRiemannSolvePerformed(faceIndex2)) {
            return true;
          }
        }
        // TODO(Dominic): Finite volumes
      }
    }
  }
  return false;
}

void exahype::Vertex::setMergePerformed(
        const tarch::la::Vector<DIMENSIONS,int>& pos1,
        const tarch::la::Vector<DIMENSIONS,int>& pos2,
        bool state) const {
  if (tarch::la::countEqualEntries(pos1,pos2)!=1) {
    return; // We only consider faces; no corners.
  }

  const int pos1Scalar = peano::utils::dLinearisedWithoutLookup(pos1,2);
  const int pos2Scalar = peano::utils::dLinearisedWithoutLookup(pos2,2);
  const int cellDescriptionsIndex1 =
      _vertexData._persistentRecords._CellDescriptionsIndex[pos1Scalar];
  const int cellDescriptionsIndex2 =
      _vertexData._persistentRecords._CellDescriptionsIndex[pos2Scalar];

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(pos1, pos2);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex1 = 2 * normalOfExchangedFace +
      (pos2(normalOfExchangedFace) > pos1(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!
  const int faceIndex2 = 2 * normalOfExchangedFace +
      (pos1(normalOfExchangedFace) > pos2(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

  if (exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex1)) {
    for (auto& p1 : exahype::solvers::
        ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex1)) {
      p1.setRiemannSolvePerformed(faceIndex1,state);
    }
  }
  // TODO(Dominic): Finite volumes are missing.

  if (exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex2)) {
    for (auto& p2 : exahype::solvers::
        ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex2)) {
      p2.setRiemannSolvePerformed(faceIndex2,state);
    }
  }
  // TODO(Dominic): Finite volumes are missing.
}

#ifdef Parallel
const int exahype::Vertex::InvalidMetadataEntry   = -1;

const int exahype::Vertex::MetadataPerSolver      =  1;

exahype::MetadataHeap::HeapEntries exahype::Vertex::encodeMetadata(int cellDescriptionsIndex) {
  assertion1(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);

  exahype::MetadataHeap::HeapEntries encodedMetaData(
      exahype::solvers::RegisteredSolvers.size()*MetadataPerSolver,
      exahype::solvers::RegisteredSolvers.size()*MetadataPerSolver);
  std::fill_n(encodedMetaData.begin(),encodedMetaData.size(),InvalidMetadataEntry); // Implicit conversion.

  // ADER-DG
  for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex)) {
    encodedMetaData[p.getSolverNumber()] = static_cast<int>(p.getType()); // Implicit conversion.
  }
  // FV
  for (auto& p : exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex)) {
    encodedMetaData[p.getSolverNumber()] = static_cast<int>(p.getType()); // Implicit conversion.
  }
  return encodedMetaData;
}

bool exahype::Vertex::hasToSendMetadata(
  const tarch::la::Vector<DIMENSIONS,int>& src,
  const tarch::la::Vector<DIMENSIONS,int>& dest,
  const int toRank) {
  const int srcScalar  = peano::utils::dLinearisedWithoutLookup(src,2);
  const int destScalar = peano::utils::dLinearisedWithoutLookup(dest,2);
  const tarch::la::Vector<TWO_POWER_D,int> adjacentRanks = getAdjacentRanks();

  return tarch::la::countEqualEntries(dest, src)  == 1 &&
         adjacentRanks(destScalar)   == toRank &&
         (adjacentRanks(srcScalar)   == tarch::parallel::Node::getInstance().getRank() ||
         State::isForkTriggeredForRank(adjacentRanks(srcScalar)));
}


bool exahype::Vertex::hasToSendMetadataIgnoreForksJoins(
  const tarch::la::Vector<DIMENSIONS,int>& src,
  const tarch::la::Vector<DIMENSIONS,int>& dest,
  const int toRank) {
  const int srcScalar  = peano::utils::dLinearisedWithoutLookup(src,2);
  const int destScalar = peano::utils::dLinearisedWithoutLookup(dest,2);
  const tarch::la::Vector<TWO_POWER_D,int> adjacentRanks = getAdjacentRanks();

  return tarch::la::countEqualEntries(dest, src)  == 1 &&
         adjacentRanks(destScalar) == toRank &&
         adjacentRanks(srcScalar)  == tarch::parallel::Node::getInstance().getRank();
}

bool exahype::Vertex::hasToReceiveMetadata(
  const tarch::la::Vector<DIMENSIONS,int>& src,
  const tarch::la::Vector<DIMENSIONS,int>& dest,
  const int fromRank) {
  const int srcScalar  = peano::utils::dLinearisedWithoutLookup(src,2);
  const int destScalar = peano::utils::dLinearisedWithoutLookup(dest,2);
  const tarch::la::Vector<TWO_POWER_D,int> adjacentRanks = getAdjacentRanks();

  return tarch::la::countEqualEntries(dest, src) == 1 &&
      adjacentRanks(srcScalar)    == fromRank &&
      (adjacentRanks(destScalar)  == tarch::parallel::Node::getInstance().getRank() ||
       State::isForkingRank(adjacentRanks(destScalar)));
}

bool exahype::Vertex::hasToReceiveMetadataIgnoreForksJoins(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const int fromRank) {
  const int srcScalar  = peano::utils::dLinearisedWithoutLookup(src,2);
  const int destScalar = peano::utils::dLinearisedWithoutLookup(dest,2);
  const tarch::la::Vector<TWO_POWER_D,int> adjacentRanks = getAdjacentRanks();

  return tarch::la::countEqualEntries(dest, src) == 1 &&
      adjacentRanks(srcScalar)   == fromRank &&
      adjacentRanks(destScalar)  == tarch::parallel::Node::getInstance().getRank();
}

bool exahype::Vertex::hasToSendDataToNeighbour(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest) const {
  const int srcScalar  = peano::utils::dLinearisedWithoutLookup(src,2);
  const int srcCellDescriptionsIndex =
      _vertexData._persistentRecords._CellDescriptionsIndex[srcScalar];

  if (!exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionsIndex) ||
      exahype::solvers::ADERDGSolver::Heap::getInstance().getData(srcCellDescriptionsIndex).empty()) {
    return false;
  }

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!

  // ADER-DG
  for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(srcCellDescriptionsIndex)) {
    if (
        #if !defined(PeriodicBC)
        !p.getIsInside(faceIndex) ||
        #endif
        p.getFaceDataExchangeCounter(faceIndex)!=0) {
      return false;
    }
  }
  //  // TODO(Dominic): Introduce counters to FiniteVolumesCellDescription.
  //  // FV
  return true;
}

bool exahype::Vertex::hasToMergeWithNeighbourData(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest) const {
  const int destScalar  = peano::utils::dLinearisedWithoutLookup(dest,2);
  const int destCellDescriptionsIndex =
      _vertexData._persistentRecords._CellDescriptionsIndex[destScalar];

  if (!exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionsIndex) ||
      exahype::solvers::ADERDGSolver::Heap::getInstance().getData(destCellDescriptionsIndex).empty()) {
    return false;
  }

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) > dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

  // ADER-DG
  for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(destCellDescriptionsIndex)) {
    if (
        #if !defined(PeriodicBC)
        !p.getIsInside(faceIndex) ||
        #endif
        p.getFaceDataExchangeCounter(faceIndex)!=0) {
      return false;
    }
  }
  //  // TODO(Dominic): Introduce counters to FiniteVolumesCellDescription.
  //  // FV
  return true;
}

void exahype::Vertex::tryDecrementFaceDataExchangeCountersOfSource(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest) const {
  const int srcScalar  = peano::utils::dLinearisedWithoutLookup(src,2);
  const int srcCellDescriptionsIndex =
      _vertexData._persistentRecords._CellDescriptionsIndex[srcScalar];

  if (!exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionsIndex) ||
      exahype::solvers::ADERDGSolver::Heap::getInstance().getData(srcCellDescriptionsIndex).empty())
    return;
  assertion1(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().
      isValidIndex(srcCellDescriptionsIndex),
      srcCellDescriptionsIndex);

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!

  // ADER-DG
  for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(srcCellDescriptionsIndex)) {
    int newCounterValue = p.getFaceDataExchangeCounter(faceIndex)-1;
    assertion1(newCounterValue<TWO_POWER_D,newCounterValue);
    p.setFaceDataExchangeCounter(faceIndex,newCounterValue);
  }

//  // TODO(Dominic): Introduce counters to FiniteVolumesCellDescription.
//  // FV
}

void exahype::Vertex::setFaceDataExchangeCountersOfDestination(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const int value) const {
  const int destScalar  = peano::utils::dLinearisedWithoutLookup(dest,2);
  const int destCellDescriptionsIndex =
      _vertexData._persistentRecords._CellDescriptionsIndex[destScalar];

  assertion1(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionsIndex),destCellDescriptionsIndex);
  assertion1(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(destCellDescriptionsIndex),destCellDescriptionsIndex);

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) > dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

  // ADER-DG
  for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(destCellDescriptionsIndex)) {
    p.setFaceDataExchangeCounter(faceIndex,value);
  }

//  // TODO(Dominic): Introduce counters to FiniteVolumesCellDescription.
//  // FV
}

#endif
