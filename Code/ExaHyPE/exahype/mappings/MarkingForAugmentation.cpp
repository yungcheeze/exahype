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
 
#include "exahype/mappings/MarkingForAugmentation.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "exahype/VertexOperations.h"
#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/solvers/Solver.h"
#include "kernels/KernelCalls.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::MarkingForAugmentation::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::MarkingForAugmentation::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::MarkingForAugmentation::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::MarkingForAugmentation::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::MarkingForAugmentation::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::MarkingForAugmentation::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::MarkingForAugmentation::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::MarkingForAugmentation::_log(
    "exahype::mappings::MarkingForAugmentation");

void exahype::mappings::MarkingForAugmentation::beginIteration(
  exahype::State& solverState) {
}

void exahype::mappings::MarkingForAugmentation::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
//#ifdef Parallel
//  return;
//#endif

  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (fineGridCell.isInitialised() &&
      multiscalelinkedcell::adjacencyInformationIsConsistent(
          VertexOperations::readCellDescriptionsIndex(fineGridVerticesEnumerator, fineGridVertices))
  ) {
    const tarch::la::Vector<THREE_POWER_D, int>
        neighbourCellDescriptionIndices =
            multiscalelinkedcell::getIndicesAroundCell(
                VertexOperations::readCellDescriptionsIndex(
                    fineGridVerticesEnumerator, fineGridVertices));

    bool refineFineGridCell = false;
    for (auto& pFine : ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())) {
      const AugmentationControl augmentationControl =
          augmentationCriterion(pFine.getSolverNumber(), pFine.getType(),
                                pFine.getLevel(),
                                neighbourCellDescriptionIndices);

      // 1. Check if ancestors and descendants need to hold
      // data or not based on virtual refinement criterion.
      // Then, allocate the necessary memory and deallocate the unnecessary memory.
      switch (pFine.getType()) {
        case exahype::records::ADERDGCellDescription::Ancestor:
        case exahype::records::ADERDGCellDescription::EmptyAncestor:
          pFine.setType(exahype::records::ADERDGCellDescription::EmptyAncestor);
          #ifdef Parallel
          if (pFine.getHelperCellNeedsToStoreFaceData()) { // TODO(Dominic): Add to docu.
            pFine.setType(exahype::records::ADERDGCellDescription::Ancestor);
          }
          #endif
          if (fineGridCell.isAssignedToRemoteRank() ||
              pFine.getParentIndex()==multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex) { // TODO(Dominic): Add to docu.
            pFine.setType(exahype::records::ADERDGCellDescription::Ancestor);
          }

          if (augmentationControl==AugmentationControl::NextToCell ||
              augmentationControl==AugmentationControl::NextToCellAndAncestor) {
              pFine.setType(exahype::records::ADERDGCellDescription::Ancestor);
          }
          fineGridCell.ensureNoUnnecessaryMemoryIsAllocated(pFine.getSolverNumber());
          fineGridCell.ensureNecessaryMemoryIsAllocated(pFine.getSolverNumber());
          break;
        case exahype::records::ADERDGCellDescription::Descendant:
        case exahype::records::ADERDGCellDescription::EmptyDescendant:
          pFine.setType(exahype::records::ADERDGCellDescription::EmptyDescendant);
          #ifdef Parallel
          if (pFine.getHelperCellNeedsToStoreFaceData()) { // TODO(Dominic): Add to docu.
            pFine.setType(exahype::records::ADERDGCellDescription::Descendant);
          }
          #endif
          if (fineGridCell.isAssignedToRemoteRank() ||
              pFine.getParentIndex()==multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex) { // TODO(Dominic): Add to docu.
            pFine.setType(exahype::records::ADERDGCellDescription::Descendant);
          }
          if (augmentationControl==AugmentationControl::NextToCell ||
              augmentationControl==AugmentationControl::NextToCellAndAncestor) {
            pFine.setType(exahype::records::ADERDGCellDescription::Descendant);
          }
          fineGridCell.ensureNecessaryMemoryIsAllocated(pFine.getSolverNumber());
          fineGridCell.ensureNoUnnecessaryMemoryIsAllocated(pFine.getSolverNumber());
          break;
        default:
          break;
      }

      // 2. Further augment or deaugment cells and descendants if no other event
      // or an augmentation event has been triggered.
      switch (pFine.getRefinementEvent()) {
        case exahype::records::ADERDGCellDescription::AugmentingRequested:
          refineFineGridCell = true;
          break;
        case exahype::records::ADERDGCellDescription::None:
          switch (pFine.getType()) {
            case exahype::records::ADERDGCellDescription::Cell:
              switch (augmentationControl) {
                case AugmentationControl::NextToAncestor:
                case AugmentationControl::NextToCellAndAncestor:
                  pFine.setRefinementEvent(exahype::records::ADERDGCellDescription::AugmentingRequested);
                  refineFineGridCell = true;
                  break;
                default:
                  break;
              }
              break;
            case exahype::records::ADERDGCellDescription::Descendant:
            case exahype::records::ADERDGCellDescription::EmptyDescendant:
              switch (augmentationControl) {
                case AugmentationControl::NextToAncestor:
                case AugmentationControl::NextToCellAndAncestor:
                  pFine.setRefinementEvent(exahype::records::ADERDGCellDescription::AugmentingRequested);
                  refineFineGridCell = true;
                  break;
                case AugmentationControl::Default:
                  pFine.setRefinementEvent(exahype::records::ADERDGCellDescription::DeaugmentingRequested);
                  break;
                default:
                  break;
              }
              break;
            default:
              break;
          }
          break;
        default:
          break;
      }

      // TODO(Dominic): Add to docu why we reset this flag here, and where it initialised (Cell.addNew...).
      #ifdef Parallel
      pFine.setHelperCellNeedsToStoreFaceData(false);
      #endif
    }

    // Refine the vertices
//    if (refineFineGridCell && _state->refineInitialGridInTouchVertexLastTime()) { // TODO what is going on here?
    if (refineFineGridCell) {
      dfor2(v)
        if (fineGridVertices[ fineGridVerticesEnumerator(v) ].getRefinementControl()==
            exahype::Vertex::Records::RefinementControl::Unrefined) {
          fineGridVertices[ fineGridVerticesEnumerator(v) ].refine();
        }
      enddforx
    }
  }

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

exahype::mappings::MarkingForAugmentation::AugmentationControl
exahype::mappings::MarkingForAugmentation::augmentationCriterion(
    const int solverNumber,
    const exahype::records::ADERDGCellDescription::Type type, const int level,
    const tarch::la::Vector<THREE_POWER_D, int>&
        neighbourCellDescriptionIndices) {
// left,right,front,back,(front,back)
#if DIMENSIONS == 2
  constexpr int neighbourPositions[4] = {3, 5, 1, 7};
#else
  constexpr int neighbourPositions[6] = {12, 14, 10, 16, 4, 22};
#endif
  bool nextToAncestor = false;
  bool nextToCell = false;

  for (int i = 0; i < DIMENSIONS_TIMES_TWO; i++) {
    const int neighbourCellDescriptionIndex =
        neighbourCellDescriptionIndices[neighbourPositions[i]];
    if (DataHeap::getInstance().isValidIndex(neighbourCellDescriptionIndex)) {
      for (auto& pNeighbour : ADERDGCellDescriptionHeap::getInstance().getData(
               neighbourCellDescriptionIndex)) {
        if (pNeighbour.getSolverNumber() == solverNumber &&
            pNeighbour.getLevel() == level) {
          switch (pNeighbour.getType()) {
            case exahype::records::ADERDGCellDescription::Ancestor:
            case exahype::records::ADERDGCellDescription::EmptyAncestor:
              nextToAncestor = true;
              break;
            case exahype::records::ADERDGCellDescription::Cell:
              nextToCell = true;
              break;
            default:
              break;
          }
        }
      }
    }
  }

  // NOTE: The order below is important.
  if (nextToCell && nextToAncestor) {
    return AugmentationControl::NextToCellAndAncestor;
  }
  if (nextToAncestor) {
    return AugmentationControl::NextToAncestor;
  }
  if (nextToCell) {
    return AugmentationControl::NextToCell;
  }
  // Erase otherwise.
  return AugmentationControl::Default;
}

#ifdef Parallel
void exahype::mappings::MarkingForAugmentation::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  logTraceInWith5Arguments("prepareSendToNeighbour(...)", vertex,
                           toRank, x, h, level);

  // TODO(Dominic): remove; We should be able to remove this and vertex.isInside()
  // We should further use (first touch decrements; decrement only once) the same counters
  // as for the spaceTimePredictor to reduce the number of long messages!={0.0} by a factor of four.
  #if !defined(PeriodicBC)
  if (vertex.isBoundary()) return;
  #endif

  tarch::la::Vector<TWO_POWER_D, int>& adjacentADERDGCellDescriptionsIndices =
      vertex.getCellDescriptionsIndex();

  if (vertex.isInside()) { // TODO(Dominic): Add to docu: This prevents metadata exchange with rank 0.
    dfor2(dest)
      dfor2(src)
      if (vertex.hasToSendMetadata(src,dest,toRank)) {
          const int srcCellDescriptionIndex = adjacentADERDGCellDescriptionsIndices(srcScalar);
          if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(srcCellDescriptionIndex)) {
            logDebug("prepareSendToNeighbour(...)","[data] sent to rank "<<toRank<<", x:"<<
                x.toString() << ", level=" <<level << ", vertex.adjacentRanks: "
                << vertex.getAdjacentRanks()
                << ", src forking: "
                << State::isForkingRank(vertex.getAdjacentRanks()(srcScalar)));

            exahype::MetadataHeap::HeapEntries metadata =
                exahype::Cell::encodeMetadata(srcCellDescriptionIndex);
            MetadataHeap::getInstance().sendData(
                metadata, toRank, x, level,
                peano::heap::MessageType::NeighbourCommunication);
          } else {
            if (srcCellDescriptionIndex==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex) {
              assertion2(false,"Should not happen!",srcCellDescriptionIndex);
            } // Dead code elimination will remove this code section if Asserts/Debug flag not set.
            if (srcCellDescriptionIndex>0) {
              // TODO(Dominic): Encountered situation where srcCellDescriptionsIndex > 0 but
              // index was not valid. I assume this happens if a fork was triggered. In this
              // case I sent out an empty message. Add to docu.
              // Just sending empty data out seems to work fine so far.
            }

            logDebug("prepareSendToNeighbour(...)","[empty] sent to rank "<<toRank<<", x:"<<
                x.toString() << ", level=" <<level << ", vertex.adjacentRanks: "
                << vertex.getAdjacentRanks()
                << ", src forking: "
                << State::isForkingRank(vertex.getAdjacentRanks()(srcScalar)));

            exahype::MetadataHeap::HeapEntries metadata =
                exahype::Cell::createEncodedMetadataSequenceForInvalidCellDescriptionsIndex();

            MetadataHeap::getInstance().sendData(
                metadata, toRank, x, level,
                peano::heap::MessageType::NeighbourCommunication);
          }
        }
      enddforx
    enddforx
  }

  logTraceOut("prepareSendToNeighbour(...)");
}

void exahype::mappings::MarkingForAugmentation::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments("mergeWithNeighbour(...)", vertex, neighbour,
                           fromRank, fineGridX, fineGridH, level);

  // TODO(Dominic): remove
//  return;

  // TODO(Dominic): AMR + MPI
  // 1. Get metadata,
  // 2. Set AugmentationRequest if neighbour is of type Ancestor
  // 3. Change cell type of local cell description if neighbour is of type Cell.
  // 4. Delete metadata.
#if !defined(PeriodicBC)
  if (vertex.isBoundary()) return;
#endif
  // TODO(Dominic): Add to docu: mergeWithNeighbour(..) happens before vertex creation events.
  // TODO(Dominic): Might need to change this (Periodic BC). Might need to send data between
  // rank 0 and other ranks too. So vertex.isOutside()/Boundary() must be possible too.
  if (vertex.isInside()) {
    assertion1(vertex.isInside(),vertex.toString());
    assertion1(neighbour.isInside(),neighbour.toString());

    tarch::la::Vector<TWO_POWER_D, int>& adjacentADERDGCellDescriptionsIndices =
          vertex.getCellDescriptionsIndex();

    dfor2(myDest)
      dfor2(mySrc)
        // TODO(Dominic): Add to docu why we invert the order:
        // MPI message order: Stack principle.
        tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
        tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;
        int destScalar = TWO_POWER_D - myDestScalar - 1;

        if (vertex.hasToReceiveMetadata(src,dest,fromRank)) {
          logDebug("mergeWithNeighbour(...)","[pre] rec. from rank "<<fromRank<<", x:"<<
                                          fineGridX.toString() << ", level=" <<level << ", vertex.adjacentRanks: "
                                          << vertex.getAdjacentRanks());

          int receivedMetadataIndex = MetadataHeap::getInstance().createData(0,0);
          assertion(MetadataHeap::getInstance().getData(receivedMetadataIndex).empty());
          MetadataHeap::getInstance().receiveData(
              receivedMetadataIndex,
              fromRank, fineGridX, level,
              peano::heap::MessageType::NeighbourCommunication);

          const int destCellDescriptionIndex = adjacentADERDGCellDescriptionsIndices(destScalar);
          if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(destCellDescriptionIndex)) {
            assertion(FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(destCellDescriptionIndex));

            // We do not have to invert the order here since we receive the solver metadata in correct
            // order and perform local actions according to these metadata. No further messages need
            // to be handled by this vertex.
            receiveADERDGMetadataInMergeWithNeigbour       (destCellDescriptionIndex,receivedMetadataIndex);
//            receiveFiniteVolumesMetadataInMergeWithNeigbour(destCellDescriptionIndex,receivedMetadataIndex);
          }
          // Clean up
          MetadataHeap::getInstance().deleteData(receivedMetadataIndex);
        }
      enddforx
    enddforx
  }

  logTraceOut("mergeWithNeighbour(...)");
}


void exahype::mappings::MarkingForAugmentation::receiveADERDGMetadataInMergeWithNeigbour(
    const int destCellDescriptionIndex,
    const int receivedMetadataIndex) {
  exahype::MetadataHeap::HeapEntries::const_iterator metadataIterator=
      MetadataHeap::getInstance().getData(receivedMetadataIndex).begin();
  const int nADERDG = metadataIterator->getU(); ++metadataIterator;

  if (nADERDG > 0) {
    logDebug("mergeWithNeighbour(...)","nADERDG: " << nADERDG);
    while (metadataIterator!=MetadataHeap::getInstance().getData(receivedMetadataIndex).begin()+1+2*nADERDG) {
      const int solverNumber = metadataIterator->getU(); ++metadataIterator;
      const int typeAsInt    = metadataIterator->getU(); ++metadataIterator;
      exahype::records::ADERDGCellDescription::Type type =
            static_cast<exahype::records::ADERDGCellDescription::Type>(typeAsInt);

      logDebug("mergeWithNeighbour(...)","solverNumber: " << solverNumber);
      logDebug("mergeWithNeighbour(...)","typeAsInt: "    << typeAsInt);

      for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(destCellDescriptionIndex)) {
        if (p.getSolverNumber()==solverNumber) {
          switch(p.getType()) {
          case exahype::records::ADERDGCellDescription::Cell:
            if (p.getRefinementEvent()==exahype::records::ADERDGCellDescription::None &&
                (type==exahype::records::ADERDGCellDescription::Ancestor ||
                 type==exahype::records::ADERDGCellDescription::EmptyAncestor)) {
              p.setRefinementEvent(exahype::records::ADERDGCellDescription::AugmentingRequested);
            }
            break;
          case exahype::records::ADERDGCellDescription::Descendant:
          case exahype::records::ADERDGCellDescription::EmptyDescendant:
            // TODO(Dominic): Add to docu what we do here.
            if (type==exahype::records::ADERDGCellDescription::Cell) {
              p.setHelperCellNeedsToStoreFaceData(true);
            }

            // 2. Request further augmentation if necessary (this might get reseted if the traversal
            // is able to descend and finds existing descendants).
            if (p.getRefinementEvent()==exahype::records::ADERDGCellDescription::None &&
                (type==exahype::records::ADERDGCellDescription::Ancestor ||
                type==exahype::records::ADERDGCellDescription::EmptyAncestor)) {
              p.setRefinementEvent(exahype::records::ADERDGCellDescription::AugmentingRequested);
            }
            break;
          case exahype::records::ADERDGCellDescription::Ancestor:
          case exahype::records::ADERDGCellDescription::EmptyAncestor:
            // TODO(Dominic): Add to docu what we do here.
            if (type==exahype::records::ADERDGCellDescription::Cell) {
              p.setHelperCellNeedsToStoreFaceData(true);
            }
            break;
          default:
            assertionMsg(false,"Should never be entered!");
            break;
          }
        }
      }
    }
  }
}

void exahype::mappings::MarkingForAugmentation::receiveFiniteVolumesMetadataInMergeWithNeigbour(
    const int destCellDescriptionIndex,
    const int receivedMetadataIndex) {
// TODO(Dominic): Implement.
  assertionMsg(false,"Dominic please implement!");
}



//
// All functions below are nop.
//




void exahype::mappings::MarkingForAugmentation::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
 // do nothing
}

bool exahype::mappings::MarkingForAugmentation::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  // cell-based
  return true;
}

void exahype::mappings::MarkingForAugmentation::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::mergeWithMaster(
    const exahype::Cell& workerGridCell,
    exahype::Vertex* const workerGridVertices,
    const peano::grid::VertexEnumerator& workerEnumerator,
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker, const exahype::State& workerState,
    exahype::State& masterState) {
}

// cell-based
void exahype::mappings::MarkingForAugmentation::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

// cell-based
void exahype::mappings::MarkingForAugmentation::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

// cell-based
void exahype::mappings::MarkingForAugmentation::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::MarkingForAugmentation::MarkingForAugmentation() {
  // do nothing
}

exahype::mappings::MarkingForAugmentation::~MarkingForAugmentation() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::MarkingForAugmentation::MarkingForAugmentation(
    const MarkingForAugmentation& masterThread) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::mergeWithWorkerThread(
    const MarkingForAugmentation& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::MarkingForAugmentation::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::endIteration(
    exahype::State& solverState) {
  // do nothing
}
