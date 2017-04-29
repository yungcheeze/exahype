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
 
#include "exahype/mappings/MeshRefinement.h"

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "tarch/la/VectorScalarOperations.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#ifdef Parallel
bool exahype::mappings::MeshRefinement::FirstIteration = true;
#endif

exahype::mappings::MeshRefinement::RefinementMode
exahype::mappings::MeshRefinement::Mode = exahype::mappings::MeshRefinement::RefinementMode::APriori;

tarch::logging::Log exahype::mappings::MeshRefinement::_log("exahype::mappings::MeshRefinement");

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::MeshRefinement::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::MeshRefinement::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      #ifdef Parallel
      peano::MappingSpecification::WholeTree,
      #else
      peano::MappingSpecification::Nop,
      #endif
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification
exahype::mappings::MeshRefinement::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification
exahype::mappings::MeshRefinement::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial,true);
}
peano::MappingSpecification
exahype::mappings::MeshRefinement::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial,true);
}
peano::MappingSpecification
exahype::mappings::MeshRefinement::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::MeshRefinement::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::MeshRefinement::MeshRefinement(const MeshRefinement& masterThread):
  _localState(masterThread._localState) {
}
#endif

void exahype::mappings::MeshRefinement::beginIteration(
  exahype::State& solverState
) {
  exahype::solvers::ADERDGSolver::Heap::getInstance().setName("ADERDGCellDescriptionHeap");
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().setName("FiniteVolumesCellDescriptionHeap");
  DataHeap::getInstance().setName("DataHeap");

  // TODO(Dominic):
//  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
//    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
//
//    solver->zeroTimeStepSizes();
//    assertion1(!solver->getNextGridUpdateRequested(),solver->toString());
//  } // Dead code elimination will get rid of this loop in Asserts and Debug mode.

  #ifdef Parallel
  exahype::solvers::ADERDGSolver::Heap::getInstance().finishedToSendSynchronousData();
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().finishedToSendSynchronousData();
  DataHeap::getInstance().finishedToSendSynchronousData();

  MetadataHeap::getInstance().setName("MetadataHeap");
  MetadataHeap::getInstance().finishedToSendSynchronousData();
  if (! MetadataHeap::getInstance().validateThatIncomingJoinBuffersAreEmpty() ) {
      exit(-1);
  }

  exahype::solvers::ADERDGSolver::Heap::getInstance().startToSendSynchronousData();
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().startToSendSynchronousData();
  DataHeap::getInstance().startToSendSynchronousData();
  MetadataHeap::getInstance().startToSendSynchronousData();
  #endif

  _localState = solverState;
}

void exahype::mappings::MeshRefinement::endIteration(exahype::State& solverState) {
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

//      solver->zeroTimeStepSizes();
      solver->setNextGridUpdateRequested();
  }

  #ifdef Parallel
  exahype::mappings::MeshRefinement::FirstIteration = false;
  #endif
}

void exahype::mappings::MeshRefinement::refineVertexIfNecessary(
  exahype::Vertex&                              fineGridVertex,
  const tarch::la::Vector<DIMENSIONS, double>&  fineGridH,
  int                                           fineGridLevel,
  bool                                          isCalledByCreationalEvent
) const {
  bool refine = false;
  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    refine |= tarch::la::allGreater(fineGridH,p->getMaximumMeshSize());
  }

  if (
    refine
    &&
    fineGridVertex.getRefinementControl()==Vertex::Records::Unrefined
  ) {
    switch ( _localState.mayRefine(isCalledByCreationalEvent,fineGridLevel) ) {
      case State::RefinementAnswer::DontRefineYet:
        break;
      case State::RefinementAnswer::Refine:
        fineGridVertex.refine();
        break;
      case State::RefinementAnswer::EnforceRefinement:
        fineGridVertex.enforceRefine();
        break;
    }
  }
}


void exahype::mappings::MeshRefinement::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  refineVertexIfNecessary(fineGridVertex,fineGridH,coarseGridVerticesEnumerator.getLevel()+1,false);
}


void exahype::mappings::MeshRefinement::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("createBoundaryVertex(...)", fineGridVertex,
                           fineGridX, fineGridH,
                           coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  refineVertexIfNecessary(fineGridVertex,fineGridH,coarseGridVerticesEnumerator.getLevel()+1,true);

  logTraceOutWith1Argument("createBoundaryVertex(...)", fineGridVertex);
}


void exahype::mappings::MeshRefinement::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("createInnerVertex(...)", fineGridVertex, fineGridX,
                           fineGridH, coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  refineVertexIfNecessary(fineGridVertex,fineGridH,coarseGridVerticesEnumerator.getLevel()+1,true);

  logTraceOutWith1Argument("createInnerVertex(...)", fineGridVertex);
}


void exahype::mappings::MeshRefinement::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("createCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);
  fineGridCell.getCellData().setCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);

  logTraceOutWith1Argument("createCell(...)", fineGridCell);
}

void exahype::mappings::MeshRefinement::mergeNeighboursLimiterStatus(exahype::Vertex& fineGridVertex) {
  dfor2(pos1)
    dfor2(pos2)
      if (fineGridVertex.hasToMergeNeighbours(pos1,pos2)) { // Assumes that we have to valid indices
        auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().
            parallelise(solvers::RegisteredSolvers.size(), peano::datatraversal::autotuning::MethodTrace::UserDefined15);
        pfor(solverNumber, 0, static_cast<int>(solvers::RegisteredSolvers.size()),grainSize.getGrainSize())
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

        if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
          const int cellDescriptionsIndex1 = fineGridVertex.getCellDescriptionsIndex()[pos1Scalar];
          const int cellDescriptionsIndex2 = fineGridVertex.getCellDescriptionsIndex()[pos2Scalar];
          const int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
          const int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
          if (element2>=0 && element1>=0) {
            static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
                mergeNeighboursLimiterStatus(
                    cellDescriptionsIndex1,element1,
                    cellDescriptionsIndex2,element2,
                    pos1,pos2);
          }
        }
        endpfor
        grainSize.parallelSectionHasTerminated();

        fineGridVertex.setMergePerformed(pos1,pos2,true);
      }
    enddforx
  enddforx
}

void exahype::mappings::MeshRefinement::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  dfor2(pos1)
    dfor2(pos2)
      if (fineGridVertex.hasToMergeNeighbours(pos1,pos2)) { // Assumes that we have to valid indices
        auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().
            parallelise(solvers::RegisteredSolvers.size(), peano::datatraversal::autotuning::MethodTrace::UserDefined15);
        pfor(solverNumber, 0, static_cast<int>(solvers::RegisteredSolvers.size()),grainSize.getGrainSize())
          auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

          if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
            const int cellDescriptionsIndex1 = fineGridVertex.getCellDescriptionsIndex()[pos1Scalar];
            const int cellDescriptionsIndex2 = fineGridVertex.getCellDescriptionsIndex()[pos2Scalar];
            const int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
            const int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
            if (element2>=0 && element1>=0) {
              static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
                  mergeNeighboursLimiterStatus(
                  cellDescriptionsIndex1,element1,
                  cellDescriptionsIndex2,element2,
                  pos1,pos2);
            }
          }
        endpfor
        grainSize.parallelSectionHasTerminated();

        fineGridVertex.setMergePerformed(pos1,pos2,true);
      }
    enddforx
  enddforx
}

void exahype::mappings::MeshRefinement::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);
  bool refineFineGridCell = false;

  for (unsigned int solverNumber=0; solverNumber<exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    refineFineGridCell |=
        solver->markForRefinement( // TODO(Dominic): Consider the maximum refinement depth here
            fineGridCell,
            fineGridVertices,
            fineGridVerticesEnumerator,
            coarseGridCell,
            coarseGridVertices,
            coarseGridVerticesEnumerator,
            fineGridPositionOfCell,
            MeshRefinement::Mode==RefinementMode::Initial,
            solverNumber);

    refineFineGridCell |=
        solver->updateStateInEnterCell(
            fineGridCell,
            fineGridVertices,
            fineGridVerticesEnumerator,
            coarseGridCell,
            coarseGridVertices,
            coarseGridVerticesEnumerator,
            fineGridPositionOfCell,
            MeshRefinement::Mode==RefinementMode::Initial,
            solverNumber);

    const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
    if (element!=exahype::solvers::Solver::NotFound) {
      solver->prepareNextNeighbourMerging(
          fineGridCell.getCellDescriptionsIndex(),element,
          fineGridVertices,fineGridVerticesEnumerator);
//      // TODO(Dominic):
//      solver->zeroTimeStepSizes(fineGridCell.getCellDescriptionsIndex(),element);
    }
  }

  // Refine all adjacent vertices if necessary and possible.
  //    if (refineFineGridCell && _state.refineInitialGridInTouchVertexLastTime()) {
  if (refineFineGridCell) {
    dfor2(v)
      if (fineGridVertices[ fineGridVerticesEnumerator(v) ].getRefinementControl()==
              exahype::Vertex::Records::RefinementControl::Unrefined) {
        fineGridVertices[ fineGridVerticesEnumerator(v) ].refine();
      }
    enddforx
  }

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::MeshRefinement::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("leaveCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);
  bool eraseFineGridCell = true;

  for (unsigned int solverNumber=0; solverNumber<exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    eraseFineGridCell &=
        solver->updateStateInLeaveCell(
            fineGridCell,
            fineGridVertices,
            fineGridVerticesEnumerator,
            coarseGridCell,
            coarseGridVertices,
            coarseGridVerticesEnumerator,
            fineGridPositionOfCell,
            solverNumber);

    const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
    if (element!=exahype::solvers::Solver::NotFound) {
      bool isStable = solver->attainedStableState(
          fineGridCell,
          fineGridVertices,
          fineGridVerticesEnumerator,
          solverNumber);
      solver->updateNextGridUpdateRequested(!isStable);
    }
  }

  // We assume that the solvers have all removed
  // their data from the heap arrays associated with this
  // cell.
  if (eraseFineGridCell) {
    fineGridCell.shutdownMetaData();

    dfor2(v)
    if (coarseGridVertices[ coarseGridVerticesEnumerator(v) ].getRefinementControl()==
        exahype::Vertex::Records::RefinementControl::Refined
    ) {
      coarseGridVertices[ coarseGridVerticesEnumerator(v) ].erase();
    }
    enddforx
  }

  logTraceOutWith1Argument("leaveCell(...)", fineGridCell);
}

#ifdef Parallel
bool exahype::mappings::MeshRefinement::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  return true;
}

void exahype::mappings::MeshRefinement::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments("mergeWithNeighbour(...)", vertex, neighbour,
                           fromRank, fineGridX, fineGridH, level);

  if (exahype::mappings::MeshRefinement::FirstIteration) {
    return;
  }
  if (tarch::la::allGreater(fineGridH,exahype::solvers::Solver::getCoarsestMeshSizeOfAllSolvers())) {
    return;
  }

  // Keep this for the grid setup
  #if !defined(PeriodicBC)
    if (vertex.isBoundary()) return;
  #endif

  dfor2(myDest)
    dfor2(mySrc)
      tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
      tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;
      int destScalar = TWO_POWER_D - myDestScalar - 1;

      if (vertex.hasToReceiveMetadata(src,dest,fromRank)) {
        logDebug("mergeWithNeighbour(...)","[pre] rec. from rank "<<fromRank<<", x:"<<
                 fineGridX.toString() << ", level=" <<level << ", vertex.adjacentRanks: "
                 << vertex.getAdjacentRanks());

        const int receivedMetadataIndex = MetadataHeap::getInstance().
            createData(0,exahype::MetadataPerSolver*exahype::solvers::RegisteredSolvers.size());
        assertion(MetadataHeap::getInstance().getData(receivedMetadataIndex).empty());
        MetadataHeap::getInstance().receiveData(
            receivedMetadataIndex,
            fromRank, fineGridX, level,
            peano::heap::MessageType::NeighbourCommunication);
        MetadataHeap::HeapEntries& receivedMetadata = MetadataHeap::getInstance().
            getData(receivedMetadataIndex);

        // Work with the neighbour cell type
        for(unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0;) {
          auto* solver = solvers::RegisteredSolvers[solverNumber];

          const int offset  = exahype::MetadataPerSolver*solverNumber;
          if (receivedMetadata[offset].getU()!=exahype::InvalidMetadataEntry) {
            const int element = solver->tryGetElement(
                vertex.getCellDescriptionsIndex()[destScalar],solverNumber);
            if (element!=exahype::solvers::Solver::NotFound) {

              exahype::MetadataHeap::HeapEntries metadataPortion(
                        receivedMetadata.begin()+offset,
                        receivedMetadata.begin()+offset+exahype::MetadataPerSolver);

              solver->mergeWithNeighbourMetadata(
                  metadataPortion,
                  src, dest,
                  vertex.getCellDescriptionsIndex()[destScalar],element);
            }
          }

          logDebug("mergeWithNeighbour(...)","solverNumber: " << solverNumber);
          logDebug("mergeWithNeighbour(...)","neighbourTypeAsInt: "
                   << receivedMetadata[solverNumber].getU());
        }

        // Clean up
        MetadataHeap::getInstance().deleteData(receivedMetadataIndex);
      }
    enddforx
  enddforx

  logTraceOut("mergeWithNeighbour(...)");
}

void exahype::mappings::MeshRefinement::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  logTraceInWith5Arguments("prepareSendToNeighbour(...)", vertex,
                           toRank, x, h, level);

  if (tarch::la::allGreater(h,exahype::solvers::Solver::getCoarsestMeshSizeOfAllSolvers())) {
    return;
  }

  // Keep this for the grid setup
  #if !defined(PeriodicBC)
  if (vertex.isBoundary()) return;
  #endif

  tarch::la::Vector<TWO_POWER_D, int> adjacentADERDGCellDescriptionsIndices =
      vertex.getCellDescriptionsIndex();

  dfor2(dest)
    dfor2(src)
      if (vertex.hasToSendMetadata(src,dest,toRank)) {
        const int srcCellDescriptionIndex = adjacentADERDGCellDescriptionsIndices(srcScalar);
        if (exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionIndex)) {
          logDebug("prepareSendToNeighbour(...)","[metadata] sent to rank "<<toRank<<", x:"<<
                   x.toString() << ", level=" <<level << ", vertex.adjacentRanks: " << vertex.getAdjacentRanks()
                   << ", src forking: " << State::isForkingRank(vertex.getAdjacentRanks()(srcScalar)));

          exahype::Vertex::sendEncodedMetadata(
              toRank,srcCellDescriptionIndex,peano::heap::MessageType::NeighbourCommunication,x,level);
        } else {
          logDebug("prepareSendToNeighbour(...)","[empty] sent to rank "<<toRank<<", x:"<<
                   x.toString() << ", level=" <<level << ", vertex.adjacentRanks: " << vertex.getAdjacentRanks()
                   << ", src forking: " << State::isForkingRank(vertex.getAdjacentRanks()(srcScalar)));

          exahype::Vertex::sendEncodedMetadataSequenceWithInvalidEntries(
              toRank,peano::heap::MessageType::NeighbourCommunication,x,level);
        }
      }
    enddforx
  enddforx

  logTraceOut("prepareSendToNeighbour(...)");
}

void exahype::mappings::MeshRefinement::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  if (localCell.isInside() && localCell.isInitialised()) {
    exahype::solvers::ADERDGSolver::sendCellDescriptions(toRank,localCell.getCellDescriptionsIndex(),
        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);
//    exahype::solvers::FiniteVolumesSolver::sendCellDescriptions(toRank,localCell.getCellDescriptionsIndex(),
//        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);

    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      int element = solver->tryGetElement(localCell.getCellDescriptionsIndex(),solverNumber);
      if(element!=exahype::solvers::Solver::NotFound) {
        solver->sendDataToWorkerOrMasterDueToForkOrJoin(
            toRank,localCell.getCellDescriptionsIndex(),element,cellCentre,level);
      } else {
        solver->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(toRank,cellCentre,level);
      }
    }
  } else if (localCell.isInside() && !localCell.isInitialised()){
    exahype::solvers::ADERDGSolver::sendEmptyCellDescriptions(toRank,
        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);
    exahype::solvers::FiniteVolumesSolver::sendEmptyCellDescriptions(toRank,
        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);

    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
       auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      solver->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(toRank,cellCentre,level);
    }
  }
}

// TODO(Dominic): How to deal with cell descriptions index that
// is copied from the remote rank but is a valid index on the local
// remote rank? Currently use geometryInfoDoesMatch!
void exahype::mappings::MeshRefinement::mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  if (localCell.isInside()) {
    if (!geometryInfoDoesMatch(localCell.getCellDescriptionsIndex(),cellCentre,cellSize,level)) {
      localCell.setCellDescriptionsIndex(
          multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
    }

    exahype::solvers::ADERDGSolver::mergeCellDescriptionsWithRemoteData(
        fromRank,localCell,
        peano::heap::MessageType::ForkOrJoinCommunication,
        cellCentre,level);
//    exahype::solvers::FiniteVolumesSolver::mergeCellDescriptionsWithRemoteData(
//        fromRank,localCell,
//        peano::heap::MessageType::ForkOrJoinCommunication,
//        cellCentre,level);

    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      int element = solver->tryGetElement(localCell.getCellDescriptionsIndex(),solverNumber);
      if (element!=exahype::solvers::Solver::NotFound) {
        solver->mergeWithWorkerOrMasterDataDueToForkOrJoin(
            fromRank,localCell.getCellDescriptionsIndex(),element,cellCentre,level);
      } else {
        solver->dropWorkerOrMasterDataDueToForkOrJoin(fromRank,cellCentre,level);
      }
    }
  }
}

bool exahype::mappings::MeshRefinement::geometryInfoDoesMatch(
    const int cellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>& cellSize,
    const int level) {
  int solverNumber=0;
  for (auto& solver : exahype::solvers::RegisteredSolvers) {
    int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);
    if (element!=exahype::solvers::Solver::NotFound) {
      if (solver->getType()==exahype::solvers::Solver::Type::ADERDG) {
        exahype::solvers::ADERDGSolver::CellDescription& cellDescription =
            exahype::solvers::ADERDGSolver::getCellDescription(
                cellDescriptionsIndex,element);

        if (!tarch::la::equals(
            cellCentre,cellDescription.getOffset()+0.5*cellDescription.getSize()) ||
            cellDescription.getLevel()!=level) {
          return false;
        }
      } else if (solver->getType()==exahype::solvers::Solver::Type::FiniteVolumes) {
        exahype::solvers::FiniteVolumesSolver::CellDescription& cellDescription =
            exahype::solvers::FiniteVolumesSolver::getCellDescription(
                cellDescriptionsIndex,element);

        if (!tarch::la::equals(
            cellCentre,cellDescription.getOffset()+0.5*cellDescription.getSize()) ||
            cellDescription.getLevel()!=level) {
          return false;
        }
      }
    }
    ++solverNumber;
  }

  return true;
}


void exahype::mappings::MeshRefinement::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      solver->sendGridUpdateFlagsToMaster(
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          verticesEnumerator.getCellCenter(),
          verticesEnumerator.getLevel());
  }
}

void exahype::mappings::MeshRefinement::mergeWithMaster(
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
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    solver->mergeWithWorkerGridUpdateFlags(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());
  }
}


//
// All methods below are nop,
//
// ==================================



void exahype::mappings::MeshRefinement::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MeshRefinement::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MeshRefinement::receiveDataFromMaster(
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

void exahype::mappings::MeshRefinement::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::MeshRefinement::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::MeshRefinement::MeshRefinement() {
  // do nothing
}

exahype::mappings::MeshRefinement::~MeshRefinement() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
void exahype::mappings::MeshRefinement::mergeWithWorkerThread(
    const MeshRefinement& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::MeshRefinement::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MeshRefinement::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MeshRefinement::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MeshRefinement::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}


void exahype::mappings::MeshRefinement::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::MeshRefinement::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
