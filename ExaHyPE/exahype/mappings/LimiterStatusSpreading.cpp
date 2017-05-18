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
 
#include "exahype/mappings/LimiterStatusSpreading.h"

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "tarch/la/VectorScalarOperations.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#ifdef Parallel
bool exahype::mappings::LimiterStatusSpreading::FirstIteration = true;
#endif

tarch::logging::Log exahype::mappings::LimiterStatusSpreading::_log("exahype::mappings::LimiterStatusSpreading");

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::LimiterStatusSpreading::communicationSpecification() const {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::touchVertexFirstTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::enterCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::Serial,true);
}
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::LimiterStatusSpreading::LimiterStatusSpreading(const LimiterStatusSpreading& masterThread) {
}
#endif

void exahype::mappings::LimiterStatusSpreading::beginIteration(
  exahype::State& solverState
) {
  // We memorise the previous request per solver
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

//      solver->zeroTimeStepSizes();
      solver->updateNextMeshUpdateRequest(solver->getMeshUpdateRequest());
  }

  // TODO(Dominic): Prepare variables for multithreading

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
}

void exahype::mappings::LimiterStatusSpreading::endIteration(exahype::State& solverState) {
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

//      solver->zeroTimeStepSizes();
      solver->setNextMeshUpdateRequest();
  }

  #ifdef Parallel
  exahype::mappings::LimiterStatusSpreading::FirstIteration = false;
  #endif
}

void exahype::mappings::LimiterStatusSpreading::mergeNeighboursLimiterStatus(exahype::Vertex& fineGridVertex) {
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

void exahype::mappings::LimiterStatusSpreading::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("touchVertexFirstTime(...)", fineGridVertex,
                           fineGridX, fineGridH,
                           coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  mergeNeighboursLimiterStatus(fineGridVertex);

  logTraceOutWith1Argument("touchVertexFirstTime(...)", fineGridVertex);
}

void exahype::mappings::LimiterStatusSpreading::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  for (unsigned int solverNumber=0; solverNumber<exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    const int element =
        solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
    if (element!=exahype::solvers::Solver::NotFound) {
      if (
          solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
          &&
          static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
          !=exahype::solvers::LimiterDomainChange::Regular
      ) {
        auto* limitingADERDG = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);

        limitingADERDG->updateLimiterStatus(fineGridCell.getCellDescriptionsIndex(),element);
        limitingADERDG->deallocateLimiterPatchOnHelperCell(fineGridCell.getCellDescriptionsIndex(),element);
        limitingADERDG->ensureRequiredLimiterPatchIsAllocated(fineGridCell.getCellDescriptionsIndex(),element);

        bool meshUpdateRequest =
            limitingADERDG->
              evaluateLimiterStatusBasedRefinementCriterion(
                  fineGridCell.getCellDescriptionsIndex(),element);
        // TODO(Dominic): Enable multithreading for this; have value per solver; reduce in endIteration
        limitingADERDG->updateNextMeshUpdateRequest(meshUpdateRequest);
      }

      solver->prepareNextNeighbourMerging(
          fineGridCell.getCellDescriptionsIndex(),element,
          fineGridVertices,fineGridVerticesEnumerator);
    }
  }

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::LimiterStatusSpreading::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
bool exahype::mappings::LimiterStatusSpreading::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  return true;
}

void exahype::mappings::LimiterStatusSpreading::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments("mergeWithNeighbour(...)", vertex, neighbour,
                           fromRank, fineGridX, fineGridH, level);

  if (exahype::mappings::LimiterStatusSpreading::FirstIteration) {
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

              // TODO(Dominic): Refinement event is set directly according to
              // solver metadata. This should not happen.
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

void exahype::mappings::LimiterStatusSpreading::prepareSendToNeighbour(
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

void exahype::mappings::LimiterStatusSpreading::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}


void exahype::mappings::LimiterStatusSpreading::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      // TODO(Dominic): Restrict limiter status too
      solver->sendMeshUpdateFlagsToMaster(
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          verticesEnumerator.getCellCenter(),
          verticesEnumerator.getLevel());
  }
}

void exahype::mappings::LimiterStatusSpreading::mergeWithMaster(
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

    // TODO(Dominic): Merge restricted limiter status
    solver->mergeWithWorkerMeshUpdateFlags(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());
  }
}


//
// All methods below are nop,
//
// ==================================



void exahype::mappings::LimiterStatusSpreading::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::receiveDataFromMaster(
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

void exahype::mappings::LimiterStatusSpreading::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::LimiterStatusSpreading::LimiterStatusSpreading() {
  // do nothing
}

exahype::mappings::LimiterStatusSpreading::~LimiterStatusSpreading() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
void exahype::mappings::LimiterStatusSpreading::mergeWithWorkerThread(
    const LimiterStatusSpreading& workerThread) {
  // do nothing
}
#endif



void exahype::mappings::LimiterStatusSpreading::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}


void exahype::mappings::LimiterStatusSpreading::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}


void exahype::mappings::LimiterStatusSpreading::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}


void exahype::mappings::LimiterStatusSpreading::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
 // do nothing
}

void exahype::mappings::LimiterStatusSpreading::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}


void exahype::mappings::LimiterStatusSpreading::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
