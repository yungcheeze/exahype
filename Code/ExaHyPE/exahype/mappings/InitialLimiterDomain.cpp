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

#include "exahype/mappings/InitialLimiterDomain.h"

#include <cmath>

#include "peano/utils/Globals.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "tarch/multicore/Loop.h"

#include "peano/utils/Loop.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

peano::CommunicationSpecification
exahype::mappings::InitialLimiterDomain::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
      MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
      MaskOutWorkerMasterDataAndStateExchange,
      true);
}

// Everything below is nop.
peano::MappingSpecification
exahype::mappings::InitialLimiterDomain::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces);
}

peano::MappingSpecification
exahype::mappings::InitialLimiterDomain::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial);
}

peano::MappingSpecification
exahype::mappings::InitialLimiterDomain::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

peano::MappingSpecification
exahype::mappings::InitialLimiterDomain::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}

peano::MappingSpecification
exahype::mappings::InitialLimiterDomain::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

peano::MappingSpecification
exahype::mappings::InitialLimiterDomain::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::InitialLimiterDomain::_log(
    "exahype::mappings::InitialLimiterDomain");

void exahype::mappings::InitialLimiterDomain::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (fineGridCell.isInitialised()) {
    const int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
    auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfSolvers, peano::datatraversal::autotuning::MethodTrace::UserDefined7);
    pfor(i, 0, numberOfSolvers, grainSize.getGrainSize())
    auto solver = exahype::solvers::RegisteredSolvers[i];

    const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),i);
    if (element!=exahype::solvers::Solver::NotFound
        && solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
      auto limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);

      limitingADERDGSolver->synchroniseTimeStepping(fineGridCell.getCellDescriptionsIndex(),element);
      limitingADERDGSolver->rollbackToPreviousTimeStep(fineGridCell.getCellDescriptionsIndex(),element);
      limitingADERDGSolver->initialiseLimiter(fineGridCell.getCellDescriptionsIndex(),element,
                                              fineGridCell,fineGridVertices,fineGridVerticesEnumerator);
    }
    endpfor
    grainSize.parallelSectionHasTerminated();
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::InitialLimiterDomain::touchVertexFirstTime(
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
  dfor2(pos1)
    dfor2(pos2)
      if (fineGridVertex.hasToMergeNeighbours(pos1,pos2)) { // Assumes that we have to valid indices // TODO(Dominic): Probably have to consider Voronoi neighbours later on when we use high order schemes
        auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().
            parallelise(solvers::RegisteredSolvers.size(), peano::datatraversal::autotuning::MethodTrace::UserDefined2);
        pfor(solverNumber, 0, static_cast<int>(solvers::RegisteredSolvers.size()),grainSize.getGrainSize())
          auto solver = exahype::solvers::RegisteredSolvers[solverNumber];

          if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
            const int cellDescriptionsIndex1 = fineGridVertex.getCellDescriptionsIndex()[pos1Scalar];
            const int cellDescriptionsIndex2 = fineGridVertex.getCellDescriptionsIndex()[pos2Scalar];
            const int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
            const int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
            if (element2>=0 && element1>=0) {
              auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
              limitingADERDGSolver->mergeLimiterStatusOfNeighbours(
                  cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
            }
          }

          #ifdef Debug // TODO(Dominic)
          _interiorFaceMerges++;
          #endif
        endpfor
        grainSize.parallelSectionHasTerminated();
        fineGridVertex.setMergePerformed(pos1,pos2,true);
      }
    enddforx
  enddforx
}

#ifdef Parallel
void exahype::mappings::InitialLimiterDomain::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
// TODO(Dominic): Merge my limiter status with the neighbour's
//  dfor2(myDest)
//    dfor2(mySrc)
//      tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
//      tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;
//
//      int destScalar = TWO_POWER_D - myDestScalar - 1;
//      int srcScalar  = TWO_POWER_D - mySrcScalar  - 1;
//
//      if (vertex.hasToReceiveMetadata(src,dest,fromRank)) {
//        int receivedMetadataIndex = MetadataHeap::getInstance().
//            createData(0,exahype::solvers::RegisteredSolvers.size());
//        MetadataHeap::getInstance().receiveData(
//            receivedMetadataIndex,
//            fromRank, fineGridX, level,
//            peano::heap::MessageType::NeighbourCommunication);
//        exahype::MetadataHeap::HeapEntries& receivedMetadata = MetadataHeap::getInstance().getData(receivedMetadataIndex);
//        assertion(receivedMetadata.size()==solvers::RegisteredSolvers.size());
//
//        if(vertex.hasToMergeWithNeighbourData(src,dest)) {
//          //
//
//          vertex.setFaceDataExchangeCountersOfDestination(src,dest,TWO_POWER_D);
//          vertex.setMergePerformed(src,dest,true);
//        }
//      }
//    enddforx
//  enddforx
}

void exahype::mappings::InitialLimiterDomain::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // TODO(Dominic): Send my limiter status to neighbour
  // TODO(Dominic): Discussion: Send status always or use counters?
//  dfor2(dest)
//    dfor2(src)
//    if (vertex.hasToSendMetadata(src,dest,toRank)) {
//      vertex.tryDecrementFaceDataExchangeCountersOfSource(src,dest);
//      if (vertex.hasToSendDataToNeighbour(src,dest)) {
//        sendSolverDataToNeighbour(
//            toRank,src,dest,
//            vertex.getCellDescriptionsIndex()[srcScalar],
//            vertex.getCellDescriptionsIndex()[destScalar],
//            x,level);
//      } else {
//        sendEmptySolverDataToNeighbour(toRank,src,dest,x,level);
//      }
//    }
//    enddforx
//  enddforx
}

//
// Below all methods are nop.
//
//===================================

void exahype::mappings::InitialLimiterDomain::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::InitialLimiterDomain::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  // do nothing
  return false;
}

void exahype::mappings::InitialLimiterDomain::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::mergeWithMaster(
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
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::receiveDataFromMaster(
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

void exahype::mappings::InitialLimiterDomain::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::InitialLimiterDomain::InitialLimiterDomain() {
  // do nothing
}

exahype::mappings::InitialLimiterDomain::~InitialLimiterDomain() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::InitialLimiterDomain::InitialLimiterDomain(
    const InitialLimiterDomain& masterThread) {
  // do nothing
}
void exahype::mappings::InitialLimiterDomain::mergeWithWorkerThread(
    const InitialLimiterDomain& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::InitialLimiterDomain::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::beginIteration(
    exahype::State& solverState) {
}

void exahype::mappings::InitialLimiterDomain::endIteration(
    exahype::State& solverState) {
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
      auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
      limitingADERDGSolver->rollbackToPreviousTimeStep();
    }
  }
}

void exahype::mappings::InitialLimiterDomain::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::InitialLimiterDomain::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
