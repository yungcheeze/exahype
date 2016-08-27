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
 
#include "exahype/mappings/SpaceTimePredictor.h"
#include "exahype/solvers/ADERDGSolver.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/utils/Globals.h"

#include "tarch/multicore/Loop.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

peano::CommunicationSpecification
exahype::mappings::SpaceTimePredictor::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}

peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

// The remainder specs all are nop
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::SpaceTimePredictor::_log(
    "exahype::mappings::SpaceTimePredictor");

int exahype::mappings::SpaceTimePredictor::_mpiTag =
    tarch::parallel::Node::reserveFreeTag(
        "exahype::mappings::SpaceTimePredictor");

#if defined(SharedMemoryParallelisation)
exahype::mappings::SpaceTimePredictor::SpaceTimePredictor(
    const SpaceTimePredictor& masterThread) {
  // TODO(Dominic): Merge from other branch.
}
void exahype::mappings::SpaceTimePredictor::mergeWithWorkerThread(
    const SpaceTimePredictor& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::SpaceTimePredictor::enterCell(
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
    const int numberOfADERDGCellDescriptions = static_cast<int>(
        ADERDGCellDescriptionHeap::getInstance()
            .getData(fineGridCell.getCellDescriptionsIndex())
            .size());

    // please use a different UserDefined per mapping/event
    const peano::datatraversal::autotuning::MethodTrace methodTrace =
        peano::datatraversal::autotuning::UserDefined7;
    const int grainSize =
        peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
            numberOfADERDGCellDescriptions, methodTrace);
    pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
      auto& p = fineGridCell.getADERDGCellDescription(i);

      exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
          exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);

      // Reset helper variables
      for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
        p.setRiemannSolvePerformed(faceIndex,false);

        // TODO(Dominic): Add to docu.
        #ifdef Parallel
        int listingsOfRemoteRank = countListingsOfRemoteRankAtFace(faceIndex,fineGridVertices,fineGridVerticesEnumerator);
        if (listingsOfRemoteRank==0) {
          listingsOfRemoteRank = TWO_POWER_D;
        }
        p.setFaceDataExchangeCounter(faceIndex,listingsOfRemoteRank);
        #endif
      }

      // space-time DoF (basisSize**(DIMENSIONS+1))
      double* lQi = 0;
      double* lFi = 0;

      // volume DoF (basisSize**(DIMENSIONS))
      double* luh = 0;
      double* lduh = 0;
      double* lQhi = 0;
      double* lFhi = 0;

      // face DoF (basisSize**(DIMENSIONS-1))
      double* lQhbnd = 0;
      double* lFhbnd = 0;

      switch (p.getType()) {
        case exahype::records::ADERDGCellDescription::Cell:
          assertion1(p.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,p.toString());
          assertion1(DataHeap::getInstance().isValidIndex(p.getSpaceTimePredictor()),p.toString());
          assertion1(DataHeap::getInstance().isValidIndex(p.getSpaceTimeVolumeFlux()),p.toString());
          assertion1(DataHeap::getInstance().isValidIndex(p.getSolution()),p.toString());
          assertion1(DataHeap::getInstance().isValidIndex(p.getUpdate()),p.toString());
          assertion1(DataHeap::getInstance().isValidIndex(p.getPredictor()),p.toString());
          assertion1(DataHeap::getInstance().isValidIndex(p.getVolumeFlux()),p.toString());
          assertion1(DataHeap::getInstance().isValidIndex(p.getExtrapolatedPredictor()),p.toString());
          assertion1(DataHeap::getInstance().isValidIndex(p.getFluctuation()),p.toString());

          lQi = DataHeap::getInstance().getData(p.getSpaceTimePredictor()).data();
          lFi = DataHeap::getInstance().getData(p.getSpaceTimeVolumeFlux()).data();

          luh  = DataHeap::getInstance().getData(p.getSolution()).data();
          lduh = DataHeap::getInstance().getData(p.getUpdate()).data();
          lQhi = DataHeap::getInstance().getData(p.getPredictor()).data();
          lFhi = DataHeap::getInstance().getData(p.getVolumeFlux()).data();

          lQhbnd = DataHeap::getInstance().getData(p.getExtrapolatedPredictor()).data();
          lFhbnd = DataHeap::getInstance().getData(p.getFluctuation()).data();

          fineGridCell.validateNoNansInADERDGSolver(i,fineGridCell,fineGridVerticesEnumerator,"exahype::mappings::SpaceTimePredictor::enterCell[pre]");

          solver->spaceTimePredictor(
              lQi, lFi, lQhi, lFhi,
              lQhbnd,
              lFhbnd,
              luh, fineGridVerticesEnumerator.getCellSize(),
              p.getPredictorTimeStepSize());

          // Perform volume integral
          // TODO(Dominic): This should move into touchVertexLastTime/prepareSendToNeighbour.
          // Here it should be performed directly after we have sent out all the face data
          // to give the network some time to deliver the MPI messages (theoretically).
          solver->volumeIntegral(lduh, lFhi, fineGridVerticesEnumerator.getCellSize());

          fineGridCell.validateNoNansInADERDGSolver(i,fineGridCell,fineGridVerticesEnumerator,"exahype::mappings::SpaceTimePredictor::enterCell[post]");

          break;
        default:
          break;
      }
    endpfor

    peano::datatraversal::autotuning::Oracle::getInstance()
        .parallelSectionHasTerminated(methodTrace);
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

#ifdef Parallel
int exahype::mappings::SpaceTimePredictor::countListingsOfRemoteRankAtFace(
    const int faceIndex,
    exahype::Vertex* const verticesAroundCell,
    const peano::grid::VertexEnumerator& verticesEnumerator) {
  int result = 0;

  const int f = faceIndex % 2;   // "0" indicates a left face, "1" indicates a right face.
  const int d = (faceIndex-f)/2; // The normal direction: 0: x, 1: y, 1: z.

  tarch::la::Vector<DIMENSIONS,int> pos(1); // This is now the center, i.e., (1,1,...,1).
  pos(d) = 2*f;                             // This is a shift from the center by one unit in direction d.

  int faceNeighbourRank = -1; // This variable is introduced to make sure that the adjacent remote rank is unique.
  // TODO(Dominic): Uniqueness is probably guaranteed by the SFC based DD.
  dfor2(v) // Loop over vertices.
    if (verticesAroundCell[ verticesEnumerator(v) ].isAdjacentToRemoteRank()) {
      dfor2(a) // Loop over adjacent ranks. Does also include own rank.
        if (tarch::la::equals(v+a,pos) &&
            verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar]!=
            tarch::parallel::Node::getInstance().getRank()) {
          // Increment
          if (faceNeighbourRank==-1) {
            faceNeighbourRank = verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar];
          }
          if (verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar]==faceNeighbourRank) {
            result++;
          }
        }
      enddforx // a
    }
  enddforx // v

  assertion2(result==0||result==TWO_POWER_D_DIVIDED_BY_TWO,result,faceIndex);

  return result;
}

void exahype::mappings::SpaceTimePredictor::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // TODO(Dominic): Add to docu.
  //
  // Sender:
  // 1. Send out face data if cell type
  // is Cell/Ancestor/Descendant AND if the counter says yes
  // 2. Send out metadata: Send out "empty" metadata if the counters says no
  // or the src cellDescriptionsIndex is not valid.
  //
  // Receiver:
  // 1. Receive metadata
  // 2. If neighbour cell type is Cell/Ancestor/Descendant and the metadata is not "empty"
  // receive face data for cell types Cell/Ancestor/Descendant.
  // Drop messages for EmptyAncestor/EmptyDescendant and if the
  // destCellDescriptionsIndex is not valid.

// TODO(Dominic): This does prevent couting down the vertices.
//#if !defined(PeriodicBC)
//  if (vertex.isBoundary()) return;
//#endif

  // TODO(Dominic): Add to docu why we remove the vertex.isInside() constraint here.
  // We might also consider to remove it from the grid setup mapping functions.
  // fineGridCell.isInside does not imply that all adjacent vertices are
  // inside. If we count down the counter only on
  // fineGridVertices that are inside we might not send out all faces
  // of a cell that is close to the boundary.

  tarch::la::Vector<TWO_POWER_D, int>& adjacentADERDGCellDescriptionsIndices =
      vertex.getCellDescriptionsIndex();

  dfor2(dest)
    dfor2(src)
      if (vertex.hasToSendMetadata(_state,src,dest,toRank)) {
        // we are solely exchanging faces
        const int srcCellDescriptionIndex = adjacentADERDGCellDescriptionsIndices(srcScalar);

        if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(srcCellDescriptionIndex)) {
          assertion1(FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(srcCellDescriptionIndex),srcCellDescriptionIndex);

          decrementCounters(src,dest,srcCellDescriptionIndex);
          if (needToSendFaceData(src,dest,srcCellDescriptionIndex)) {
            sentADERDGFaceData(toRank,x,level,src,dest,srcCellDescriptionIndex,adjacentADERDGCellDescriptionsIndices(destScalar));
            //            sentFiniteVolumesFaceData(toRank,x,level,src,dest,srcCellDescriptionIndex,adjacentADERDGCellDescriptionsIndices(destScalar));
            auto encodedMetadata = exahype::Cell::encodeMetadata(srcCellDescriptionIndex);
            MetadataHeap::getInstance().sendData(
                encodedMetadata, toRank, x, level,
                peano::heap::MessageType::NeighbourCommunication);
          } else {
            logDebug("prepareSendToNeighbour(...)","[empty] sent to rank "<<toRank<<", x:"<<
                x.toString() << ", level=" <<level << ", vertex.adjacentRanks: "
                << vertex.getAdjacentRanks());
            auto encodedMetadata = exahype::Cell::createEncodedMetadataSequenceForInvalidCellDescriptionsIndex();
            MetadataHeap::getInstance().sendData(
                encodedMetadata, toRank, x, level,
                peano::heap::MessageType::NeighbourCommunication);
          }
        } else {
          logDebug("prepareSendToNeighbour(...)","[empty] sent to rank "<<toRank<<", x:"<<
              x.toString() << ", level=" <<level << ", vertex.adjacentRanks: "
              << vertex.getAdjacentRanks());
          auto encodedMetadata = exahype::Cell::createEncodedMetadataSequenceForInvalidCellDescriptionsIndex();
          MetadataHeap::getInstance().sendData(
              encodedMetadata, toRank, x, level,
              peano::heap::MessageType::NeighbourCommunication);
        }
      }
    enddforx
  enddforx
}

bool exahype::mappings::SpaceTimePredictor::needToSendFaceData(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    int cellDescriptionsIndex) {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!

  // ADER-DG
  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)) {
    if (p.getFaceDataExchangeCounter(faceIndex)!=0) {
      return false;
    }
  }

//  // TODO(Dominic): Introduce counters to FiniteVolumesCellDescription.
//  // FV
//  for (auto& p : FiniteVolumesCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)) {
//    if (p.getFaceDataExchangeCounter(faceIndex)!=0) {
//      return false;
//    }
//  }

  return true;
}

void exahype::mappings::SpaceTimePredictor::decrementCounters(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    int cellDescriptionsIndex) {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!

  // ADER-DG
  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)) {
    int newCounterValue = p.getFaceDataExchangeCounter(faceIndex)-1;
    assertion1(newCounterValue<TWO_POWER_D,newCounterValue);
    p.setFaceDataExchangeCounter(faceIndex,newCounterValue);
  }

//  // TODO(Dominic): Introduce counters to FiniteVolumesCellDescription.
//  // FV
//  for (auto& p : FiniteVolumesCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)) {
//    int newCounterValue = p.getFaceDataExchangeCounter(faceIndex)-1;
//    assertion1(newCounterValue<TWO_POWER_D,newCounterValue);
//    p.setFaceDataExchangeCounter(faceIndex,newCounterValue);
//  }
}

void exahype::mappings::SpaceTimePredictor::sentADERDGFaceData(
    int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    int level,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    int srcCellDescriptionIndex,
    int destCellDescriptionIndex) {
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!

  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(srcCellDescriptionIndex)) {
    exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
        exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);
    assertion1(p.getFaceDataExchangeCounter(faceIndex)==0,p.getFaceDataExchangeCounter(faceIndex));

    if (p.getType() == exahype::records::ADERDGCellDescription::Cell     ||
        p.getType() == exahype::records::ADERDGCellDescription::Ancestor ||
        p.getType() == exahype::records::ADERDGCellDescription::Descendant) {

      assertion(DataHeap::getInstance().isValidIndex(p.getExtrapolatedPredictor()));
      assertion(DataHeap::getInstance().isValidIndex(p.getFluctuation()));

      const int numberOfFaceDof = solver->getUnknownsPerFace();
      const double* lQhbnd = DataHeap::getInstance().getData(
          p.getExtrapolatedPredictor()).data() +
          (faceIndex * numberOfFaceDof);
      const double* lFhbnd = DataHeap::getInstance().getData(
          p.getFluctuation()).data() +
          (faceIndex * numberOfFaceDof);

      // TODO(Dominic) This can only happen if we consider
      // vertices that are outside in the MPI communication,
      // i.e. they belong to everyone else. Need to check if we
      // can drop isInside() in some of the mappings.
      // Q: If rank 0 is neighbour to every other rank that is on the boundary,
      // can we use rank 0 to handle the periodic boundary conditions?
      // If rank 0 is forwarding this will however take one extra iteration.
      // Can we resolve the periodic neighbour rank and alter the vertex position+-domain size(x,y,z) before/after we push
      // the metadata and facedata to the heap? The state knows all active ranks? Does he know
      // their respective subdomains?
      // We further need augmentation on the periodic neighbour rank. This can
      // be done by the same periodic neighbour rank resolving and exchanging metadata.
      // PeriodicBC might be easier to implement in a MPI scenario than
      // in a non-MPI scenario.
      if (destCellDescriptionIndex == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex) {
#ifdef PeriodicBC
        assertionMsg(false, "Vasco, we have to implement this");
        DataHeap::getInstance().sendData(
            sentMinMax, toRank, x, level, peano::heap::MessageType::NeighbourCommunication);
        DataHeap::getInstance().sendData(
            lQhbnd, numberOfFaceDof, toRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
        DataHeap::getInstance().sendData(
            lFhbnd, numberOfFaceDof, toRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
#else
        assertionMsg(false, "should never been entered");
#endif
      } else {
        logInfo( //  TODO(Dominic): Change to logDebug.
            "sentADERDGFaceData(...)",
            "send three arrays to rank " <<
            toRank << " for vertex x=" << x << ", level=" << level <<
            ", dest type=" << multiscalelinkedcell::indexToString(destCellDescriptionIndex) <<
            ", src=" << src << ", dest=" << dest <<
            ", counter=" << p.getFaceDataExchangeCounter(faceIndex)
          );

        std::vector<double> sentMinMax(2);
        sentMinMax[0] = p.getSolutionMin(faceIndex);
        sentMinMax[1] = p.getSolutionMax(faceIndex);

        // Send order: minMax,lQhbnd,lFhbnd
        // Receive order: lFhbnd,lQhbnd,minMax
        DataHeap::getInstance().sendData(
            sentMinMax, toRank, x, level, peano::heap::MessageType::NeighbourCommunication);
        DataHeap::getInstance().sendData(
            lQhbnd, numberOfFaceDof, toRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
        DataHeap::getInstance().sendData(
            lFhbnd, numberOfFaceDof, toRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);

        // TODO(Dominic): If anarchic time stepping send the time step over too.
      }
    }

  }
}

bool exahype::mappings::SpaceTimePredictor::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  for (auto& p : exahype::solvers::RegisteredSolvers) {
    p->sendToRank(worker, _mpiTag);
  }

  return true;
}

void exahype::mappings::SpaceTimePredictor::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  for (auto& p : exahype::solvers::RegisteredSolvers) {
    p->receiveFromMasterRank(tarch::parallel::NodePool::getInstance().getMasterRank(),
                       _mpiTag);
  }
}



//
// Below all methods are nop.
//
// ====================================



void exahype::mappings::SpaceTimePredictor::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithMaster(
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

void exahype::mappings::SpaceTimePredictor::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif


exahype::mappings::SpaceTimePredictor::SpaceTimePredictor() {
  // do nothing
}

exahype::mappings::SpaceTimePredictor::~SpaceTimePredictor() {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}


void exahype::mappings::SpaceTimePredictor::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::beginIteration(
    exahype::State& solverState) {
  #ifdef Parallel
  _state = &solverState;
  #endif
}

void exahype::mappings::SpaceTimePredictor::endIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
