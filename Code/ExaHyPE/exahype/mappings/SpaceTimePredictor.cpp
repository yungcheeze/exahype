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

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::SpaceTimePredictor::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
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

exahype::mappings::SpaceTimePredictor::SpaceTimePredictor() {
  // do nothing
}

exahype::mappings::SpaceTimePredictor::~SpaceTimePredictor() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::SpaceTimePredictor::SpaceTimePredictor(
    const SpaceTimePredictor& masterThread) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithWorkerThread(
    const SpaceTimePredictor& workerThread) {
  // do nothing
}
#endif

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

#ifdef Parallel
void exahype::mappings::SpaceTimePredictor::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
#if !defined(PeriodicBC)
  if (vertex.isBoundary()) return;
#endif

  tarch::la::Vector<TWO_POWER_D, int>& adjacentADERDGCellDescriptionsIndices =
      vertex.getCellDescriptionsIndex();

  // @todo Hier stimmen die Abfolgen nicht!

  dfor2(dest)
      dfor2(src) if (vertex.getAdjacentRanks()(destScalar) == toRank &&
                     vertex.getAdjacentRanks()(srcScalar) ==
                         tarch::parallel::Node::getInstance().getRank() &&
                     tarch::la::countEqualEntries(dest, src) ==
                         1  // we are solely exchanging faces
                     ) {
    const int srcCellDescriptionIndex =
        adjacentADERDGCellDescriptionsIndices(srcScalar);
    if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
            srcCellDescriptionIndex)) {
      std::vector<records::ADERDGCellDescription>& cellDescriptions =
          ADERDGCellDescriptionHeap::getInstance().getData(
              srcCellDescriptionIndex);

      for (int currentSolver = 0;
           currentSolver < static_cast<int>(cellDescriptions.size());
           currentSolver++) {
        if (cellDescriptions[currentSolver].getType() ==
            exahype::records::ADERDGCellDescription::Cell) {
          exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
              exahype::solvers::RegisteredSolvers
              [cellDescriptions[currentSolver].getSolverNumber()]);

          const int numberOfFaceDof = solver->getUnknownsPerFace();
          const int normalOfExchangedFace =
              tarch::la::equalsReturnIndex(src, dest);
          assertion(normalOfExchangedFace >= 0 &&
                    normalOfExchangedFace < DIMENSIONS);
          const int offsetInFaceArray =
              2 * normalOfExchangedFace +
              (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1
                                                                        : 0);

          assertion(DataHeap::getInstance().isValidIndex(
              cellDescriptions[currentSolver].getExtrapolatedPredictor()));
          assertion(DataHeap::getInstance().isValidIndex(
              cellDescriptions[currentSolver].getFluctuation()));

          const double* lQhbnd = DataHeap::getInstance()
                                     .getData(cellDescriptions[currentSolver]
                                                  .getExtrapolatedPredictor())
                                     .data() +
                                 (offsetInFaceArray * numberOfFaceDof);
          const double* lFhbnd =
              DataHeap::getInstance()
                  .getData(cellDescriptions[currentSolver].getFluctuation())
                  .data() +
              (offsetInFaceArray * numberOfFaceDof);

          if (adjacentADERDGCellDescriptionsIndices(destScalar) ==
              multiscalelinkedcell::HangingVertexBookkeeper::
                  DomainBoundaryAdjacencyIndex) {
#ifdef PeriodicBC
            assertionMsg(false, "Vasco, we have to implement this");
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
            logDebug(
                "prepareSendToNeighbour(...)",
                "send two arrays to rank "
                    << toRank << " for vertex " << vertex.toString()
                    << ", dest type="
                    << multiscalelinkedcell::indexToString(
                           adjacentADERDGCellDescriptionsIndices(destScalar))
                    << ", src=" << src << ", dest=" << dest);
            DataHeap::getInstance().sendData(
                lQhbnd, numberOfFaceDof, toRank, x, level,
                peano::heap::MessageType::NeighbourCommunication);
            DataHeap::getInstance().sendData(
                lFhbnd, numberOfFaceDof, toRank, x, level,
                peano::heap::MessageType::NeighbourCommunication);
          }
        } else {
          assertionMsg(false, "Dominic, please implement");
        }
      }
    }
  }
  enddforx enddforx
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
    p->receiveFromRank(tarch::parallel::NodePool::getInstance().getMasterRank(),
                       _mpiTag);
  }
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

  if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
          fineGridCell.getCellDescriptionsIndex())) {
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

      // space-time DoF (basisSize**(DIMENSIONS+1))
      double* lQi = 0;
      double* lFi = 0;

      // volume DoF (basisSize**(DIMENSIONS))
      double* luh = 0;
      double* lQhi = 0;
      double* lFhi = 0;

      // face DoF (basisSize**(DIMENSIONS-1))
      double* lQhbnd = 0;
      double* lFhbnd = 0;

      switch (p.getType()) {
        case exahype::records::ADERDGCellDescription::Cell:
          assertion1(p.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,p.toString());
          lQi = DataHeap::getInstance().getData(p.getSpaceTimePredictor()).data();
          lFi = DataHeap::getInstance().getData(p.getSpaceTimeVolumeFlux()).data();

          luh = DataHeap::getInstance().getData(p.getSolution()).data();
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
  // do nothing
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
