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
 
#include "exahype/mappings/TimeStepComputation.h"

#include "peano/utils/Globals.h"

#include "tarch/multicore/Loop.h"
#include "tarch/parallel/Node.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#include <limits>


bool exahype::mappings::TimeStepComputation::SkipReductionInBatchedTimeSteps = false;



peano::CommunicationSpecification
exahype::mappings::TimeStepComputation::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterProcessingOfLocalSubtree,
      true);
}

peano::MappingSpecification
exahype::mappings::TimeStepComputation::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
/**
 * Nop.
 */
peano::MappingSpecification exahype::mappings::TimeStepComputation::
    touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification exahype::mappings::TimeStepComputation::
    touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::TimeStepComputation::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::TimeStepComputation::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::TimeStepComputation::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::TimeStepComputation::_log(
    "exahype::mappings::TimeStepComputation");
int exahype::mappings::TimeStepComputation::_mpiTag =
    tarch::parallel::Node::reserveFreeTag(
        "exahype::mappings::TimeStepComputation");

void exahype::mappings::TimeStepComputation::
    prepareEmptyLocalTimeStepData() {
  _minTimeStepSizes.resize(exahype::solvers::RegisteredSolvers.size());

  for (int i = 0;
       i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    _minTimeStepSizes[i] = std::numeric_limits<double>::max();
  }
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::TimeStepComputation::TimeStepComputation(
    const TimeStepComputation& masterThread) {
  prepareEmptyLocalTimeStepData();
}

// Merge over threads
void exahype::mappings::TimeStepComputation::mergeWithWorkerThread(
    const TimeStepComputation& workerThread) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    _minTimeStepSizes[i] =
        std::min(_minTimeStepSizes[i], workerThread._minTimeStepSizes[i]);
  }
}
#endif

void exahype::mappings::TimeStepComputation::enterCell(
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
    // ADER-DG
    const int numberOfADERDGCellDescriptions = static_cast<int>(ADERDGCellDescriptionHeap::getInstance()
            .getData(fineGridCell.getCellDescriptionsIndex()).size());
    // please use a different UserDefined per mapping/event
    peano::datatraversal::autotuning::MethodTrace methodTrace = peano::datatraversal::autotuning::UserDefined2;
    int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfADERDGCellDescriptions, methodTrace);

    pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
      records::ADERDGCellDescription& p =
          ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[i];
      exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
          exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);

      double* luh = 0;
      double admissibleTimeStepSize;

      if (p.getType()==exahype::records::ADERDGCellDescription::Cell) {
         assertion1(p.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,p.toString());
         luh = DataHeap::getInstance().getData(p.getSolution()).data();

         admissibleTimeStepSize = solver->stableTimeStepSize(
             luh, fineGridVerticesEnumerator.getCellSize());

         assertionEquals(admissibleTimeStepSize,admissibleTimeStepSize);  // assert no nan

         // direct update of the cell description time steps
         p.setCorrectorTimeStamp(p.getPredictorTimeStamp());
         p.setCorrectorTimeStepSize(p.getPredictorTimeStepSize());
         p.setPredictorTimeStamp(p.getPredictorTimeStamp() +
                                 admissibleTimeStepSize);
         p.setPredictorTimeStepSize(admissibleTimeStepSize);

         // indirect update of the solver time step sizes
         _minTimeStepSizes[p.getSolverNumber()] = std::min(
             admissibleTimeStepSize, _minTimeStepSizes[p.getSolverNumber()]);
      }
    endpfor peano::datatraversal::autotuning::Oracle::getInstance()
        .parallelSectionHasTerminated(methodTrace);

    // FINITE VOLUMES
    const int numberOfFiniteVolumesCellDescriptions = static_cast<int>(FiniteVolumesCellDescriptionHeap::getInstance()
            .getData(fineGridCell.getCellDescriptionsIndex()).size());
    // please use a different UserDefined per mapping/event
    methodTrace = peano::datatraversal::autotuning::UserDefined2;
    grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfFiniteVolumesCellDescriptions, methodTrace);

    pfor(i, 0, numberOfFiniteVolumesCellDescriptions, grainSize)
      records::FiniteVolumesCellDescription& p =
          FiniteVolumesCellDescriptionHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[i];
      exahype::solvers::FiniteVolumesSolver* solver = static_cast<exahype::solvers::FiniteVolumesSolver*>(
          exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);

      if (p.getType()==exahype::records::FiniteVolumesCellDescription::Cell) {
//         assertion1(p.getRefinementEvent()==exahype::records::FiniteVolumesCellDescription::None,p.toString()); // todo refine
        assertion1(multiscalelinkedcell::HangingVertexBookkeeper::allAdjacencyInformationIsAvailable(
            VertexOperations::readCellDescriptionsIndex(fineGridVerticesEnumerator, fineGridVertices)),fineGridVerticesEnumerator.toString());

        const tarch::la::Vector<THREE_POWER_D, int> neighbourCellDescriptionsIndices = multiscalelinkedcell::getIndicesAroundCell(
            VertexOperations::readCellDescriptionsIndex(fineGridVerticesEnumerator, fineGridVertices));

        double* finiteVolumesSolutions[THREE_POWER_D];
        for (int nScalar=0; nScalar<THREE_POWER_D; ++nScalar) {
          if (FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(neighbourCellDescriptionsIndices[nScalar])) {
            exahype::records::FiniteVolumesCellDescription& pNeighbour =
                FiniteVolumesCellDescriptionHeap::getInstance().getData(neighbourCellDescriptionsIndices[nScalar])[p.getSolverNumber()]; // todo assumes same number of patches per cell
            finiteVolumesSolutions[nScalar] = DataHeap::getInstance().getData(pNeighbour.getSolution()).data();
          } else {
            finiteVolumesSolutions[nScalar] = DataHeap::getInstance().getData(p.getSolution()).data();
          }
        }

        double admissibleTimeStepSize = solver->stableTimeStepSize(
           finiteVolumesSolutions, fineGridVerticesEnumerator.getCellSize());

        assertion(!std::isnan(admissibleTimeStepSize));

        p.setTimeStamp(p.getTimeStamp()+p.getTimeStepSize());
        p.setTimeStepSize(admissibleTimeStepSize);

        _minTimeStepSizes[p.getSolverNumber()] = std::min(
           admissibleTimeStepSize, _minTimeStepSizes[p.getSolverNumber()]); // todo MPI
      }
    endpfor peano::datatraversal::autotuning::Oracle::getInstance()
        .parallelSectionHasTerminated(methodTrace);
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::TimeStepComputation::beginIteration(
    exahype::State& solverState) {
  prepareEmptyLocalTimeStepData();
}

void exahype::mappings::TimeStepComputation::endIteration(
    exahype::State& solverState) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[i];

    logDebug("mergeLocalTimeStepDataIntoSolvers()",
             "solver " << i << " is updated with time step size "
                       << _minTimeStepSizes[i]);
    solver->updateNextTimeStepSize(_minTimeStepSizes[i]);
    solver->startNewTimeStep();
  }
}

#ifdef Parallel
// Worker-master comm. Send to master.
void exahype::mappings::TimeStepComputation::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // TODO(Dominic): Merge with other Master-Worker exchanges: FaceUnknownsProjection
  DataHeap::getInstance().startToSendSynchronousData(); // See method documentation.

  for (auto dt : _minTimeStepSizes) {
    assertion1(dt>0,dt);
  }

  DataHeap::getInstance().sendData(
      _minTimeStepSizes.data(),_minTimeStepSizes.size(),
      tarch::parallel::NodePool::getInstance().getMasterRank(),
      verticesEnumerator.getCellCenter(),
      verticesEnumerator.getLevel(),
      peano::heap::MessageType::MasterWorkerCommunication);
}

// Worker-master comm. Receive from worker.
void exahype::mappings::TimeStepComputation::mergeWithMaster(
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
  for (auto dt : _minTimeStepSizes) {
    assertion1(dt>0.0,dt);
  } // Dead code elimination will get rid of this loop if Debug/Asserts flag is not set.

  std::vector<double> receivedMinTimeStepSizes(_minTimeStepSizes.size());
  DataHeap::getInstance().receiveData(
      receivedMinTimeStepSizes.data(),receivedMinTimeStepSizes.size(),
      worker,
      fineGridVerticesEnumerator.getCellCenter(),
      fineGridVerticesEnumerator.getLevel(),
      peano::heap::MessageType::MasterWorkerCommunication);

  DataHeap::getInstance().finishedToSendSynchronousData(); // See method documentation.

  for (int i = 0; i < static_cast<int>(_minTimeStepSizes.size()); i++) {
    _minTimeStepSizes[i] =
        std::min(_minTimeStepSizes[i], receivedMinTimeStepSizes[i]);
  }

  for (auto dt : receivedMinTimeStepSizes) {
    assertion1(dt>0.0,dt);
  } // Dead code elimination will get rid of this.
  for (auto dt : _minTimeStepSizes) {
    assertion1(dt>0.0,dt);
  } // Dead code elimination will get rid of this loop if Debug/Asserts flag is not set.
}


bool exahype::mappings::TimeStepComputation::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  if (
    !peano::parallel::loadbalancing::Oracle::getInstance().isLoadBalancingActivated()
    &&
    SkipReductionInBatchedTimeSteps
  ) {
    return false;
  }
  else return true;
}


void exahype::mappings::TimeStepComputation::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::receiveDataFromMaster(
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

void exahype::mappings::TimeStepComputation::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::TimeStepComputation::TimeStepComputation() {
  // do nothing
}

exahype::mappings::TimeStepComputation::~TimeStepComputation() {
  // do nothing
}

void exahype::mappings::TimeStepComputation::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::TimeStepComputation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
