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
 
#include "exahype/mappings/TimeStepSizeComputation.h"

#include "peano/utils/Globals.h"

#include "tarch/multicore/Loop.h"
#include "tarch/parallel/Node.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#include "peano/utils/UserInterface.h"

#include <limits>

peano::CommunicationSpecification
exahype::mappings::TimeStepSizeComputation::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterProcessingOfLocalSubtree,
      true);
}

peano::MappingSpecification
exahype::mappings::TimeStepSizeComputation::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
/**
 * Nop.
 */
peano::MappingSpecification exahype::mappings::TimeStepSizeComputation::
    touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification exahype::mappings::TimeStepSizeComputation::
    touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::TimeStepSizeComputation::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::TimeStepSizeComputation::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::TimeStepSizeComputation::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::TimeStepSizeComputation::_log(
    "exahype::mappings::TimeStepSizeComputation");

void exahype::mappings::TimeStepSizeComputation::prepareLocalTimeStepVariables(){
  int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _minTimeStepSizes.resize(numberOfSolvers);
  _minCellSizes.resize(numberOfSolvers);
  _maxCellSizes.resize(numberOfSolvers);

  for (unsigned int solverNumber=0;
      solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _minTimeStepSizes[solverNumber] = std::numeric_limits<double>::max();
    _minCellSizes    [solverNumber] = std::numeric_limits<double>::max();
    _maxCellSizes    [solverNumber] = -std::numeric_limits<double>::max(); // "-", min
  }
}

void exahype::mappings::TimeStepSizeComputation::prepareTemporaryVariables() {
  if (_tempEigenValues==nullptr) {
    int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
    _tempEigenValues = new double*[numberOfSolvers];

    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      _tempEigenValues[solverNumber]  = new double[solver->getNumberOfVariables()]; // TOOD(Dominic): Check if we need number of parameters too
      ++solverNumber;
    }
  }
}

void exahype::mappings::TimeStepSizeComputation::deleteTemporaryVariables() {
  if(_tempEigenValues!=nullptr) {
    for (unsigned int solverNumber=0;
        solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
      delete[] _tempEigenValues[solverNumber];
      _tempEigenValues[solverNumber] = nullptr;
    }

    delete[] _tempEigenValues;
    _tempEigenValues = nullptr;
  }
}

exahype::mappings::TimeStepSizeComputation::~TimeStepSizeComputation() {
  deleteTemporaryVariables();
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::TimeStepSizeComputation::TimeStepSizeComputation(const TimeStepSizeComputation& masterThread) {
  prepareLocalTimeStepVariables();
  prepareTemporaryVariables();
}

// Merge over threads
void exahype::mappings::TimeStepSizeComputation::mergeWithWorkerThread(
    const TimeStepSizeComputation& workerThread) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    _minTimeStepSizes[i] =
        std::min(_minTimeStepSizes[i], workerThread._minTimeStepSizes[i]);
    _minCellSizes[i] =
        std::min(_minCellSizes[i], workerThread._minCellSizes[i]);
    _maxCellSizes[i] =
        std::max(_maxCellSizes[i], workerThread._maxCellSizes[i]);
  }
}
#endif

void exahype::mappings::TimeStepSizeComputation::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  prepareLocalTimeStepVariables();
  prepareTemporaryVariables();

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::TimeStepSizeComputation::endIteration(
    exahype::State& state) {
  logTraceInWith1Argument("endIteration(State)", state);

  deleteTemporaryVariables();

  state.setStabilityConditionOfOneSolverWasViolated(false);

  int solverNumber=0;
  for (auto solver : exahype::solvers::RegisteredSolvers) {
    assertion1(std::isfinite(_minTimeStepSizes[solverNumber]),_minTimeStepSizes[solverNumber]);
    assertion1(_minTimeStepSizes[solverNumber]>0.0,_minTimeStepSizes[solverNumber]);

    logDebug("endIteration(state)","_minCellSizes[solverNumber]="<<_minCellSizes[solverNumber]<<
             ",_minCellSizes[solverNumber]="<<_maxCellSizes[solverNumber])

    solver->updateNextMinCellSize(_minCellSizes[solverNumber]);
    solver->updateNextMaxCellSize(_maxCellSizes[solverNumber]);

    if (tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
      assertion2(solver->getNextMinCellSize()<std::numeric_limits<double>::max(),solver->getNextMinCellSize(),_minCellSizes[solverNumber]);
      assertion3(solver->getNextMaxCellSize()>0,solver->getNextMaxCellSize(),_maxCellSizes[solverNumber],_minCellSizes[solverNumber]);
    }

    solver->updateMinNextTimeStepSize(_minTimeStepSizes[solverNumber]);
    if (
        (solver->getType()==exahype::solvers::Solver::Type::ADERDG
            || solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) &&
        state.fuseADERDGPhases()
        #ifdef Parallel
        && tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getGlobalMasterRank()
        #endif
    ) {
      auto aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
      const double stableTimeStepSize = aderdgSolver->getMinNextPredictorTimeStepSize();
      const double usedTimeStepSize   = aderdgSolver->getMinPredictorTimeStepSize();
      bool useTimeStepSizeWasInstable = usedTimeStepSize > stableTimeStepSize;

      if (useTimeStepSizeWasInstable) {
        state.setStabilityConditionOfOneSolverWasViolated(true);

        const double timeStepSizeWeight = state.getTimeStepSizeWeightForPredictionRerun();
        aderdgSolver->updateMinNextPredictorTimeStepSize(
            timeStepSizeWeight * stableTimeStepSize);
        aderdgSolver->setMinPredictorTimeStepSize(
            timeStepSizeWeight * stableTimeStepSize);
      } else {
        aderdgSolver->updateMinNextPredictorTimeStepSize(
            0.5 * (stableTimeStepSize + usedTimeStepSize));
      }
    }

    if (state.reinitTimeStepData()) {
      solver->reinitTimeStepData();
    }
    solver->startNewTimeStep();

    ++solverNumber;
  }

  logTraceOutWith1Argument("endIteration(State)", state);
}

// TODO(Dominic): remove the commented out code.
//static double startNewTimeStepFV(
//    const int cellDescriptionsIndex,
//    const int element,
//    exahype::Vertex* const fineGridVertices,
//    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
//  exahype::solvers::FiniteVolumesSolver::CellDescription& p =
//      exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex)[element];
//  exahype::solvers::FiniteVolumesSolver* solver = static_cast<exahype::solvers::FiniteVolumesSolver*>(
//      exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);
//
//  if (p.getType()==exahype::records::FiniteVolumesCellDescription::Cell) {
////         assertion1(p.getRefinementEvent()==exahype::records::FiniteVolumesCellDescription::None,p.toString()); // todo refine
//    const tarch::la::Vector<THREE_POWER_D, int> neighbourCellDescriptionsIndices = multiscalelinkedcell::getIndicesAroundCell(
//                    exahype::VertexOperations::readCellDescriptionsIndex(fineGridVerticesEnumerator, fineGridVertices));
//    assertion1(multiscalelinkedcell::HangingVertexBookkeeper::allAdjacencyInformationIsAvailable(
//        exahype::VertexOperations::readCellDescriptionsIndex(fineGridVerticesEnumerator, fineGridVertices)),
//               fineGridVerticesEnumerator.toString());
//
//    double* finiteVolumesSolutions[THREE_POWER_D];
//    for (int nScalar=0; nScalar<THREE_POWER_D; ++nScalar) {
//      if (exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(neighbourCellDescriptionsIndices[nScalar])) {
//        int element = solver->tryGetElement(cellDescriptionsIndex,p.getSolverNumber());
//        exahype::records::FiniteVolumesCellDescription& pNeighbour =
//           solver->getCellDescription(neighbourCellDescriptionsIndices[nScalar],element);
//        finiteVolumesSolutions[nScalar] = exahype::DataHeap::getInstance().getData(pNeighbour.getSolution()).data();
//      } else {
//        finiteVolumesSolutions[nScalar] = exahype::DataHeap::getInstance().getData(p.getSolution()).data();
//      }
//    }
//
//    double admissibleTimeStepSize = solver->stableTimeStepSize(
//        finiteVolumesSolutions, p.getSize());
//
//    assertion(!std::isnan(admissibleTimeStepSize));
//
//    p.setTimeStamp(p.getTimeStamp()+p.getTimeStepSize());
//    p.setTimeStepSize(admissibleTimeStepSize);
//
//    return admissibleTimeStepSize;
//  }
//
//  return std::numeric_limits<double>::max();
//}

void exahype::mappings::TimeStepSizeComputation::enterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
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
    const int numberOfSolvers = static_cast<int>(exahype::solvers::RegisteredSolvers.size());
    // please use a different UserDefined per mapping/event
    peano::datatraversal::autotuning::MethodTrace methodTrace = peano::datatraversal::autotuning::UserDefined8;
    #ifdef SharedMemoryParallelisation
    int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfSolvers, methodTrace);
    #endif

    pfor(solverNumber, 0, numberOfSolvers, grainSize)
      exahype::solvers::Solver* solver =
          exahype::solvers::RegisteredSolvers[solverNumber];
      int element = exahype::solvers::RegisteredSolvers[solverNumber]->tryGetElement(
          fineGridCell.getCellDescriptionsIndex(),solverNumber);

      if (element!=exahype::solvers::Solver::NotFound) {
        double admissibleTimeStepSize =
            solver->startNewTimeStep(
                fineGridCell.getCellDescriptionsIndex(),element,
                _tempEigenValues[solverNumber]);

        _minTimeStepSizes[solverNumber] = std::min(
            admissibleTimeStepSize, _minTimeStepSizes[solverNumber]);
        _minCellSizes[solverNumber] = std::min(
            fineGridVerticesEnumerator.getCellSize()[0],_minCellSizes[solverNumber]);
        _maxCellSizes[solverNumber] = std::max(
            fineGridVerticesEnumerator.getCellSize()[0],_maxCellSizes[solverNumber]);
      }
    endpfor peano::datatraversal::autotuning::Oracle::getInstance()
        .parallelSectionHasTerminated(methodTrace);
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}


//
// All methods below are nop.
//
//=====================================


#ifdef Parallel

void exahype::mappings::TimeStepSizeComputation::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::mergeWithMaster(
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

bool exahype::mappings::TimeStepSizeComputation::prepareSendToWorker(
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


void exahype::mappings::TimeStepSizeComputation::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::receiveDataFromMaster(
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

void exahype::mappings::TimeStepSizeComputation::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::TimeStepSizeComputation::TimeStepSizeComputation() {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
