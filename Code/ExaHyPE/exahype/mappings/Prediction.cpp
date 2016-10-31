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
 
#include "exahype/mappings/Prediction.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/utils/Globals.h"

#include "tarch/multicore/Loop.h"
#include "tarch/multicore/Lock.h"

#include "peano/utils/Loop.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/solvers/ADERDGSolver.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"

#include "peano/utils/UserInterface.h"

#include <algorithm>

peano::CommunicationSpecification
exahype::mappings::Prediction::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}

peano::MappingSpecification
exahype::mappings::Prediction::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::Prediction::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

// The remainder specs all are nop
peano::MappingSpecification
exahype::mappings::Prediction::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::Prediction::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::Prediction::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::Prediction::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::Prediction::_log(
    "exahype::mappings::Prediction");


void exahype::mappings::Prediction::prepareTemporaryVariables() {
  assertion(_tempSpaceTimeUnknowns    ==nullptr);
  assertion(_tempSpaceTimeFluxUnknowns==nullptr);
  assertion(_tempUnknowns             ==nullptr);
  assertion(_tempFluxUnknowns         ==nullptr);

  int numberOfSolvers        = exahype::solvers::RegisteredSolvers.size();
  _tempSpaceTimeUnknowns     = new double**[numberOfSolvers];
  _tempSpaceTimeFluxUnknowns = new double**[numberOfSolvers];
  _tempUnknowns              = new double* [numberOfSolvers];
  _tempFluxUnknowns          = new double* [numberOfSolvers];

  int solverNumber=0;
  for (auto solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::ADER_DG) {
      auto aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);

      _tempSpaceTimeUnknowns[solverNumber] = new double*[4];
      for (int i=0; i<4; ++i) { // max; see spaceTimePredictorNonlinear
        _tempSpaceTimeUnknowns[solverNumber][i] =
            new double[aderdgSolver->getSpaceTimeUnknownsPerCell()];
      }
      //
      _tempSpaceTimeFluxUnknowns[solverNumber] = new double*[2];
      for (int i=0; i<2; ++i) { // max; see spaceTimePredictorNonlinear
        _tempSpaceTimeFluxUnknowns[solverNumber][i] =
            new double[aderdgSolver->getSpaceTimeFluxUnknownsPerCell()];
      }
      //
      _tempUnknowns    [solverNumber] = new double[aderdgSolver->getUnknownsPerCell()];
      //
      _tempFluxUnknowns[solverNumber] = new double[aderdgSolver->getFluxUnknownsPerCell()];
    } else {
      _tempSpaceTimeUnknowns    [solverNumber] = nullptr;
      _tempSpaceTimeFluxUnknowns[solverNumber] = nullptr;
      _tempUnknowns             [solverNumber] = nullptr;
      _tempFluxUnknowns         [solverNumber] = nullptr;
    }

    ++solverNumber;
  }
}

void exahype::mappings::Prediction::deleteTemporaryVariables() {
  if (_tempSpaceTimeUnknowns!=nullptr) {
    assertion(_tempSpaceTimeFluxUnknowns!=nullptr);
    assertion(_tempUnknowns             !=nullptr);
    assertion(_tempFluxUnknowns         !=nullptr);

    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      if (solver->getType()==exahype::solvers::Solver::Type::ADER_DG) {
        //
        for (int i=0; i<4; ++i) {
          delete[] _tempSpaceTimeUnknowns[solverNumber][i];
        }
        delete[] _tempSpaceTimeUnknowns[solverNumber];
        _tempSpaceTimeUnknowns[solverNumber] = nullptr;
        //
        for (int i=0; i<2; ++i) {
          delete[] _tempSpaceTimeFluxUnknowns[solverNumber][i];
        }
        delete[] _tempSpaceTimeFluxUnknowns[solverNumber];
        _tempSpaceTimeFluxUnknowns[solverNumber] = nullptr;
        //
        delete[] _tempUnknowns[solverNumber];
        _tempUnknowns[solverNumber] = nullptr;
        //
        delete[] _tempFluxUnknowns[solverNumber];
        _tempFluxUnknowns[solverNumber] = nullptr;
      }

      ++solverNumber;
    }

    delete[] _tempSpaceTimeUnknowns;
    delete[] _tempSpaceTimeFluxUnknowns;
    delete[] _tempUnknowns;
    delete[] _tempFluxUnknowns;
    _tempSpaceTimeUnknowns     = nullptr;
    _tempSpaceTimeFluxUnknowns = nullptr;
    _tempUnknowns              = nullptr;
    _tempFluxUnknowns          = nullptr;
  }
}

exahype::mappings::Prediction::Prediction() {
}

exahype::mappings::Prediction::~Prediction() {
  deleteTemporaryVariables();
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::Prediction::Prediction(
    const Prediction& masterThread) :
  _tempSpaceTimeUnknowns(nullptr),
  _tempSpaceTimeFluxUnknowns(nullptr),
  _tempUnknowns(nullptr),
  _tempFluxUnknowns(nullptr) {
  prepareTemporaryVariables();
}

void exahype::mappings::Prediction::mergeWithWorkerThread(
    const Prediction& workerThread) {
}
#endif

void exahype::mappings::Prediction::beginIteration(
    exahype::State& solverState) {
  prepareTemporaryVariables();
}

void exahype::mappings::Prediction::endIteration(
    exahype::State& solverState) {
  deleteTemporaryVariables();
}

void exahype::mappings::Prediction::performPredictionAndVolumeIntegral(exahype::solvers::ADERDGSolver* solver,
                                        exahype::solvers::ADERDGSolver::CellDescription& cellDescription,
                                        exahype::Vertex* const fineGridVertices,
                                        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  if (cellDescription.getType()==exahype::records::ADERDGCellDescription::Cell) {
    assertion1(cellDescription.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,cellDescription.toString());
    solver->validateNoNansInADERDGSolver(cellDescription,fineGridVerticesEnumerator,"exahype::mappings::Prediction::enterCell[pre]");

    solver->performPredictionAndVolumeIntegral(
        cellDescription,
        _tempSpaceTimeUnknowns    [cellDescription.getSolverNumber()],
        _tempSpaceTimeFluxUnknowns[cellDescription.getSolverNumber()],
        _tempUnknowns             [cellDescription.getSolverNumber()],
        _tempFluxUnknowns         [cellDescription.getSolverNumber()]); // todo remove this argument

    solver->validateNoNansInADERDGSolver(cellDescription,fineGridVerticesEnumerator,"exahype::mappings::Prediction::enterCell[post]");
  }
}

void exahype::mappings::Prediction::enterCell(
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
    const int numberOfADERDGCellDescriptions = static_cast<int>(
        exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
            fineGridCell.getCellDescriptionsIndex()).size());
    // please use a different UserDefined per mapping/event
    const peano::datatraversal::autotuning::MethodTrace methodTrace =
        peano::datatraversal::autotuning::UserDefined4;
    const int grainSize =
        peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
            numberOfADERDGCellDescriptions, methodTrace);
    pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
      auto& cellDescription = exahype::solvers::ADERDGSolver::getCellDescription(
              fineGridCell.getCellDescriptionsIndex(),i);

      // TODO(Dominic):
      switch (exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]->getType()) {
      case exahype::solvers::Solver::Type::ADER_DG: {
        exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
            exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);
        solver->synchroniseTimeStepping(fineGridCell.getCellDescriptionsIndex(),i); // Time step synchr. might be done multiple times per traversal; but this is no issue.
        exahype::Cell::resetNeighbourMergeHelperVariables(
            cellDescription,fineGridVertices,fineGridVerticesEnumerator);

        performPredictionAndVolumeIntegral(solver,cellDescription,fineGridVertices,fineGridVerticesEnumerator);
        }
        break;
      case exahype::solvers::Solver::Type::LimitedADERDG: {
        exahype::solvers::LimitedADERDGSolver* solver = static_cast<exahype::solvers::LimitedADERDGSolver*>(
            exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);
        solver->synchroniseTimeStepping(fineGridCell.getCellDescriptionsIndex(),i); // Time step synchr. might be done multiple times per traversal; but this is no issue.
        exahype::Cell::resetNeighbourMergeHelperVariables(
            cellDescription,fineGridVertices,fineGridVerticesEnumerator);

        performPredictionAndVolumeIntegral(solver,cellDescription,fineGridVertices,fineGridVerticesEnumerator);
        }
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
  // TODO(Dominic): Add getters


//
// Below all methods are nop.
//
// ====================================

#ifdef Parallel
void exahype::mappings::Prediction::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

bool exahype::mappings::Prediction::prepareSendToWorker(
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

void exahype::mappings::Prediction::receiveDataFromMaster(
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

void exahype::mappings::Prediction::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::Prediction::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Prediction::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Prediction::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithMaster(
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

void exahype::mappings::Prediction::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::Prediction::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Prediction::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}


void exahype::mappings::Prediction::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Prediction::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::Prediction::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
