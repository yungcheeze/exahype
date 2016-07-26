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
 
#include "exahype/mappings/MarkingForRefinement.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "tarch/la/VectorScalarOperations.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"
#include "kernels/KernelCalls.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::MarkingForRefinement::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::MarkingForRefinement::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::MarkingForRefinement::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::MarkingForRefinement::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::MarkingForRefinement::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::MarkingForRefinement::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::MarkingForRefinement::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::MarkingForRefinement::_log(
    "exahype::mappings::MarkingForRefinement");

exahype::mappings::MarkingForRefinement::MarkingForRefinement() {
  // do nothing
}

exahype::mappings::MarkingForRefinement::~MarkingForRefinement() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::MarkingForRefinement::MarkingForRefinement(
    const MarkingForRefinement& masterThread) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::mergeWithWorkerThread(
    const MarkingForRefinement& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::MarkingForRefinement::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::MarkingForRefinement::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::MarkingForRefinement::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  // do nothing
  return true;
}

void exahype::mappings::MarkingForRefinement::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::mergeWithMaster(
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

void exahype::mappings::MarkingForRefinement::receiveDataFromMaster(
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

void exahype::mappings::MarkingForRefinement::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::MarkingForRefinement::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if ( fineGridCell.isInitialised() ) {
    logInfo( "enterCell(...)", "enter with " << fineGridCell.getNumberOfADERDGCellDescriptions() << " ADER-DG solver(s) and " << fineGridCell.getNumberOfFiniteVolumeCellDescriptions() << " FV solver(s)")

    for (int i=0; i<fineGridCell.getNumberOfADERDGCellDescriptions(); i++) {
      auto* uncastedSolver = exahype::solvers::RegisteredSolvers[ fineGridCell.getADERDGCellDescription(i).getSolverNumber() ];
      assertion( uncastedSolver->getType()==exahype::solvers::Solver::Type::ADER_DG );
      exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(uncastedSolver);

      if (
        tarch::la::allSmallerEquals(fineGridVerticesEnumerator.getCellSize(),solver->getMaximumMeshSize())
        &&
        fineGridCell.getADERDGCellDescription(i).getType()==exahype::records::ADERDGCellDescription::Cell
        &&
        fineGridCell.getADERDGCellDescription(i).getRefinementEvent()==exahype::records::ADERDGCellDescription::None
      ) {
        double*  solution = DataHeap::getInstance().getData(fineGridCell.getADERDGCellDescription(i).getSolution()).data();

        exahype::solvers::Solver::RefinementControl refinementControl = solver->refinementCriterion(
            solution, fineGridVerticesEnumerator.getCellCenter(),
            fineGridVerticesEnumerator.getCellSize(),
            fineGridCell.getADERDGCellDescription(i).getCorrectorTimeStamp(), fineGridCell.getADERDGCellDescription(i).getLevel());

        switch (refinementControl) {
          case exahype::solvers::Solver::RefinementControl::Refine:
            fineGridCell.getADERDGCellDescription(i).setRefinementEvent(exahype::records::ADERDGCellDescription::RefiningRequested);

            dfor2(v)
              if (fineGridVertices[ fineGridVerticesEnumerator(v) ].getRefinementControl()==exahype::Vertex::Records::RefinementControl::Unrefined) {
                fineGridVertices[ fineGridVerticesEnumerator(v) ].refine();
              }
            enddforx
          break;
          case exahype::solvers::Solver::RefinementControl::Erase:
            fineGridCell.getADERDGCellDescription(i).setRefinementEvent(exahype::records::ADERDGCellDescription::ErasingRequested);
          break;
          default:
            break;
        }
      }
    }

    logInfo( "enterCell(...)", "continue  with " << fineGridCell.getNumberOfFiniteVolumeCellDescriptions() << " FV solver(s)")
    for (int i=0; i<fineGridCell.getNumberOfFiniteVolumeCellDescriptions(); i++) {
      auto* uncastedSolver = exahype::solvers::RegisteredSolvers[ fineGridCell.getFiniteVolumesCellDescription(i).getSolverNumber() ];
      assertion( uncastedSolver->getType()==exahype::solvers::Solver::Type::FiniteVolumes );
      exahype::solvers::FiniteVolumesSolver* solver = static_cast<exahype::solvers::FiniteVolumesSolver*>(uncastedSolver);

      if (
        tarch::la::allSmallerEquals(fineGridVerticesEnumerator.getCellSize(),solver->getMaximumMeshSize())
        &&
        fineGridCell.getFiniteVolumesCellDescription(i).getType()==exahype::records::FiniteVolumesCellDescription::Cell
/*
        &&
        fineGridCell.FiniteVolumesSolver(i).getRefinementEvent()==exahype::records::ADERDGCellDescription::None
*/
      ) {
        double*  solution = DataHeap::getInstance().getData(fineGridCell.getFiniteVolumesCellDescription(i).getSolution()).data();

        exahype::solvers::Solver::RefinementControl refinementControl = solver->refinementCriterion(
            solution, fineGridVerticesEnumerator.getCellCenter(),
            fineGridVerticesEnumerator.getCellSize(),
            fineGridCell.getADERDGCellDescription(i).getCorrectorTimeStamp(), fineGridCell.getADERDGCellDescription(i).getLevel());

        switch (refinementControl) {
          case exahype::solvers::Solver::RefinementControl::Refine:
//            fineGridCell.getADERDGCellDescription(i).setRefinementEvent(exahype::records::ADERDGCellDescription::RefiningRequested);

            dfor2(v)
              if (fineGridVertices[ fineGridVerticesEnumerator(v) ].getRefinementControl()==exahype::Vertex::Records::RefinementControl::Unrefined) {
                fineGridVertices[ fineGridVerticesEnumerator(v) ].refine();
              }
            enddforx
          break;
          case exahype::solvers::Solver::RefinementControl::Erase:
            assertionMsg( false, "not implemented yet" );
            //fineGridCell.getADERDGCellDescription(i).setRefinementEvent(exahype::records::ADERDGCellDescription::ErasingRequested);
          break;
          default:
            break;
        }
      }
    }

    logInfo( "enterCell(...)", "left with " << fineGridCell.getNumberOfADERDGCellDescriptions() << " ADER-DG solver(s) and " << fineGridCell.getNumberOfFiniteVolumeCellDescriptions() << " FV solver(s)")
  }

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::MarkingForRefinement::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::beginIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::endIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
