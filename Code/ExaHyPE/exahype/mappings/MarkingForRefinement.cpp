#include "exahype/mappings/MarkingForRefinement.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "kernels/KernelCalls.h"
#include "exahype/solvers/Solver.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::MarkingForRefinement::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterLastTouchVertexLastTime,
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
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
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
      peano::MappingSpecification::WholeTree,
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
exahype::mappings::MarkingForRefinement::MarkingForRefinement(const MarkingForRefinement& masterThread) {
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

void exahype::mappings::MarkingForRefinement::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::mergeWithRemoteDataDueToForkOrJoin(
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
  if (ADERDGCellDescriptionHeap::getInstance().
      isValidIndex(fineGridCell.getADERDGCellDescriptionsIndex())) {

    for (std::vector<exahype::records::ADERDGCellDescription>::
        iterator pFine = ADERDGCellDescriptionHeap::getInstance().getData(
            fineGridCell.getADERDGCellDescriptionsIndex()).begin();
        pFine != ADERDGCellDescriptionHeap::getInstance().getData(
            fineGridCell.getADERDGCellDescriptionsIndex()).end();
        ++pFine) {
      switch (pFine->getRefinementEvent()) {
      case exahype::records::ADERDGCellDescription::None:
        switch (pFine->getType()) {
        case exahype::records::ADERDGCellDescription::Cell:
          double* solution = DataHeap::getInstance().getData(
              pFine->getSolution()).data();

          exahype::solvers::Solver* solver =
              exahype::solvers::RegisteredSolvers[pFine->getSolverNumber()];
          exahype::solvers::Solver::RefinementControl refinementControl =
              solver->refinementCriterion(solution,
                  fineGridVerticesEnumerator.getCellCenter(),
                  fineGridVerticesEnumerator.getCellSize(),
                  pFine->getCorrectorTimeStamp(), // todo careful with the time stamps
                  pFine->getLevel());

          switch (refinementControl) {
          case exahype::solvers::Solver::Refine:
            pFine->setRefinementEvent(exahype::records::ADERDGCellDescription::Refining);
            break;
          case exahype::solvers::Solver::Coarsen:
            pFine->setRefinementEvent(exahype::records::ADERDGCellDescription::ErasingRequested);
            break;
          }
          break;
        }
        break;
      }
    }
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

void exahype::mappings::MarkingForRefinement::endIteration(exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  logTraceInWith2Arguments( "descend(...)", coarseGridCell.toString(), coarseGridVerticesEnumerator.toString() );

  if (ADERDGCellDescriptionHeap::getInstance().
      isValidIndex(coarseGridCell.getADERDGCellDescriptionsIndex())) {
    for (std::vector<exahype::records::ADERDGCellDescription>::
        iterator pCoarse = ADERDGCellDescriptionHeap::getInstance().getData(
            coarseGridCell.getADERDGCellDescriptionsIndex()).begin();
        pCoarse != ADERDGCellDescriptionHeap::getInstance().getData(
            coarseGridCell.getADERDGCellDescriptionsIndex()).end();
        ++pCoarse) {
      switch (pCoarse->getType()) {
      case exahype::records::ADERDGCellDescription::EmptyAncestor:
      case exahype::records::ADERDGCellDescription::Ancestor:
        switch (pCoarse->getRefinementEvent()) {
        case exahype::records::ADERDGCellDescription::None:
          bool eraseChildren = true;

          dfor3(k)
          for (std::vector<exahype::records::ADERDGCellDescription>::
              iterator pFine = ADERDGCellDescriptionHeap::getInstance().
              getData(fineGridCells[kScalar].
                  getADERDGCellDescriptionsIndex()).begin();
              pFine != ADERDGCellDescriptionHeap::getInstance().
                  getData(fineGridCells[kScalar].
                      getADERDGCellDescriptionsIndex()).end();
              ++pFine) {
            if (pCoarse->getSolverNumber()==pFine->getSolverNumber()) {
              eraseChildren = eraseChildren && pFine->getRefinementEvent()==
                  exahype::records::ADERDGCellDescription::ErasingRequested;
            }
          }
          enddforx

          if (eraseChildren) {
            pCoarse->setRefinementEvent(exahype::records::ADERDGCellDescription::ErasingChildren);
            // @todo: 16/04/27:Dominic Etienne Charier:
            // pCoarse->setType(exahype::records::ADERDGCellDescription::Cell);
            // Clean the solution field of the newly initialised cell.
            coarseGridCell.initialiseCellDescription(pCoarse->getSolverNumber());

            dfor3(k)
            for (std::vector<exahype::records::ADERDGCellDescription>::
                iterator pFine = ADERDGCellDescriptionHeap::getInstance().
                getData(fineGridCells[kScalar].
                    getADERDGCellDescriptionsIndex()).begin();
                pFine != ADERDGCellDescriptionHeap::getInstance().
                    getData(fineGridCells[kScalar].
                        getADERDGCellDescriptionsIndex()).end();
                ++pFine) {
              if (pCoarse->getSolverNumber()==pFine->getSolverNumber()) {
                assertion1(pFine->getRefinementEvent()==exahype::records::ADERDGCellDescription::ErasingRequested,
                           toString());
                pFine->setRefinementEvent(exahype::records::ADERDGCellDescription::Restricting);
                break;
                }
              }
            }
            enddforx
          }
          break;
        }
        break;
      }
    }

    logTraceOut( "descend(...)" );
  }

  void exahype::mappings::MarkingForRefinement::ascend(
      exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell) {
    // do nothing
  }
