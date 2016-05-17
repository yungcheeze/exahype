#include "exahype/mappings/Augmentation.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "exahype/solvers/Solver.h"
#include "kernels/KernelCalls.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::Augmentation::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::Augmentation::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Augmentation::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Augmentation::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Augmentation::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Augmentation::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Augmentation::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial);
}

tarch::logging::Log exahype::mappings::Augmentation::_log(
    "exahype::mappings::Augmentation");

exahype::mappings::Augmentation::Augmentation() {
  // do nothing
}

exahype::mappings::Augmentation::~Augmentation() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::Augmentation::Augmentation(
    const Augmentation& masterThread) {
  // do nothing
}

void exahype::mappings::Augmentation::mergeWithWorkerThread(
    const Augmentation& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::Augmentation::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Augmentation::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Augmentation::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Augmentation::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Augmentation::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Augmentation::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Augmentation::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::Augmentation::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::Augmentation::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Augmentation::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Augmentation::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Augmentation::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Augmentation::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::Augmentation::prepareSendToWorker(
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

void exahype::mappings::Augmentation::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Augmentation::mergeWithMaster(
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

void exahype::mappings::Augmentation::receiveDataFromMaster(
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

void exahype::mappings::Augmentation::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Augmentation::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::Augmentation::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Augmentation::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Augmentation::enterCell(
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
      coarseGridCell.getADERDGCellDescriptionsIndex())) {
    for (std::vector<exahype::records::ADERDGCellDescription>::iterator
        pCoarse =
            ADERDGCellDescriptionHeap::getInstance()
        .getData(coarseGridCell.getADERDGCellDescriptionsIndex())
        .begin();
        pCoarse !=
            ADERDGCellDescriptionHeap::getInstance()
        .getData(coarseGridCell.getADERDGCellDescriptionsIndex())
        .end();
        ++pCoarse) {
      bool solverNumberFound = false;

      switch (pCoarse->getRefinementEvent()) {
      case exahype::records::ADERDGCellDescription::Augmenting:
        if (fineGridCell.getADERDGCellDescriptionsIndex() !=
            exahype::Cell::InvalidCellDescriptionsIndex) {
          for (std::vector<exahype::records::ADERDGCellDescription>::iterator
              pFine = ADERDGCellDescriptionHeap::getInstance()
              .getData(fineGridCell.getADERDGCellDescriptionsIndex()).begin();
              pFine != ADERDGCellDescriptionHeap::getInstance()
              .getData(fineGridCell.getADERDGCellDescriptionsIndex()).end();
              ++pFine) {
            assertion(pCoarse->getType() ==
                exahype::records::ADERDGCellDescription::Cell ||
                exahype::records::ADERDGCellDescription::Descendant ||
                exahype::records::ADERDGCellDescription::EmptyDescendant
            );

            if (pFine->getSolverNumber() == pCoarse->getSolverNumber()) {
              assertion(
                  exahype::records::ADERDGCellDescription::Descendant ||
                  exahype::records::ADERDGCellDescription::EmptyDescendant);
              solverNumberFound = true;
            }
          }
        }

        // We set the default type of the new cell description
        // to EmptyDescendant.
        //
        // !!! Rationale
        //
        // It's more likely that a descendant does not need to hold data.
        // The multiscalelinkedcell is furthermore a little slower at
        // subgrid boundaries. Therefore, we would need one extra
        // iteration.
        if (!solverNumberFound) {
          fineGridCell.addNewCellDescription(
              pCoarse->getSolverNumber(),
              exahype::records::ADERDGCellDescription::EmptyDescendant,
              exahype::records::ADERDGCellDescription::None,
              fineGridVerticesEnumerator.getLevel(),
              coarseGridCell.getADERDGCellDescriptionsIndex(),
              fineGridVerticesEnumerator.getCellSize(),
              fineGridVerticesEnumerator.getCellCenter());
        }
        break;
      default:
        break;
      }
    }
  }

  logTraceOut("enterCell(...)");
}

void exahype::mappings::Augmentation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Augmentation::beginIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::Augmentation::endIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::Augmentation::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  logTraceInWith2Arguments("descend(...)", coarseGridCell.toString(),
      coarseGridVerticesEnumerator.toString());

  if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      coarseGridCell.getADERDGCellDescriptionsIndex())) {
    for (std::vector<exahype::records::ADERDGCellDescription>::iterator
        pCoarse = ADERDGCellDescriptionHeap::getInstance()
        .getData(coarseGridCell.getADERDGCellDescriptionsIndex())
        .begin();
        pCoarse != ADERDGCellDescriptionHeap::getInstance()
        .getData(coarseGridCell.getADERDGCellDescriptionsIndex())
        .end();
        ++pCoarse) {
      switch (pCoarse->getRefinementEvent()) {
      case exahype::records::ADERDGCellDescription::AugmentingRequested:
        pCoarse->setRefinementEvent(
            exahype::records::ADERDGCellDescription::Augmenting);
        break;
      default:
        break;
      }
    }
  }

  logTraceOut("descend(...)");
}

void exahype::mappings::Augmentation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  logTraceInWith2Arguments("ascend(...)", coarseGridCell.toString(),
      coarseGridVerticesEnumerator.toString());

  if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      coarseGridCell.getADERDGCellDescriptionsIndex())) {
    for (std::vector<exahype::records::ADERDGCellDescription>::iterator
        pCoarse =
            ADERDGCellDescriptionHeap::getInstance()
        .getData(coarseGridCell.getADERDGCellDescriptionsIndex())
        .begin();
        pCoarse !=
            ADERDGCellDescriptionHeap::getInstance()
        .getData(coarseGridCell.getADERDGCellDescriptionsIndex())
        .end();
        ++pCoarse) {
      bool augmentingDone = true;

      switch (pCoarse->getType()) {
      case exahype::records::ADERDGCellDescription::Cell:
      case exahype::records::ADERDGCellDescription::Descendant:
      case exahype::records::ADERDGCellDescription::EmptyDescendant:
        switch (pCoarse->getRefinementEvent()) {
        case exahype::records::ADERDGCellDescription::Augmenting:
          augmentingDone = true;
          // clang-format off
          dfor3(k)
          if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
              fineGridCells[kScalar].getADERDGCellDescriptionsIndex())) {
            for (std::vector<exahype::records::ADERDGCellDescription>::
                iterator pFine =
                    ADERDGCellDescriptionHeap::getInstance()
                .getData(
                    fineGridCells[kScalar]
                                  .getADERDGCellDescriptionsIndex())
                                  .begin();
                pFine !=
                    ADERDGCellDescriptionHeap::getInstance()
                .getData(fineGridCells[kScalar]
                                       .getADERDGCellDescriptionsIndex())
                                       .end();
                ++pFine) {
              if (pCoarse->getSolverNumber() == pFine->getSolverNumber()) {
                augmentingDone =
                    augmentingDone &&
                    (pFine->getRefinementEvent() ==
                        exahype::records::ADERDGCellDescription::None
                                            ||
                    pFine->getRefinementEvent() ==
                        exahype::records::ADERDGCellDescription::DeaugmentingRequested);
              }
            }
          } else {
            augmentingDone = false;
          }
          enddforx

          if (augmentingDone) {
            pCoarse->setRefinementEvent(
                exahype::records::ADERDGCellDescription::None);
          }
          break;
        default:
          break;
        }
        break;
        default:
          break;
      }
    }
  }

  logTraceOut("ascend(...)");
}
