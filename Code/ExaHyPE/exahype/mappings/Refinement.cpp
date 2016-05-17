#include "exahype/mappings/Refinement.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "exahype/solvers/Solver.h"
#include "kernels/KernelCalls.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::Refinement::communicationSpecification() {
    return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::Refinement::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Refinement::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Refinement::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Refinement::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Refinement::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Refinement::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,
                                     peano::MappingSpecification::Serial);
}

tarch::logging::Log exahype::mappings::Refinement::_log(
    "exahype::mappings::Refinement");

exahype::mappings::Refinement::Refinement() {
  // do nothing
}

exahype::mappings::Refinement::~Refinement() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::Refinement::Refinement(const Refinement& masterThread) {
  // do nothing
}

void exahype::mappings::Refinement::mergeWithWorkerThread(
    const Refinement& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::Refinement::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Refinement::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Refinement::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Refinement::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Refinement::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Refinement::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("createCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  fineGridCell.getCellData().setADERDGCellDescriptionsIndex(
      exahype::Cell::InvalidCellDescriptionsIndex);

  logTraceOutWith1Argument("createCell(...)", fineGridCell);
}

void exahype::mappings::Refinement::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::Refinement::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::Refinement::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Refinement::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Refinement::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Refinement::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Refinement::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::Refinement::prepareSendToWorker(
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

void exahype::mappings::Refinement::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Refinement::mergeWithMaster(
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

void exahype::mappings::Refinement::receiveDataFromMaster(
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

void exahype::mappings::Refinement::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Refinement::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::Refinement::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Refinement::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("touchVertexLastTime(...)", fineGridVertex,
                           fineGridX, fineGridH,
                           coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  tarch::la::Vector<TWO_POWER_D, int>& cellDescriptions =
      fineGridVertex.getADERDGCellDescriptionsIndex();
  bool refineFineGridVertex = false;
  // clang-format off
  // Loop over the 2^d adjacent cells and check if refinement is necessary.
  // If so, refine vertex
  dfor2(c)
  if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      cellDescriptions[cScalar])) {
    for (std::vector<exahype::records::ADERDGCellDescription>::iterator pFine =
        ADERDGCellDescriptionHeap::getInstance()
        .getData(cellDescriptions[cScalar])
        .begin();
        pFine !=
            ADERDGCellDescriptionHeap::getInstance()
        .getData(cellDescriptions[cScalar])
        .end();
        ++pFine) {
      assertion3(static_cast<unsigned int>(pFine->getSolverNumber()) <
          solvers::RegisteredSolvers.size(),
          pFine->getSolverNumber(), solvers::RegisteredSolvers.size(),
          toString(coarseGridVerticesEnumerator.getCellFlags()));

      switch (pFine->getType()) {
      case exahype::records::ADERDGCellDescription::Cell:
        switch (pFine->getRefinementEvent()) {
        case exahype::records::ADERDGCellDescription::RefiningRequested:
          refineFineGridVertex = true;
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
  enddforx  // clang-format on
      // refine vertex
      if (refineFineGridVertex) {
    if (!fineGridVertex.isHangingNode()  // todo discuss with Tobias
        && !fineGridVertex.isRefinedOrRefining()) {
      fineGridVertex.refine();
    }
  }
  logTraceOutWith1Argument("touchVertexLastTime(...)", fineGridVertex);
}

void exahype::mappings::Refinement::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
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
        case exahype::records::ADERDGCellDescription::Refining:
            if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(fineGridCell.getADERDGCellDescriptionsIndex())) {
              for (std::vector<exahype::records::ADERDGCellDescription>::iterator
                       pFine = ADERDGCellDescriptionHeap::getInstance()
                                   .getData(fineGridCell
                                                .getADERDGCellDescriptionsIndex())
                                   .begin();
                   pFine !=
                   ADERDGCellDescriptionHeap::getInstance()
                       .getData(fineGridCell
                                    .getADERDGCellDescriptionsIndex())
                       .end();
                   ++pFine) {
                assertion(pCoarse->getType() ==
                          exahype::records::ADERDGCellDescription::Cell);

                if (pFine->getSolverNumber() == pCoarse->getSolverNumber()) {
                  assertion(
                      exahype::records::ADERDGCellDescription::Descendant ||
                      exahype::records::ADERDGCellDescription::EmptyDescendant);

                  pFine->setType(exahype::records::ADERDGCellDescription::Cell);
                  fineGridCell.ensureNecessaryMemoryIsAllocated(
                      pFine->getSolverNumber());

                  solverNumberFound = true;
                }
              }
            }

            if (!solverNumberFound) {
              fineGridCell.addNewCellDescription(
                  pCoarse->getSolverNumber(),
                  exahype::records::ADERDGCellDescription::Cell,
                  exahype::records::ADERDGCellDescription::Prolongating,
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
}

void exahype::mappings::Refinement::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Refinement::beginIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::Refinement::endIteration(exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::Refinement::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::Refinement::descend(
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
             pCoarse =
                 ADERDGCellDescriptionHeap::getInstance()
                     .getData(coarseGridCell.getADERDGCellDescriptionsIndex())
                     .begin();
         pCoarse !=
         ADERDGCellDescriptionHeap::getInstance()
             .getData(coarseGridCell.getADERDGCellDescriptionsIndex())
             .end();
         ++pCoarse) {
      switch (pCoarse->getRefinementEvent()) {
        case exahype::records::ADERDGCellDescription::RefiningRequested:
          pCoarse->setRefinementEvent(exahype::records::ADERDGCellDescription::Refining);
          break;
        default:
          break;
      }
    }
  }
  logTraceOut("descend(...)");
}
