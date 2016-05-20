#include "exahype/mappings/MarkingForAugmentation.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"
#include "exahype/VertexOperations.h"

#include "kernels/KernelCalls.h"
#include "exahype/solvers/Solver.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::MarkingForAugmentation::communicationSpecification() {
    return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::MarkingForAugmentation::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::MarkingForAugmentation::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::MarkingForAugmentation::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::MarkingForAugmentation::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::MarkingForAugmentation::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::MarkingForAugmentation::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::MarkingForAugmentation::_log(
    "exahype::mappings::MarkingForAugmentation");

exahype::mappings::MarkingForAugmentation::MarkingForAugmentation() {
  // do nothing
}

exahype::mappings::MarkingForAugmentation::~MarkingForAugmentation() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::MarkingForAugmentation::MarkingForAugmentation(const MarkingForAugmentation& masterThread) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::mergeWithWorkerThread(
    const MarkingForAugmentation& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::MarkingForAugmentation::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::MarkingForAugmentation::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::MarkingForAugmentation::prepareSendToWorker(
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

void exahype::mappings::MarkingForAugmentation::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::mergeWithMaster(
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

void exahype::mappings::MarkingForAugmentation::receiveDataFromMaster(
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

void exahype::mappings::MarkingForAugmentation::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::MarkingForAugmentation::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::enterCell(
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
      isValidIndex(fineGridCell.getADERDGCellDescriptionsIndex())
      &&
      multiscalelinkedcell::HangingVertexBookkeeper::
                allAdjacencyInformationIsAvailable(VertexOperations::readADERDGCellDescriptionsIndex(
                fineGridVerticesEnumerator,fineGridVertices))
  ){
    const tarch::la::Vector<THREE_POWER_D,int> neighbourCellDescriptionIndices =
        multiscalelinkedcell::getIndicesAroundCell(
            VertexOperations::readADERDGCellDescriptionsIndex(
                fineGridVerticesEnumerator,fineGridVertices));

    for (std::vector<exahype::records::ADERDGCellDescription>::
        iterator pFine = ADERDGCellDescriptionHeap::getInstance().getData(
            fineGridCell.getADERDGCellDescriptionsIndex()).begin();
        pFine != ADERDGCellDescriptionHeap::getInstance().getData(
            fineGridCell.getADERDGCellDescriptionsIndex()).end();
        ++pFine) {
      const exahype::solvers::Solver::AugmentationControl augmentationControl =
          augmentationCriterion(
              pFine->getSolverNumber(),
              pFine->getType(),
              pFine->getLevel(),
              neighbourCellDescriptionIndices);

      // 1. Check if ancestors and descendants need to hold
      // data or not based on virtual refinement criterion.
      switch (pFine->getType()) {
      case exahype::records::ADERDGCellDescription::Ancestor:
      case exahype::records::ADERDGCellDescription::EmptyAncestor:
        switch (augmentationControl) {
        case exahype::solvers::Solver::NextToCell:
        case exahype::solvers::Solver::NextToCellAndAncestor:
          pFine->setType(exahype::records::ADERDGCellDescription::Ancestor);
          break;
        default:
          pFine->setType(exahype::records::ADERDGCellDescription::EmptyAncestor);
          break;
        }
        break;
      case exahype::records::ADERDGCellDescription::Descendant:
      case exahype::records::ADERDGCellDescription::EmptyDescendant:
        switch (augmentationControl) {
        case exahype::solvers::Solver::NextToCell:
        case exahype::solvers::Solver::NextToCellAndAncestor:
          pFine->setType(exahype::records::ADERDGCellDescription::Descendant);
          break;
        default:
          pFine->setType(exahype::records::ADERDGCellDescription::EmptyDescendant);
          break;
        }
        break;
      default:
        break;
      }

      // 2. Further augment or deaugment cells and descendants if no other event or an
      // augmentation event has been
      // triggered.
      switch (pFine->getRefinementEvent()) {
      case exahype::records::ADERDGCellDescription::None:
      case exahype::records::ADERDGCellDescription::AugmentingRequested:
      case exahype::records::ADERDGCellDescription::DeaugmentingRequested:
        switch (pFine->getType()) {
        case exahype::records::ADERDGCellDescription::Cell:
        case exahype::records::ADERDGCellDescription::Descendant:
        case exahype::records::ADERDGCellDescription::EmptyDescendant:
          switch (augmentationControl) {
          case exahype::solvers::Solver::NextToAncestor:
          case exahype::solvers::Solver::NextToCellAndAncestor:
            pFine->setRefinementEvent(exahype::records::ADERDGCellDescription::AugmentingRequested);
            break;
          case exahype::solvers::Solver::Unncessary:
            pFine->setRefinementEvent(exahype::records::ADERDGCellDescription::DeaugmentingRequested);
            break;
          default:
            break;
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

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

exahype::solvers::Solver::AugmentationControl
exahype::mappings::MarkingForAugmentation::augmentationCriterion(
    const int solverNumber,
    const exahype::records::ADERDGCellDescription::Type type,
    const int level,
    const tarch::la::Vector<THREE_POWER_D, int>& neighbourCellDescriptionIndices) const {
  // left,right,front,back,(front,back)
#if DIMENSIONS == 2
  constexpr int neighbourPositions[4] = {3,5,1,7};
#else
  constexpr int neighbourPositions[6] = {12,14,10,16,4,22};
#endif
  bool nextToAncestor = false;
  bool nextToCell     = false;

  for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
    const int neighbourCellDescriptionIndex =
        neighbourCellDescriptionIndices[neighbourPositions[i]];
    if (DataHeap::getInstance().isValidIndex(neighbourCellDescriptionIndex)) {
      // NOTE: The order of the following two for loops is important.
      // First, check if we are next to a ancestor.
      for (std::vector<exahype::records::ADERDGCellDescription>::
          iterator pNeighbour = ADERDGCellDescriptionHeap::getInstance().getData(
              neighbourCellDescriptionIndex).begin();
          pNeighbour != ADERDGCellDescriptionHeap::getInstance().getData(
              neighbourCellDescriptionIndex).end();
          ++pNeighbour) {
        if (pNeighbour->getSolverNumber()==solverNumber &&
            pNeighbour->getLevel()==level) {
          switch (pNeighbour->getType()) { // What type has the neighbour?
          case exahype::records::ADERDGCellDescription::Ancestor:
          case exahype::records::ADERDGCellDescription::EmptyAncestor:
            switch (type) { // What type am I?
            case exahype::records::ADERDGCellDescription::Cell:
            case exahype::records::ADERDGCellDescription::EmptyDescendant:
            case exahype::records::ADERDGCellDescription::Descendant:
              nextToAncestor = true;
              break;
            default:
              break;
            }
            break;
          case exahype::records::ADERDGCellDescription::Cell:
            switch (type) { // What type am I?
            case exahype::records::ADERDGCellDescription::EmptyAncestor:
            case exahype::records::ADERDGCellDescription::Ancestor:
            case exahype::records::ADERDGCellDescription::EmptyDescendant:
            case exahype::records::ADERDGCellDescription::Descendant:
              nextToCell = true;
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
    }
  }

  // NOTE: The order is important below
  if (nextToCell && nextToAncestor) {
    return exahype::solvers::Solver::NextToCellAndAncestor;
  }
  if (nextToAncestor) {
    return exahype::solvers::Solver::NextToAncestor;
  }
  if (nextToCell) {
    return exahype::solvers::Solver::NextToCell;
  }
  // Erase otherwise.
  return exahype::solvers::Solver::Unncessary;
}

void exahype::mappings::MarkingForAugmentation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::beginIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::endIteration(exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::MarkingForAugmentation::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  logTraceInWith2Arguments( "descend(...)", coarseGridCell.toString(), coarseGridVerticesEnumerator.toString() );

//  bool unsetAugmentation = false;
//  bool eraseChildren     = true;
//
//  if (ADERDGCellDescriptionHeap::getInstance().
//      isValidIndex(coarseGridCell.getADERDGCellDescriptionsIndex())) {
//    for (std::vector<exahype::records::ADERDGCellDescription>::
//        iterator pCoarse = ADERDGCellDescriptionHeap::getInstance().getData(
//            coarseGridCell.getADERDGCellDescriptionsIndex()).begin();
//        pCoarse != ADERDGCellDescriptionHeap::getInstance().getData(
//            coarseGridCell.getADERDGCellDescriptionsIndex()).end();
//        ++pCoarse) {
//      switch (pCoarse->getType()) {
//      case exahype::records::ADERDGCellDescription::Cell:
//      case exahype::records::ADERDGCellDescription::Descendant:
//      case exahype::records::ADERDGCellDescription::EmptyDescendant:
//        switch (pCoarse->getRefinementEvent()) {
//        case exahype::records::ADERDGCellDescription::Augmenting:
//        case exahype::records::ADERDGCellDescription::DeaugmentingRequested:
//        case exahype::records::ADERDGCellDescription::None:
//          unsetAugmentation = false;
//          eraseChildren     = true;
//
//          // 1. Loop over fine grid cells. Determine if erasing of
//          // fine grid cells is necessary and if coarse grid cell
//          // should be augmented or not.
//          dfor3(k)
//          if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(fineGridCells[kScalar].
//              getADERDGCellDescriptionsIndex())) {
//            for (std::vector<exahype::records::ADERDGCellDescription>::
//                iterator pFine = ADERDGCellDescriptionHeap::getInstance().
//                getData(fineGridCells[kScalar].
//                    getADERDGCellDescriptionsIndex()).begin();
//                pFine != ADERDGCellDescriptionHeap::getInstance().
//                    getData(fineGridCells[kScalar].
//                        getADERDGCellDescriptionsIndex()).end();
//                ++pFine) {
//              if (pCoarse->getSolverNumber()==pFine->getSolverNumber()) {
//                assertion1(pFine->getType()==
//                    exahype::records::ADERDGCellDescription::Descendant
//                    ||
//                    exahype::records::ADERDGCellDescription::EmptyDescendant,
//                    toString(fineGridVerticesEnumerator.getCellFlags()));
//                unsetAugmentation = true;
//                eraseChildren = eraseChildren && pFine->getRefinementEvent()==
//                    exahype::records::ADERDGCellDescription::DeaugmentingRequested;
//              }
//            }
//          } else {
//            eraseChildren = false;
//          }
//          enddforx
//
//          // 2. Untrigger the augmenting and deaugmenting request events
//          // if coarse grid cell description already has descendants
//          if (unsetAugmentation) {
//            assertion1(pCoarse->getRefinementEvent()==exahype::records::ADERDGCellDescription::Augmenting ||
//                pCoarse->getRefinementEvent()==exahype::records::ADERDGCellDescription::DeaugmentingRequested
//                ,toString(fineGridVerticesEnumerator.getCellFlags()));
//            pCoarse->setRefinementEvent(exahype::records::ADERDGCellDescription::None);
//          }
//
//          // 3. Trigger erasing event on the fine grid cell description.
//          if (eraseChildren) {
//            assertion1(pCoarse->getRefinementEvent()==exahype::records::ADERDGCellDescription::None,
//                toString(fineGridVerticesEnumerator.getCellFlags()));
//            pCoarse->setRefinementEvent(exahype::records::ADERDGCellDescription::ErasingChildren);
//
//            dfor3(k)
//            for (std::vector<exahype::records::ADERDGCellDescription>::
//                iterator pFine = ADERDGCellDescriptionHeap::getInstance().
//                getData(fineGridCells[kScalar].
//                    getADERDGCellDescriptionsIndex()).begin();
//                pFine != ADERDGCellDescriptionHeap::getInstance().
//                    getData(fineGridCells[kScalar].
//                        getADERDGCellDescriptionsIndex()).end();
//                ++pFine) {
//              if (pCoarse->getSolverNumber()==pFine->getSolverNumber()) {
//                pFine->setType(exahype::records::ADERDGCellDescription::Erased);
//                pFine->setRefinementEvent(exahype::records::ADERDGCellDescription::Erasing);
//              }
//            }
//            enddforx
//          }
//          break;
//        default:
//          break;
//        }
//        break;
//        default:
//          break;
//      }
//    }
//  }

  logTraceOut( "descend(...)" );
}

void exahype::mappings::MarkingForAugmentation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
