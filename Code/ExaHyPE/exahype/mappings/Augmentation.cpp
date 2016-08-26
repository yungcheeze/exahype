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
 
#include "exahype/mappings/Augmentation.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "exahype/solvers/Solver.h"
#include "kernels/KernelCalls.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

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
  return peano::MappingSpecification(peano::MappingSpecification::Nop,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Augmentation::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Augmentation::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Augmentation::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Augmentation::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Augmentation::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,
                                     peano::MappingSpecification::Serial);
}

tarch::logging::Log exahype::mappings::Augmentation::_log(
    "exahype::mappings::Augmentation");

void exahype::mappings::Augmentation::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  if (coarseGridCell.isInitialised()) {
    for (auto& pCoarse : ADERDGCellDescriptionHeap::getInstance().getData(
        coarseGridCell.getCellDescriptionsIndex())) {
      bool solverNotFound = true;

      switch (pCoarse.getRefinementEvent()) {
        case exahype::records::ADERDGCellDescription::DeaugmentingRequested:
          /*
           * Reset the deaugmenting request if the coarse grid Descendant has children
           * (of type Descendant).
           *
           * Rationale:
           * We cannot erase a coarse grid cell that has children (of type Descendant)
           * before erasing the children.
           *
           * Note:
           * A more sophisticated procedure has to performed for the refinement event
           * AugmentationRequested. We need to use the taversal's descend event to handle
           * this event.
           * We thus do not rely on fineGridCell.isRefined() in the previous enterCell event
           * to check if we need to reset the deaugmenting request.
           *
           */
          assertion1(pCoarse.getType()==exahype::records::ADERDGCellDescription::EmptyDescendant ||
                     pCoarse.getType()==exahype::records::ADERDGCellDescription::Descendant,pCoarse.toString());
          pCoarse.setRefinementEvent(exahype::records::ADERDGCellDescription::None);
          break;
        case exahype::records::ADERDGCellDescription::AugmentingRequested:
        case exahype::records::ADERDGCellDescription::Augmenting:
          assertion1(pCoarse.getType()==exahype::records::ADERDGCellDescription::Cell ||
                     pCoarse.getType()==exahype::records::ADERDGCellDescription::EmptyDescendant ||
                     pCoarse.getType()==exahype::records::ADERDGCellDescription::Descendant,pCoarse.toString());

          solverNotFound = true;
          if (fineGridCell.isInitialised()) {
            assertion( fineGridCell.isInitialised() );
            for (auto& pFine : ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.
                                                                                getCellDescriptionsIndex())) {
              if (pCoarse.getSolverNumber() == pFine.getSolverNumber()) {
                assertion1(pFine.getType()==exahype::records::ADERDGCellDescription::Descendant ||
                           pFine.getType()==exahype::records::ADERDGCellDescription::EmptyDescendant,
                           pFine.toString());
                solverNotFound = false;
              }
            }
          }
          // Initialise uninitialised cells.
          if (solverNotFound) {
            fineGridCell.addNewCellDescription(
                pCoarse.getSolverNumber(),
                exahype::records::ADERDGCellDescription::EmptyDescendant,
                exahype::records::ADERDGCellDescription::None,
                fineGridVerticesEnumerator.getLevel(),
                coarseGridCell.getCellDescriptionsIndex(),
                fineGridVerticesEnumerator.getCellSize(),
                // We pass the lower left corner of the cell as offset.
                fineGridVerticesEnumerator.getVertexPosition());
            pCoarse.setRefinementEvent(exahype::records::ADERDGCellDescription::Augmenting);
          } else {
            pCoarse.setRefinementEvent(exahype::records::ADERDGCellDescription::None);
          }
          break;
        default:
          break;
      }
    }
  }
}

void exahype::mappings::Augmentation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  if (fineGridCell.isInitialised()) {
    for (auto& pFine : ADERDGCellDescriptionHeap::getInstance().getData(
        fineGridCell.getCellDescriptionsIndex())) {

      switch (pFine.getRefinementEvent()) {
        case exahype::records::ADERDGCellDescription::Augmenting:
          pFine.setRefinementEvent(exahype::records::ADERDGCellDescription::None);
          break;
        default:
          break;
      }
    }
  }
}

void exahype::mappings::Augmentation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  logTraceInWith2Arguments("ascend(...)", coarseGridCell.toString(),
                           coarseGridVerticesEnumerator.toString());

  if (coarseGridCell.isInitialised()) {
    for (auto& pCoarse : ADERDGCellDescriptionHeap::getInstance().getData(
             coarseGridCell.getCellDescriptionsIndex())) {
      bool eraseChildren = true;

      switch (pCoarse.getType()) {
        case exahype::records::ADERDGCellDescription::Cell:
        case exahype::records::ADERDGCellDescription::EmptyDescendant:
        case exahype::records::ADERDGCellDescription::Descendant:
          switch (pCoarse.getRefinementEvent()) {
            case exahype::records::ADERDGCellDescription::None:
              eraseChildren = true;
              dfor3(k)
                if (fineGridCells[kScalar].isInitialised()) {
                  for (auto& pFine : ADERDGCellDescriptionHeap::getInstance().getData(
                          fineGridCells[kScalar].getCellDescriptionsIndex())) {
                    if (pCoarse.getSolverNumber() == pFine.getSolverNumber()) {
                      eraseChildren = eraseChildren &&
                          pFine.getRefinementEvent() ==
                              exahype::records::ADERDGCellDescription::DeaugmentingRequested;
                    }
                  }
                }
              enddforx

              if (eraseChildren) {
                dfor3(k)
                  if (fineGridCells[kScalar].isInitialised()) {
                    auto pFine = ADERDGCellDescriptionHeap::getInstance().
                        getData(fineGridCells[kScalar].getCellDescriptionsIndex()).begin();
                    while (pFine != ADERDGCellDescriptionHeap::getInstance().
                        getData(fineGridCells[kScalar].getCellDescriptionsIndex()).end()) {
                      if (pCoarse.getSolverNumber()==pFine->getSolverNumber()) {
                        assertion1(pFine->getRefinementEvent() == exahype::records::ADERDGCellDescription::DeaugmentingRequested,
                                   pFine->toString());
                        pFine->setType(exahype::records::ADERDGCellDescription::Erased);
                        pFine->setRefinementEvent(exahype::records::ADERDGCellDescription::Erasing);
                        fineGridCells[kScalar].ensureNoUnnecessaryMemoryIsAllocated(pFine->getSolverNumber());
                        pFine = ADERDGCellDescriptionHeap::getInstance().
                            getData(fineGridCells[kScalar].getCellDescriptionsIndex()).erase(pFine);
                      } else {
                        ++pFine;
                      }
                    }

                    if (ADERDGCellDescriptionHeap::getInstance().getData(fineGridCells[kScalar].getCellDescriptionsIndex()).empty()) {
                      fineGridCells[kScalar].shutdownMetaData();
                    }
                  }
                enddforx

              // reset if not all children requested deaugmenting
              } else {
                dfor3(k)
                  if (fineGridCells[kScalar].isInitialised()) {
                  for (auto& pFine : ADERDGCellDescriptionHeap::getInstance().getData(fineGridCells[kScalar]
                                                               .getCellDescriptionsIndex())) {
                    if (pCoarse.getSolverNumber()==pFine.getSolverNumber()) {
                      if (pFine.getRefinementEvent()==exahype::records::ADERDGCellDescription::DeaugmentingRequested) {
                        pFine.setRefinementEvent(exahype::records::ADERDGCellDescription::None);
                      }
                    }
                  }
                }
                enddforx
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



//
// Below all methods are nop.
//
// ====================================




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

exahype::mappings::Augmentation::Augmentation() {
  // do nothing
}

exahype::mappings::Augmentation::~Augmentation() {
  // do nothing
}

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
 // do nothing
}
