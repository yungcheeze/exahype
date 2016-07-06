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
 
#include "exahype/mappings/Refinement.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "exahype/solvers/Solver.h"
#include "kernels/KernelCalls.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "peano/datatraversal/autotuning/Oracle.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::Refinement::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
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
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,
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
  bool eraseFineGridVertex  = true;
  /*
   * Loop over the 2^d adjacent cells and check if refinement
   * or easing is necessary.
   */
  dfor2(c)
    if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
        cellDescriptions[cScalar])) {
      eraseFineGridVertex = false;

      for (auto& pFine : ADERDGCellDescriptionHeap::getInstance()
                         .getData(cellDescriptions[cScalar])) {
        assertion3(static_cast<unsigned int>(pFine.getSolverNumber()) <
            solvers::RegisteredSolvers.size(),
            pFine.getSolverNumber(), solvers::RegisteredSolvers.size(),
            toString(coarseGridVerticesEnumerator.getCellFlags()));

        switch (pFine.getType()) {
        case exahype::records::ADERDGCellDescription::Cell:
          switch (pFine.getRefinementEvent()) {
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
    } else {
      if (cellDescriptions[cScalar]!=exahype::Cell::ErasedCellDescriptionsIndex) {
        eraseFineGridVertex = false;
      }
    }
  enddforx

  // NOTE: Order might be important.
  if (eraseFineGridVertex) {
    if (!fineGridVertex.isHangingNode() && !fineGridVertex.isRefinedOrRefining()) {
      fineGridVertex.erase();
    }
  }

  if (refineFineGridVertex) {
    if (!fineGridVertex.isHangingNode() && !fineGridVertex.isRefinedOrRefining()) { // todo discuss with Tobias
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
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                             fineGridVerticesEnumerator.toString(),
                             coarseGridCell, fineGridPositionOfCell);

    if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(fineGridCell.getADERDGCellDescriptionsIndex())) {
      if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(coarseGridCell.getADERDGCellDescriptionsIndex())) {
        // Please use a different UserDefined per mapping/event.
        const int numberOfADERDGCellDescriptions = static_cast<int>(
            ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.getADERDGCellDescriptionsIndex()).size());
        const peano::datatraversal::autotuning::MethodTrace methodTrace =peano::datatraversal::autotuning::UserDefined1;
        const int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
                numberOfADERDGCellDescriptions, methodTrace);
        pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
          records::ADERDGCellDescription& pFine = fineGridCell.getADERDGCellDescription(i);
          for (auto& pCoarse : ADERDGCellDescriptionHeap::getInstance().
              getData(coarseGridCell.getADERDGCellDescriptionsIndex())) {
            if (pFine.getSolverNumber() == pCoarse.getSolverNumber()) {
              switch (pFine.getRefinementEvent()) {
                case exahype::records::ADERDGCellDescription::Prolongating:
                  assertion3(pFine.getParentIndex() ==coarseGridCell.getADERDGCellDescriptionsIndex(),
                          pFine.getParentIndex(),coarseGridCell.getADERDGCellDescriptionsIndex(),coarseGridCell.toString());
                  assertion(pFine.getType()==exahype::records::ADERDGCellDescription::Cell);
                  assertion1(pCoarse.getType() ==exahype::records::ADERDGCellDescription::Cell,
                             toString(fineGridVerticesEnumerator.getCellFlags()));
                  assertion3(pCoarse.getRefinementEvent() == exahype::records::ADERDGCellDescription::Refining,
                             pCoarse.getType(), pCoarse.getRefinementEvent(),toString(fineGridVerticesEnumerator.getCellFlags()));
                  prolongateVolumeData(
                      pFine,
                      pCoarse,
                      fineGridPositionOfCell);
                  pFine.setRefinementEvent(exahype::records::ADERDGCellDescription::None);
                  break;
                default:
                  break;
              }
           }
         }
        endpfor peano::datatraversal::autotuning::Oracle::getInstance()
            .parallelSectionHasTerminated(methodTrace);
      }
    }
    logTraceOutWith1Argument("enterCell(...)", fineGridCell);
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
  logTraceInWith2Arguments("ascend(...)", coarseGridCell.toString(),
                           coarseGridVerticesEnumerator.toString());

  if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
          coarseGridCell.getADERDGCellDescriptionsIndex())) {
    for (auto& pCoarse : ADERDGCellDescriptionHeap::getInstance().getData(
             coarseGridCell.getADERDGCellDescriptionsIndex())) {
      bool eraseChildren = true;

      switch (pCoarse.getType()) {
        case exahype::records::ADERDGCellDescription::Cell:
          switch (pCoarse.getRefinementEvent()) {
            case exahype::records::ADERDGCellDescription::Refining:
              pCoarse.setType(exahype::records::ADERDGCellDescription::Ancestor);
              pCoarse.setRefinementEvent(exahype::records::ADERDGCellDescription::None);
              coarseGridCell.ensureNoUnnecessaryMemoryIsAllocated(pCoarse.getSolverNumber());
              break;
            default:
              break;
          }
          break;
        case exahype::records::ADERDGCellDescription::EmptyAncestor:
        case exahype::records::ADERDGCellDescription::Ancestor:
          switch (pCoarse.getRefinementEvent()) {
            case exahype::records::ADERDGCellDescription::None:
              eraseChildren = true;

              dfor3(k)
                assertion(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                    fineGridCells[kScalar].getADERDGCellDescriptionsIndex()));
                for (auto& pFine : ADERDGCellDescriptionHeap::getInstance()
                    .getData(fineGridCells[kScalar]
                    .getADERDGCellDescriptionsIndex())) {
                  if (pCoarse.getSolverNumber() == pFine.getSolverNumber()) {
                    eraseChildren = eraseChildren &&
                        pFine.getRefinementEvent() ==
                            exahype::records::ADERDGCellDescription::
                            ErasingRequested;
                  }
                }
              enddforx

              if (eraseChildren) {
                coarseGridCell.ensureNecessaryMemoryIsAllocated(pCoarse.getSolverNumber());

                dfor3(k)
                  auto pFine = ADERDGCellDescriptionHeap::getInstance().
                                  getData(fineGridCells[kScalar].getADERDGCellDescriptionsIndex()).begin();
                  while (pFine != ADERDGCellDescriptionHeap::getInstance().
                      getData(fineGridCells[kScalar].getADERDGCellDescriptionsIndex()).end()) {
                    if (pCoarse.getSolverNumber() ==
                        pFine->getSolverNumber()) {
                      assertion1(pFine->getRefinementEvent() == exahype::records::ADERDGCellDescription::ErasingRequested,
                              toString(fineGridVerticesEnumerator.getCellFlags()));
                      exahype::Cell::SubcellPosition subcellPosition =
                          fineGridCells[kScalar].computeSubcellPositionOfCellOrAncestor(*pFine);
                      restrictVolumeData(pCoarse,(*pFine),subcellPosition.subcellIndex);
                      pFine->setType(exahype::records::ADERDGCellDescription::Erased);
                      pFine->setRefinementEvent(exahype::records::ADERDGCellDescription::Erasing);
                      fineGridCells[kScalar].ensureNoUnnecessaryMemoryIsAllocated(pFine->getSolverNumber());
                      pFine = ADERDGCellDescriptionHeap::getInstance().
                         getData(fineGridCells[kScalar].getADERDGCellDescriptionsIndex()).erase(pFine);
                    } else {
                      ++pFine;
                    }
                  }

                  if (ADERDGCellDescriptionHeap::getInstance().getData(fineGridCells[kScalar].getADERDGCellDescriptionsIndex()).empty()) {
                    fineGridCells[kScalar].shutdownMetaData();
                  }
                enddforx

              // reset if not all children requested erasing
              } else {
                dfor3(k)
                  for (auto& pFine : ADERDGCellDescriptionHeap::getInstance().getData(fineGridCells[kScalar]
                     .getADERDGCellDescriptionsIndex())) {
                    if (pCoarse.getSolverNumber() ==
                        pFine.getSolverNumber()) {
                      if (pFine.getRefinementEvent()==exahype::records::ADERDGCellDescription::ErasingRequested) {
                        pFine.setRefinementEvent(exahype::records::ADERDGCellDescription::None);
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
    for (auto& pCoarse : ADERDGCellDescriptionHeap::getInstance().getData(
             coarseGridCell.getADERDGCellDescriptionsIndex())) {
      bool solverNotFound = true;

      switch (pCoarse.getType()) {
        case exahype::records::ADERDGCellDescription::Cell:
          switch (pCoarse.getRefinementEvent()) {
            case exahype::records::ADERDGCellDescription::ErasingRequested:
              /*
               * Reset the erasing request if the coarse grid Cell has children
               * (of type Descendant).
               *
               * Rationale:
               * We cannot erase a Cell that has children (of type Descendant)
               * before erasing the children.
               *
               * Note:
               * A more sophisticated procedure has to performed for the refinement event
               * RefiningRequested. We need to use the taversal's descend event to handle
               * this event.
               * We thus do not rely on fineGridCell.isRefined() in the previous enterCell event
               * to check if we need to reset the erasing request.
               *
               */
              pCoarse.setRefinementEvent(exahype::records::ADERDGCellDescription::None);
              break;
            case exahype::records::ADERDGCellDescription::RefiningRequested:
              /*
               * If the coarse grid cell has children of type Descendant,
               * we change the type of the children to Cell.
               * We furthermore set the
               */
              dfor3(k)
                solverNotFound = true;
                if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                    fineGridCells[kScalar].getADERDGCellDescriptionsIndex())) {
                  for (auto& pFine : ADERDGCellDescriptionHeap::getInstance().getData(fineGridCells[kScalar].
                                                                                      getADERDGCellDescriptionsIndex())) {
                    if (pCoarse.getSolverNumber() == pFine.getSolverNumber()) {
                      assertion1(pFine.getType()==exahype::records::ADERDGCellDescription::Descendant ||
                                 pFine.getType()==exahype::records::ADERDGCellDescription::EmptyDescendant,pFine.toString());
                      pFine.setType(exahype::records::ADERDGCellDescription::Cell);
                      pFine.setRefinementEvent(exahype::records::ADERDGCellDescription::Prolongating);
                      fineGridCells[kScalar].ensureNecessaryMemoryIsAllocated(pFine.getSolverNumber());
                      pCoarse.setRefinementEvent(exahype::records::ADERDGCellDescription::Refining);
                      solverNotFound = false;
                    }
                  }
                }
                // Initialise uninitialised cells.
                if (solverNotFound) {
                  fineGridCells[kScalar].addNewCellDescription(
                      pCoarse.getSolverNumber(),
                      exahype::records::ADERDGCellDescription::Cell,
                      exahype::records::ADERDGCellDescription::Prolongating,
                      fineGridVerticesEnumerator.getLevel(),
                      coarseGridCell.getADERDGCellDescriptionsIndex(),
                      fineGridVerticesEnumerator.getCellSize(),
                      // We pass the lower left corner of the cell as offset.
                      fineGridVerticesEnumerator.getVertexPosition());
                  fineGridCells[kScalar].ensureNecessaryMemoryIsAllocated(pCoarse.getSolverNumber());

                  pCoarse.setRefinementEvent(exahype::records::ADERDGCellDescription::Refining);
                }
              enddforx
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

  logTraceOut("descend(...)");
}

void exahype::mappings::Refinement::prolongateVolumeData(
    exahype::records::ADERDGCellDescription& p,
    const exahype::records::ADERDGCellDescription& pCoarse,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) const {
  const int levelFine = p.getLevel();
  const int levelCoarse = pCoarse.getLevel();
  assertion(levelCoarse < levelFine);

  double* luhFine   = DataHeap::getInstance().getData(p.getSolution()).data();
  double* luhCoarse = DataHeap::getInstance().getData(pCoarse.getSolution()).data();

  exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[p.getSolverNumber()];
  solver->volumeUnknownsProlongation(luhFine, luhCoarse, levelCoarse, levelFine,
                                     subcellIndex);
}

void exahype::mappings::Refinement::restrictVolumeData(
    exahype::records::ADERDGCellDescription& pCoarse,
    const exahype::records::ADERDGCellDescription& p,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) const {
  const int levelFine = p.getLevel();
  const int levelCoarse = pCoarse.getLevel();
  assertion(levelCoarse < levelFine);

  double* luhFine   = DataHeap::getInstance().getData(p.getSolution()).data();
  double* luhCoarse = DataHeap::getInstance().getData(pCoarse.getSolution()).data();

  exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[p.getSolverNumber()];
  solver->volumeUnknownsRestriction(luhCoarse, luhFine, levelCoarse, levelFine,
                                    subcellIndex);
}
