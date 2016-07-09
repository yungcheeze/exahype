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
  return peano::MappingSpecification(peano::MappingSpecification::Nop,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Refinement::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Refinement::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Refinement::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Refinement::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,
                                     peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::Refinement::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,
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
  // do nothing
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
  if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      coarseGridCell.getADERDGCellDescriptionsIndex())) {
    for (auto& pCoarse : ADERDGCellDescriptionHeap::getInstance().getData(
        coarseGridCell.getADERDGCellDescriptionsIndex())) {
      switch (pCoarse.getType()) {
        case exahype::records::ADERDGCellDescription::Cell:
          switch (pCoarse.getRefinementEvent()) {
            case exahype::records::ADERDGCellDescription::ErasingRequested:
              /*
               * Change the erasing request to a change to descendant request if the coarse grid Cell
               * has children (of type Descendant).
               *
               * Rationale:
               * We cannot directly erase a Cell that has children (of type Descendant).
               */
              pCoarse.setRefinementEvent(exahype::records::ADERDGCellDescription::ChangeToDescendantRequested);
              break;
            case exahype::records::ADERDGCellDescription::RefiningRequested:
            case exahype::records::ADERDGCellDescription::Refining:
              /*
               * If the coarse grid cell has children of type Descendant,
               * we change the type of the children to Cell.
               * We furthermore set the
               */
              if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                  fineGridCell.getADERDGCellDescriptionsIndex())) {
                for (auto& pFine : ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.getADERDGCellDescriptionsIndex())) {
                  if (pCoarse.getSolverNumber() == pFine.getSolverNumber()) {
                    assertion1(pFine.getType()==exahype::records::ADERDGCellDescription::Descendant ||
                               pFine.getType()==exahype::records::ADERDGCellDescription::EmptyDescendant,pFine.toString());
                    pFine.setType(exahype::records::ADERDGCellDescription::Cell);
                    pFine.setRefinementEvent(exahype::records::ADERDGCellDescription::None);
                    fineGridCell.ensureNecessaryMemoryIsAllocated(pFine.getSolverNumber());
                    prolongateVolumeData(
                        pFine,
                        pCoarse,
                        fineGridPositionOfCell);
                    pCoarse.setRefinementEvent(exahype::records::ADERDGCellDescription::Refining);
                  }
                }
              } else  {
                fineGridCell.addNewCellDescription(
                    pCoarse.getSolverNumber(),
                    exahype::records::ADERDGCellDescription::Cell,
                    exahype::records::ADERDGCellDescription::None,
                    fineGridVerticesEnumerator.getLevel(),
                    coarseGridCell.getADERDGCellDescriptionsIndex(),
                    fineGridVerticesEnumerator.getCellSize(),
                    // We pass the lower left corner of the cell as offset.
                    fineGridVerticesEnumerator.getVertexPosition());
                fineGridCell.ensureNecessaryMemoryIsAllocated(pCoarse.getSolverNumber());

                for (auto& pFine : ADERDGCellDescriptionHeap::getInstance().
                    getData(fineGridCell.getADERDGCellDescriptionsIndex())) {
                  if (pCoarse.getSolverNumber() == pFine.getSolverNumber()) {
                    prolongateVolumeData(
                        pFine,
                        pCoarse,
                        fineGridPositionOfCell);
                  }
                }
                pCoarse.setRefinementEvent(exahype::records::ADERDGCellDescription::Refining);
              }
              break;
            case exahype::records::ADERDGCellDescription::None:
              if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                  fineGridCell.getADERDGCellDescriptionsIndex())) {
                /**
                 * In this case the parent did not trigger refinement.
                 */

              }
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

void exahype::mappings::Refinement::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
            fineGridCell.getADERDGCellDescriptionsIndex())) {
      for (auto& pFine : ADERDGCellDescriptionHeap::getInstance().getData(
               fineGridCell.getADERDGCellDescriptionsIndex())) {
        switch (pFine.getRefinementEvent()) {
          case exahype::records::ADERDGCellDescription::Refining:
            assertion1(pFine.getType()==exahype::records::ADERDGCellDescription::Cell,pFine.toString());
            pFine.setType(exahype::records::ADERDGCellDescription::EmptyAncestor);
            pFine.setRefinementEvent(exahype::records::ADERDGCellDescription::None);
            fineGridCell.ensureNoUnnecessaryMemoryIsAllocated(pFine.getSolverNumber());
            break;
          default:
            break;
        }
      }
    }
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
      bool eraseChildren               = true;
      bool changeChildrenToDescendants = false;

      switch (pCoarse.getType()) {
        case exahype::records::ADERDGCellDescription::EmptyAncestor:
        case exahype::records::ADERDGCellDescription::Ancestor:
          eraseChildren = true;

          dfor3(k)
            assertion(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                fineGridCells[kScalar].getADERDGCellDescriptionsIndex()));
            for (auto& pFine : ADERDGCellDescriptionHeap::getInstance()
                .getData(fineGridCells[kScalar]
                .getADERDGCellDescriptionsIndex())) {
              if (pCoarse.getSolverNumber() == pFine.getSolverNumber()) {
                eraseChildren = eraseChildren &&
                    (pFine.getRefinementEvent() == exahype::records::ADERDGCellDescription::ErasingRequested ||
                     pFine.getRefinementEvent() == exahype::records::ADERDGCellDescription::ChangeToDescendantRequested);
                changeChildrenToDescendants =
                    changeChildrenToDescendants ||
                    pFine.getRefinementEvent() == exahype::records::ADERDGCellDescription::ChangeToDescendantRequested;
              }
            }
          enddforx

          if (eraseChildren) {
            pCoarse.setType(exahype::records::ADERDGCellDescription::Cell);
            coarseGridCell.ensureNecessaryMemoryIsAllocated(pCoarse.getSolverNumber());

            dfor3(k)
              auto pFine = ADERDGCellDescriptionHeap::getInstance().
                              getData(fineGridCells[kScalar].getADERDGCellDescriptionsIndex()).begin();
              while (pFine != ADERDGCellDescriptionHeap::getInstance().
                  getData(fineGridCells[kScalar].getADERDGCellDescriptionsIndex()).end()) {
                if (pCoarse.getSolverNumber() == pFine->getSolverNumber()) {
                  exahype::Cell::SubcellPosition subcellPosition =
                      fineGridCells[kScalar].computeSubcellPositionOfCellOrAncestor(*pFine);
                  restrictVolumeData(pCoarse,(*pFine),subcellPosition.subcellIndex);
                  if (changeChildrenToDescendants) {
                    pFine->setType(exahype::records::ADERDGCellDescription::EmptyDescendant);
                    pFine->setRefinementEvent(exahype::records::ADERDGCellDescription::None);
                    fineGridCells[kScalar].ensureNoUnnecessaryMemoryIsAllocated(pFine->getSolverNumber());
                    ++pFine;
                  } else {
                    pFine->setType(exahype::records::ADERDGCellDescription::Erased);
                    fineGridCells[kScalar].ensureNoUnnecessaryMemoryIsAllocated(pFine->getSolverNumber());
                    pFine = ADERDGCellDescriptionHeap::getInstance().
                       getData(fineGridCells[kScalar].getADERDGCellDescriptionsIndex()).erase(pFine);
                  }
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
  // do nothing
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
