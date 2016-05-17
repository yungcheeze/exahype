#include "exahype/mappings/VolumeUnknownsProjection.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "exahype/solvers/Solver.h"
#include "kernels/KernelCalls.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/VertexOperations.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::VolumeUnknownsProjection::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification exahype::mappings::VolumeUnknownsProjection::
    touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification exahype::mappings::VolumeUnknownsProjection::
    touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::VolumeUnknownsProjection::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::VolumeUnknownsProjection::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::VolumeUnknownsProjection::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::VolumeUnknownsProjection::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::VolumeUnknownsProjection::_log(
    "exahype::mappings::VolumeUnknownsProjection");

exahype::mappings::VolumeUnknownsProjection::VolumeUnknownsProjection() {
  // do nothing
}

exahype::mappings::VolumeUnknownsProjection::~VolumeUnknownsProjection() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::VolumeUnknownsProjection::VolumeUnknownsProjection(
    const VolumeUnknownsProjection& masterThread) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::mergeWithWorkerThread(
    const VolumeUnknownsProjection& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::VolumeUnknownsProjection::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::VolumeUnknownsProjection::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::VolumeUnknownsProjection::prepareSendToWorker(
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

void exahype::mappings::VolumeUnknownsProjection::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::mergeWithMaster(
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

void exahype::mappings::VolumeUnknownsProjection::receiveDataFromMaster(
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

void exahype::mappings::VolumeUnknownsProjection::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::VolumeUnknownsProjection::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::enterCell(
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
          fineGridCell.getADERDGCellDescriptionsIndex())) {
    if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
            coarseGridCell.getADERDGCellDescriptionsIndex())) {
      // please use a different UserDefined per mapping/event
      const int numberOfADERDGCellDescriptions = static_cast<int>(
          ADERDGCellDescriptionHeap::getInstance()
              .getData(fineGridCell.getADERDGCellDescriptionsIndex())
              .size());
      const peano::datatraversal::autotuning::MethodTrace methodTrace =
          peano::datatraversal::autotuning::UserDefined1;
      const int grainSize =
          peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
              numberOfADERDGCellDescriptions, methodTrace);

      // clang-format off
      pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
          records::ADERDGCellDescription& pFine =
              fineGridCell.getADERDGCellDescription(i);

        for (std::vector<exahype::records::ADERDGCellDescription>::iterator
                 pCoarse = ADERDGCellDescriptionHeap::getInstance()
                               .getData(coarseGridCell.getADERDGCellDescriptionsIndex())
                               .begin();
             pCoarse !=
             ADERDGCellDescriptionHeap::getInstance()
                 .getData(coarseGridCell.getADERDGCellDescriptionsIndex())
                 .end();
             ++pCoarse) {
          if (pFine.getSolverNumber() == pCoarse->getSolverNumber()) {
            switch (pCoarse->getRefinementEvent()) {
              case exahype::records::ADERDGCellDescription::Refining:
                switch (pFine.getRefinementEvent()) {
                  case exahype::records::ADERDGCellDescription::Prolongating:
                    assertion3(
                        pFine.getParentIndex() == coarseGridCell.getADERDGCellDescriptionsIndex(),
                        pFine.getParentIndex(), coarseGridCell.getADERDGCellDescriptionsIndex(),
                        coarseGridCell.toString());
                    assertion(pFine.getType() ==
                              exahype::records::ADERDGCellDescription::Cell);
                    assertion1(
                        pCoarse->getType() ==
                            exahype::records::ADERDGCellDescription::Cell,
                        toString(fineGridVerticesEnumerator.getCellFlags()));
                    assertion3(
                        pCoarse->getRefinementEvent() ==
                            exahype::records::ADERDGCellDescription::Refining,
                        pCoarse->getType(), pCoarse->getRefinementEvent(),
                        toString(fineGridVerticesEnumerator.getCellFlags()));

                    //              prolongateVolumeData(
                    //                  p,
                    //                  *pCoarse,
                    //                  fineGridPositionOfCell);
                    pFine.setRefinementEvent(
                          exahype::records::ADERDGCellDescription::None);
                    break;
                  default:
                    break;
                }
                break;
              case exahype::records::ADERDGCellDescription::ErasingChildren:
                switch (pFine.getRefinementEvent()) {
                  case exahype::records::ADERDGCellDescription::Restricting:
                    if (pFine.getSolverNumber() == pCoarse->getSolverNumber()) {
                      assertion1(
                          pCoarse->getType() ==
                              exahype::records::ADERDGCellDescription::Cell,
                          toString(fineGridVerticesEnumerator.getCellFlags()));
                      assertion1(
                          pCoarse->getRefinementEvent() ==
                              exahype::records::ADERDGCellDescription::
                                  ErasingChildren,
                          toString(fineGridVerticesEnumerator.getCellFlags()));

                      //                  restrictVolumeData(p, *pCoarse,
                      //                  fineGridPositionOfCell);
                      pFine.setType(exahype::records::ADERDGCellDescription::Erased);
                      pFine.setRefinementEvent(exahype::records::ADERDGCellDescription::Erasing);
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
      endpfor peano::datatraversal::autotuning::Oracle::getInstance()
          .parallelSectionHasTerminated(methodTrace);
    }
  }

  // clang-format on
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::VolumeUnknownsProjection::prolongateVolumeData(
    const exahype::records::ADERDGCellDescription& p,
    const exahype::records::ADERDGCellDescription& pCoarse,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) const {
  const int levelFine = p.getLevel();
  const int levelCoarse = pCoarse.getLevel();
  assertion(levelCoarse < levelFine);

  double* luhFine = DataHeap::getInstance().getData(p.getSolution()).data();
  double* luhCoarse =
      DataHeap::getInstance().getData(pCoarse.getSolution()).data();

  exahype::solvers::Solver* solver =
      exahype::solvers::RegisteredSolvers[p.getSolverNumber()];
  solver->volumeUnknownsProlongation(luhFine, luhCoarse, levelCoarse, levelFine,
                                     subcellIndex);
}

void exahype::mappings::VolumeUnknownsProjection::restrictVolumeData(
    const exahype::records::ADERDGCellDescription& p,
    const exahype::records::ADERDGCellDescription& pCoarse,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) const {
  const int levelFine = p.getLevel();
  const int levelCoarse = pCoarse.getLevel();
  assertion(levelCoarse < levelFine);

  double* luhFine = DataHeap::getInstance().getData(p.getSolution()).data();
  double* luhCoarse =
      DataHeap::getInstance().getData(pCoarse.getSolution()).data();

  exahype::solvers::Solver* solver =
      exahype::solvers::RegisteredSolvers[p.getSolverNumber()];
  solver->volumeUnknownsRestriction(luhCoarse, luhFine, levelCoarse, levelFine,
                                    subcellIndex);
}

void exahype::mappings::VolumeUnknownsProjection::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::beginIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::endIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::ascend(
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
      bool refiningDone = true;

      switch (pCoarse->getType()) {
      case exahype::records::ADERDGCellDescription::Cell:
        switch (pCoarse->getRefinementEvent()) {
        case exahype::records::ADERDGCellDescription::Refining:
          refiningDone = true;
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
                refiningDone =
                    refiningDone &&
                    pFine->getRefinementEvent() ==
                        exahype::records::ADERDGCellDescription::None;
              }
            }
          } else {
            refiningDone = false;
          }
          enddforx

          if (refiningDone) {
            //                std::cout << "Refining done" << std::endl;

            // Rationale: It's more likely that a ancestor needs to hold data.
            pCoarse->setType(
                exahype::records::ADERDGCellDescription::EmptyAncestor);
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
