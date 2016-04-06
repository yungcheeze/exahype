#include "exahype/mappings/VirtualRefinement.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "kernels/KernelCalls.h"
#include "exahype/solvers/Solver.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::VirtualRefinement::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::
      SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::
      SendDataAndStateAfterLastTouchVertexLastTime,
      false);
}

peano::MappingSpecification
exahype::mappings::VirtualRefinement::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::VirtualRefinement::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::VirtualRefinement::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::VirtualRefinement::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::VirtualRefinement::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::VirtualRefinement::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::VirtualRefinement::_log(
    "exahype::mappings::VirtualRefinement");

exahype::mappings::VirtualRefinement::VirtualRefinement() {
  // do nothing
}

exahype::mappings::VirtualRefinement::~VirtualRefinement() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::VirtualRefinement::VirtualRefinement(const VirtualRefinement& masterThread) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::mergeWithWorkerThread(
    const VirtualRefinement& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::VirtualRefinement::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::VirtualRefinement::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::VirtualRefinement::prepareSendToWorker(
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

void exahype::mappings::VirtualRefinement::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::mergeWithMaster(
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

void exahype::mappings::VirtualRefinement::receiveDataFromMaster(
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

void exahype::mappings::VirtualRefinement::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::VirtualRefinement::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::VirtualRefinement::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);
  const tarch::la::Vector<THREE_POWER_D,int> neighbourCellDescriptionIndices =
      multiscalelinkedcell::getIndicesAroundCell(
          VertexOperations::readADERDGCellDescriptionsIndex(
              fineGridVerticesEnumerator,fineGridVertices));

  // If fine grid cell description does not exist, i.e. cell is new
  if (!ADERDGCellDescriptionHeap::getInstance().
      isValidIndex(fineGridCell.getADERDGCellDescriptionsIndex())) {
    assertion(ADERDGCellDescriptionHeap::getInstance().
           isValidIndex(coarseGridCell.getADERDGCellDescriptionsIndex()));
    int  solverNumber=0;
    for (std::vector<exahype::solvers::Solver*>::const_iterator p =
        exahype::solvers::RegisteredSolvers.begin();
        p != exahype::solvers::RegisteredSolvers.end(); p++) { // @todo replace by parloops?
      if (fineGridVerticesEnumerator.getLevel()>=(*p)->getMinimumTreeDepth()+1) {
        exahype::records::ADERDGCellDescription& cellDescriptionParent =
            fineGridCell.getADERDGCellDescription(solverNumber);
        // If coarse grid cell description requested refinement
        if (cellDescriptionParent.getVirtualRefinementNecessary()) {
          assertion(cellDescriptionParent.
                    getType()==exahype::Cell::RealCell
                    ||
                    cellDescriptionParent.
                    getType()==exahype::Cell::VirtualShell);
          assertion(cellDescriptionParent.getParent());

          fineGridCell.addNewCellDescription(
              solverNumber,
              exahype::Cell::VirtualShell,
              fineGridVerticesEnumerator.getLevel(),
              coarseGridCell.getADERDGCellDescriptionsIndex(),
              fineGridPositionOfCell,
              fineGridVerticesEnumerator.getCellSize(),
              fineGridVerticesEnumerator.getCellCenter());

          cellDescriptionParent.setVirtualRefinementNecessary(false);
        }
      }
      solverNumber++;
    }

  } else {
    bool refineFineGridCell=false;

    int  solverNumber=0;
    for (std::vector<exahype::solvers::Solver*>::const_iterator p =
        exahype::solvers::RegisteredSolvers.begin();
        p != exahype::solvers::RegisteredSolvers.end(); p++) { // @todo replace by parloops?

      assertion(static_cast<unsigned int>(solverNumber) <
                ADERDGCellDescriptionHeap::getInstance().
                getData(fineGridCell.getADERDGCellDescriptionsIndex()).size());

      exahype::records::ADERDGCellDescription& cellDescription =
          fineGridCell.getADERDGCellDescription(solverNumber);

      // check if virtual shell has real cell neighbours
      if (cellDescription.getType()==exahype::Cell::VirtualShell) {
        cellDescription.setHasNeighboursOfTypeCell(
            hasNeighboursOfType(
                solverNumber,exahype::Cell::RealCell,
                neighbourCellDescriptionIndices)
        );
      }

      // only refine non-parent real cells and virtual shells virtually
      // if they have real shells as neighbour
      if (!cellDescription.getParent()
          &&
          (cellDescription.getType()==exahype::Cell::RealCell
          ||
          cellDescription.getType()==exahype::Cell::VirtualShell)
          &&
          hasNeighboursOfType(
            solverNumber,exahype::Cell::RealShell,
            neighbourCellDescriptionIndices)) {
        cellDescription.setVirtualRefinementNecessary(true);
        cellDescription.setParent(true);

        refineFineGridCell = true;
      }
    }
    // Note that fineGridVertices->refine() refines all adjacent cell
    // not only the targeted fineGridCell (Peano Cookbook)
    if (refineFineGridCell) {
      dfor2(k) // loop over 2^d vertices and set the refinement flag
        if (fineGridVertices[kScalar].getRefinementControl()==
            Vertex::Records::Unrefined) {
          fineGridVertices->refine();
        }
      enddforx
    }
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

  bool exahype::mappings::VirtualRefinement::hasNeighboursOfType(
      const int solverNumber,
      exahype::Cell::CellDescriptionType cellType,
      const tarch::la::Vector<THREE_POWER_D, int>& neighbourCellDescriptionIndices) {
    // left,right,front,back,(front,back)
#if DIMENSIONS == 2
    constexpr int neighbourPositions[4] = {3,5,1,7};
#else
    constexpr int neighbourPositions[6] = {12,14,10,16,4,22};
#endif

    for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
      const int neighbourCellDescriptionIndex =
          neighbourCellDescriptionIndices[neighbourPositions[i]];

      if (DataHeap::getInstance().isValidIndex(neighbourCellDescriptionIndex)) {
        exahype::records::ADERDGCellDescription& neighbourCellDescription =
            ADERDGCellDescriptionHeap::getInstance().getData(
                neighbourCellDescriptionIndex)[solverNumber];

        assertion(neighbourCellDescription.getSolverNumber()==solverNumber);

        if (neighbourCellDescription.getType()==cellType) {
#if defined(Debug) || defined(Asserts)
          if (neighbourCellDescription.getType()==exahype::Cell::RealShell) {
            assertion(neighbourCellDescription.getParent());
          }
#endif
          return true;
        }
      }
    }
    return false;
  }

  void exahype::mappings::VirtualRefinement::leaveCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
    // do nothing
  }

  void exahype::mappings::VirtualRefinement::beginIteration(
      exahype::State& solverState) {
    // do nothing
  }

  void exahype::mappings::VirtualRefinement::endIteration(exahype::State& solverState) {
    // do nothing
  }

  void exahype::mappings::VirtualRefinement::descend(
      exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell) {
    // do nothing
  }

  void exahype::mappings::VirtualRefinement::ascend(
      exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell) {
    // do nothing
  }
