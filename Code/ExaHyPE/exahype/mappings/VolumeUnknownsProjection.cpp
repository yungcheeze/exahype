#include "exahype/mappings/VolumeUnknownsProjection.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "kernels/KernelCalls.h"
#include "exahype/solvers/Solver.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/VertexOperations.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::VolumeUnknownsProjection::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::
      SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::
      SendDataAndStateAfterLastTouchVertexLastTime,
      false);
}

peano::MappingSpecification
exahype::mappings::VolumeUnknownsProjection::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::VolumeUnknownsProjection::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::VolumeUnknownsProjection::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::VolumeUnknownsProjection::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::VolumeUnknownsProjection::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
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
exahype::mappings::VolumeUnknownsProjection::VolumeUnknownsProjection(const VolumeUnknownsProjection& masterThread) {
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

void exahype::mappings::VolumeUnknownsProjection::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::VolumeUnknownsProjection::mergeWithRemoteDataDueToForkOrJoin(
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
  if (!ADERDGCellDescriptionHeap::getInstance().
      isValidIndex(fineGridCell.getADERDGCellDescriptionsIndex())) {
    assertion(ADERDGCellDescriptionHeap::getInstance().
              isValidIndex(coarseGridCell.getADERDGCellDescriptionsIndex()));

    const int numberOfADERDGCellDescriptions = static_cast<int>(
        ADERDGCellDescriptionHeap::getInstance()
        .getData(fineGridCell.getADERDGCellDescriptionsIndex())
        .size());
    // please use a different UserDefined per mapping/event
    const peano::datatraversal::autotuning::MethodTrace methodTrace =
        peano::datatraversal::autotuning::UserDefined1;
    const int grainSize =
        peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
            numberOfADERDGCellDescriptions, methodTrace);
    pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
      records::ADERDGCellDescription& p =
        fineGridCell.getADERDGCellDescription(i);

      // if we have at least one parent
      if (ADERDGCellDescriptionHeap::getInstance().
          isValidIndex(p.getParentIndex())) {
        switch (p.getType()) {
        case exahype::records::ADERDGCellDescription::Cell:
          if () {
            exahype::Cell::SubcellPosition subcellPosition =
                fineGridCell.getSubcellPositionOfVirtualShell(p);
            prolongateVolumeData(
                p,
                subcellPosition.parentIndex,
                subcellPosition.subcellIndex);
            }
          break;
        case exahype::records::ADERDGCellDescription::Cell:
        case p.getType() == exahype::records::ADERDGCellDescription::Shell:
        restrictVolumeData(
            p,
            p.getParentIndex(),
            p.getFineGridPositionOfCell());
        break;
        }
      }

    endpfor peano::datatraversal::autotuning::Oracle::getInstance()
    .parallelSectionHasTerminated(methodTrace);
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}


void exahype::mappings::VolumeUnknownsProjection::prolongateVolumeData(
                const exahype::records::ADERDGCellDescription& cellDescription,
                const int parentIndex,
                const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) const {
  exahype::records::ADERDGCellDescription& cellDescriptionParent =
      ADERDGCellDescriptionHeap::getInstance().
        getData(parentIndex)[cellDescription.getSolverNumber()];

  assertion(cellDescriptionParent.getSolverNumber()==
          cellDescription.getSolverNumber());
  assertion(cellDescriptionParent.getParent());
  assertion(cellDescriptionParent.getType()==exahype::records::ADERDGCellDescription::Cell
            ||
            cellDescriptionParent.getType()==exahype::records::ADERDGCellDescription::VirtualShell);
#if defined(Debug) || defined(Asserts)
              if (cellDescriptionParent.getType()==exahype::records::ADERDGCellDescription::VirtualShell) {
                assertion(cellDescriptionParent.getHasNeighboursOfTypeCell());
              }
#endif

  const int levelFine   = cellDescription.getLevel();
  const int levelCoarse = cellDescriptionParent.getLevel();
  assertion(levelCoarse < levelFine);
  const int levelDelta  = levelFine - levelCoarse;



      double* lQhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getExtrapolatedPredictor()).
          data();
      double* lQhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getExtrapolatedPredictor()).
          data();

      exahype::solvers::Solver* solver = exahype::solvers::
          RegisteredSolvers[cellDescription.getSolverNumber()];
      solver->volumeUnknownsProlongation(
          lQhbndFine,
          lQhbndCoarse,
          levelCoarse,levelFine,
          getSubVolumeIndex(subcellIndex));
}

void exahype::mappings::VolumeUnknownsProjection::restrictVolumeData(
                const exahype::records::ADERDGCellDescription& cellDescription,
                const int parentIndex,
                const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) const {
  exahype::records::ADERDGCellDescription& cellDescriptionParent =
      ADERDGCellDescriptionHeap::getInstance().
        getData(parentIndex)[cellDescription.getSolverNumber()];

  assertion(cellDescriptionParent.getSolverNumber()==
          cellDescription.getSolverNumber());
  assertion(cellDescriptionParent.getParent());
  assertion(cellDescriptionParent.getType()==exahype::records::ADERDGCellDescription::Shell);

  const int levelFine   = cellDescription.getLevel();
  const int levelCoarse = cellDescriptionParent.getLevel();
  assertion(levelCoarse < levelFine);

  for (int d=0; d < DIMENSIONS; d++) {
    // Check if cell is at "left" or "right" d Volume of parent
    if (subcellIndex[d]==0) {
      const int VolumeIndex = 2*d;

      double* lQhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getExtrapolatedPredictor()).
          data()[VolumeIndex];
      double* lQhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getExtrapolatedPredictor()).
          data()[VolumeIndex];
      double* lFhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getFluctuation()).
          data()[VolumeIndex];
      double* lFhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getFluctuation()).
          data()[VolumeIndex];

      exahype::solvers::Solver* solver = exahype::solvers::
          RegisteredSolvers[cellDescription.getSolverNumber()];
      solver->VolumeUnknownsRestriction(
          lQhbndCoarse,lFhbndCoarse,
          lQhbndFine,lFhbndFine,
          levelCoarse,levelFine,
          getSubVolumeIndex(subcellIndex,d));

    } else if (subcellIndex[d]==2) {
      const int VolumeIndex = 2*d+1;

      double* lQhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getExtrapolatedPredictor()).
          data()[VolumeIndex];
      double* lQhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getExtrapolatedPredictor()).
          data()[VolumeIndex];
      double* lFhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getFluctuation()).
          data()[VolumeIndex];
      double* lFhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getFluctuation()).
          data()[VolumeIndex];

      exahype::solvers::Solver* solver = exahype::solvers::
          RegisteredSolvers[cellDescription.getSolverNumber()];
      solver->VolumeUnknownsRestriction(
          lQhbndCoarse,lFhbndCoarse,
          lQhbndFine,lFhbndFine,
          levelCoarse,levelFine,
          getSubVolumeIndex(subcellIndex,d));
    }
  }
}

tarch::la::Vector<DIMENSIONS-1,int> exahype::mappings::VolumeUnknownsProjection::getSubVolumeIndex(
    const tarch::la::Vector<DIMENSIONS,int>& subcellIndex,
    const int d) const {
  tarch::la::Vector<DIMENSIONS-1,int> subVolumeIndex;

  int i=0;
  for (int j=0; j<DIMENSIONS; j++) {
    if (j!=d) {
      subVolumeIndex[i] = subcellIndex[j];
      i++;
    }
  }

  return subVolumeIndex;
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

void exahype::mappings::VolumeUnknownsProjection::endIteration(exahype::State& solverState) {
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
  // do nothing
}
