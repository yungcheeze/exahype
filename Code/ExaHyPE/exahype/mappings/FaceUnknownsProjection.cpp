#include "exahype/mappings/FaceUnknownsProjection.h"

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
exahype::mappings::FaceUnknownsProjection::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::FaceUnknownsProjection::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::FaceUnknownsProjection::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::FaceUnknownsProjection::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::FaceUnknownsProjection::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::FaceUnknownsProjection::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::FaceUnknownsProjection::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::FaceUnknownsProjection::_log(
    "exahype::mappings::FaceUnknownsProjection");

exahype::mappings::FaceUnknownsProjection::FaceUnknownsProjection() {
  // do nothing
}

exahype::mappings::FaceUnknownsProjection::~FaceUnknownsProjection() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::FaceUnknownsProjection::FaceUnknownsProjection(const FaceUnknownsProjection& masterThread) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::mergeWithWorkerThread(
    const FaceUnknownsProjection& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::FaceUnknownsProjection::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::FaceUnknownsProjection::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::FaceUnknownsProjection::prepareSendToWorker(
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

void exahype::mappings::FaceUnknownsProjection::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::mergeWithMaster(
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

void exahype::mappings::FaceUnknownsProjection::receiveDataFromMaster(
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

void exahype::mappings::FaceUnknownsProjection::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::FaceUnknownsProjection::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::enterCell(
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
        for (std::vector<exahype::records::ADERDGCellDescription>::
            iterator pParent = ADERDGCellDescriptionHeap::getInstance().
            getData(p.getParentIndex()).begin();
            pParent != ADERDGCellDescriptionHeap::getInstance().
                getData(p.getParentIndex()).end();
            ++pParent) {
          if (p.getSolverNumber()==pParent->getSolverNumber()) {
            switch (p.getType()) {
              case exahype::records::ADERDGCellDescription::Descendant:
                exahype::Cell::SubcellPosition subcellPosition =
                    fineGridCell.computeSubcellPositionOfDescendant(p);

                prolongateFaceData(
                    p,
                    subcellPosition.parentIndex,
                    subcellPosition.subcellIndex);
                break;
              case exahype::records::ADERDGCellDescription::Cell:
              case exahype::records::ADERDGCellDescription::Ancestor:
                exahype::Cell::SubcellPosition subcellPosition =
                    fineGridCell.computeSubcellPositionOfCellOrAncestor(p);

                restrictFaceData(
                    p,
                    subcellPosition.parentIndex,
                    subcellPosition.subcellIndex);
                break;
            }
          }
        }
      }
    endpfor peano::datatraversal::autotuning::Oracle::getInstance()
    .parallelSectionHasTerminated(methodTrace);
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}


void exahype::mappings::FaceUnknownsProjection::prolongateFaceData(
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

  for (int d=0; d < DIMENSIONS; d++) {
    // Check if cell is at "left" or "right" d face of parent
    if (subcellIndex[d]==0) {
      const int faceIndex = 2*d;

      double* lQhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getExtrapolatedPredictor()).
          data()[faceIndex];
      double* lQhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getExtrapolatedPredictor()).
          data()[faceIndex];
      double* lFhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getFluctuation()).
          data()[faceIndex];
      double* lFhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getFluctuation()).
          data()[faceIndex];

      exahype::solvers::Solver* solver = exahype::solvers::
          RegisteredSolvers[cellDescription.getSolverNumber()];
      solver->faceUnknownsProlongation(
          lQhbndFine,lFhbndFine,
          lQhbndCoarse,lFhbndCoarse,
          levelCoarse,levelFine,
          getSubfaceIndex(subcellIndex,d));

    } else if (subcellIndex[d]==tarch::la::aPowI(levelDelta,3)-1) {
      const int faceIndex = 2*d+1;

      double* lQhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getExtrapolatedPredictor()).
          data()[faceIndex];
      double* lQhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getExtrapolatedPredictor()).
          data()[faceIndex];
      double* lFhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getFluctuation()).
          data()[faceIndex];
      double* lFhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getFluctuation()).
          data()[faceIndex];

      exahype::solvers::Solver* solver = exahype::solvers::
          RegisteredSolvers[cellDescription.getSolverNumber()];
      solver->faceUnknownsProlongation(
          lQhbndFine,lFhbndFine,
          lQhbndCoarse,lFhbndCoarse,
          levelCoarse,levelFine,
          getSubfaceIndex(subcellIndex,d));
    }
  }
}

void exahype::mappings::FaceUnknownsProjection::restrictFaceData(
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
    // Check if cell is at "left" or "right" d face of parent
    if (subcellIndex[d]==0) {
      const int faceIndex = 2*d;

      double* lQhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getExtrapolatedPredictor()).
          data()[faceIndex];
      double* lQhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getExtrapolatedPredictor()).
          data()[faceIndex];
      double* lFhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getFluctuation()).
          data()[faceIndex];
      double* lFhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getFluctuation()).
          data()[faceIndex];

      exahype::solvers::Solver* solver = exahype::solvers::
          RegisteredSolvers[cellDescription.getSolverNumber()];
      solver->faceUnknownsRestriction(
          lQhbndCoarse,lFhbndCoarse,
          lQhbndFine,lFhbndFine,
          levelCoarse,levelFine,
          getSubfaceIndex(subcellIndex,d));

    } else if (subcellIndex[d]==2) {
      const int faceIndex = 2*d+1;

      double* lQhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getExtrapolatedPredictor()).
          data()[faceIndex];
      double* lQhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getExtrapolatedPredictor()).
          data()[faceIndex];
      double* lFhbndFine = &DataHeap::getInstance().
          getData(cellDescription.getFluctuation()).
          data()[faceIndex];
      double* lFhbndCoarse = &DataHeap::getInstance().
          getData(cellDescriptionParent.getFluctuation()).
          data()[faceIndex];

      exahype::solvers::Solver* solver = exahype::solvers::
          RegisteredSolvers[cellDescription.getSolverNumber()];
      solver->faceUnknownsRestriction(
          lQhbndCoarse,lFhbndCoarse,
          lQhbndFine,lFhbndFine,
          levelCoarse,levelFine,
          getSubfaceIndex(subcellIndex,d));
    }
  }
}

tarch::la::Vector<DIMENSIONS-1,int> exahype::mappings::FaceUnknownsProjection::getSubfaceIndex(
    const tarch::la::Vector<DIMENSIONS,int>& subcellIndex,
    const int d) const {
  tarch::la::Vector<DIMENSIONS-1,int> subfaceIndex;

  int i=0;
  for (int j=0; j<DIMENSIONS; j++) {
    if (j!=d) {
      subfaceIndex[i] = subcellIndex[j];
      i++;
    }
  }

  return subfaceIndex;
}

void exahype::mappings::FaceUnknownsProjection::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::beginIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::endIteration(exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
