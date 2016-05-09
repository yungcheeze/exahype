#include "exahype/mappings/BoundaryConditions.h"

#include "peano/utils/Globals.h"

#include "tarch/multicore/Lock.h"
#include "tarch/multicore/Loop.h"
#include "peano/datatraversal/autotuning/Oracle.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/solvers/Solver.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::BoundaryConditions::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::BoundaryConditions::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::BoundaryConditions::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::BoundaryConditions::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::BoundaryConditions::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::BoundaryConditions::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::BoundaryConditions::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::BoundaryConditions::_log(
    "exahype::mappings::BoundaryConditions");

exahype::mappings::BoundaryConditions::BoundaryConditions() {
  // do nothing
}

exahype::mappings::BoundaryConditions::~BoundaryConditions() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::BoundaryConditions::BoundaryConditions(
    const BoundaryConditions& masterThread)
: _localState(masterThread._localState) {}

void exahype::mappings::BoundaryConditions::mergeWithWorkerThread(
    const BoundaryConditions& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::BoundaryConditions::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::BoundaryConditions::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::BoundaryConditions::prepareSendToWorker(
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

void exahype::mappings::BoundaryConditions::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::mergeWithMaster(
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

void exahype::mappings::BoundaryConditions::receiveDataFromMaster(
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

void exahype::mappings::BoundaryConditions::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::BoundaryConditions::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("touchVertextFirstTime(...)", fineGridVertex,
      fineGridX, fineGridH,
      coarseGridVerticesEnumerator.toString(),
      coarseGridCell, fineGridPositionOfVertex);

  tarch::la::Vector<TWO_POWER_D, int>& adjacentADERDGCellDescriptionsIndices =
      fineGridVertex.getADERDGCellDescriptionsIndex();
  /*
   * Right cell-left cell   pair indices: 0,1; 2,3;   4,5; 6;7
   * Front cell-back cell   pair indices: 0,2; 1,3;   4,6; 5;7
   * Top   cell-bottom cell pair indices: 0,4; 1,5;   2,6; 3;7
   *
   * Note that from the viewpoint of a cell, the face
   * has always the "opposite" index, i.e., we solve a Riemann
   * problem on the left face of the right cell (which
   * is the right face of the left cell).
   */
  constexpr int cellIndicesLeft[4] = {0, 2, 4, 6};
  constexpr int cellIndicesRight[4] = {1, 3, 5, 7};
  constexpr int cellIndicesFront[4] = {0, 1, 4, 5};
  constexpr int cellIndicesBack[4] = {2, 3, 6, 7};
#if DIMENSIONS == 3
  constexpr int cellIndicesBottom[4] = {0, 1, 2, 3};
  constexpr int cellIndicesTop[4] = {4, 5, 6, 7};
#endif

  // Left/right face
  for (int i = 0; i < TWO_POWER_D_DIVIDED_BY_TWO; i++) {
    applyBoundaryConditions(adjacentADERDGCellDescriptionsIndices,
        cellIndicesLeft[i], cellIndicesRight[i],
        EXAHYPE_FACE_RIGHT, EXAHYPE_FACE_LEFT, 0);

    applyBoundaryConditions(adjacentADERDGCellDescriptionsIndices,
        cellIndicesFront[i], cellIndicesBack[i],
        EXAHYPE_FACE_BACK, EXAHYPE_FACE_FRONT, 1);

#if DIMENSIONS == 3
    applyBoundaryConditions(adjacentADERDGCellDescriptionsIndices,
        cellIndicesBottom[i], cellIndicesTop[i],
        EXAHYPE_FACE_TOP, EXAHYPE_FACE_BOTTOM, 2);
#endif
  }
  logTraceOutWith1Argument("touchVertextFirstTime(...)", fineGridVertex);
}

void exahype::mappings::BoundaryConditions::applyBoundaryConditions(
    tarch::la::Vector<TWO_POWER_D, int>& adjacentADERDGCellDescriptionsIndices,
    const int cellIndexL, const int cellIndexR, const int faceIndexL,
    const int faceIndexR, const int normalNonZero) {
  int cellIndex =
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex;
  int faceIndex = -1;

  bool noBoundaryFace = true;
  bool bothAreBoundaryFaces = false;
  if (adjacentADERDGCellDescriptionsIndices[cellIndexL] ==
      multiscalelinkedcell::HangingVertexBookkeeper::
      DomainBoundaryAdjacencyIndex) {
    noBoundaryFace = false;
    bothAreBoundaryFaces = true;
    cellIndex = cellIndexR;
    faceIndex = faceIndexR;
  }

  if (adjacentADERDGCellDescriptionsIndices[cellIndexR] ==
      multiscalelinkedcell::HangingVertexBookkeeper::
      DomainBoundaryAdjacencyIndex) {
    noBoundaryFace = false;
    bothAreBoundaryFaces = bothAreBoundaryFaces & true;
    cellIndex = cellIndexL;
    faceIndex = faceIndexL;
  }

  // Only continue if this is a boundary face.
  if (!noBoundaryFace && !bothAreBoundaryFaces) {
    if (ADERDGCellDescriptionHeap::getInstance().
        isValidIndex(adjacentADERDGCellDescriptionsIndices[cellIndex])) {

      assertion1(adjacentADERDGCellDescriptionsIndices[cellIndexL] >
      multiscalelinkedcell::HangingVertexBookkeeper::
      InvalidAdjacencyIndex
      ||
      adjacentADERDGCellDescriptionsIndices[cellIndexR] >
      multiscalelinkedcell::HangingVertexBookkeeper::
      InvalidAdjacencyIndex,
      adjacentADERDGCellDescriptionsIndices.toString());

      const auto numberOfADERDGCellDescriptions =
          static_cast<int>(ADERDGCellDescriptionHeap::getInstance()
      .getData(adjacentADERDGCellDescriptionsIndices[cellIndex])
      .size());

      const peano::datatraversal::autotuning::MethodTrace methodTrace =
          peano::datatraversal::autotuning::UserDefined0;
      const int grainSize =
          peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
              numberOfADERDGCellDescriptions, methodTrace);
      pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
      records::ADERDGCellDescription& p =
          ADERDGCellDescriptionHeap::getInstance().getData(
              adjacentADERDGCellDescriptionsIndices[cellIndex])[i];

      exahype::solvers::Solver* solver =
          exahype::solvers::RegisteredSolvers[p.getSolverNumber()];
      bool riemannSolveNotPerformed = false;

      switch(p.getType()) {
      case exahype::records::ADERDGCellDescription::Cell:
        switch(p.getRefinementEvent()) {
        case exahype::records::ADERDGCellDescription::None:
          // Lock the critical multithreading area.
          // Note that two boundary vertices can operate on the same face at the same
          // time.
          riemannSolveNotPerformed = false;
          riemannSolveNotPerformed = !p.getRiemannSolvePerformed(faceIndex);
          if (riemannSolveNotPerformed) {
            p.setRiemannSolvePerformed(faceIndex, true);
          }

          // Apply the boundary conditions.
          // todo Copy and periodic boundary conditions should be configured by
          // additional mappings after
          // the initial grid refinement and the initialisation of the cell
          // descriptions.
          // The configuration should only involve an edit of the index maps generated
          // by the multiscalelinkedcell toolbox. In case of consequent mesh-refinement,
          // these indices should be propagated down to the finer cells.
          // The resulting Riemann problems are then simply solved
          // by exahype::mappings::RiemannSolver.
          if (riemannSolveNotPerformed) {
            // @todo 03/02/16:Dominic Etienne Charrier
            // Change to solver->getUnknownsPerFace()
            const int numberOfFaceDof =
                solver->getUnknownsPerFace();

            double* lQhbnd =
                DataHeap::getInstance().getData(
                    p.getExtrapolatedPredictor()).data() +
                    (faceIndex * numberOfFaceDof);
            double* lFhbnd =
                DataHeap::getInstance().getData(
                    p.getFluctuation()).data() +
                    (faceIndex * numberOfFaceDof);

            // @todo
            // timestepping::synchroniseTimeStepping(solve,*p);
            assertionEquals(lQhbnd[0],lQhbnd[0]); // assert no nan
            assertionEquals(lFhbnd[0],lFhbnd[0]); // assert no nan

            // todo Dominic Charrier2503
            // Do not solve a Riemann problem here:
            // Invoke user defined boundary condition function
            // At the moment, we simply copy the cell solution to the boundary.
            solver->riemannSolver(
                lFhbnd, lFhbnd, lQhbnd, lQhbnd,
                p.getCorrectorTimeStepSize(),  // solve.getCorrectorTimeStepSize(),//_localState.getPreviousMinTimeStepSize(),
                normalNonZero);

            assertionEquals(lFhbnd[0],lFhbnd[0]); // assert no nan
          }
          break;
        default:
          break;
        }
        break;
        default:
          break;
      }
      endpfor peano::datatraversal::autotuning::Oracle::getInstance()
      .parallelSectionHasTerminated(methodTrace);
    }
  }
}

void exahype::mappings::BoundaryConditions::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  _localState = solverState;

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::BoundaryConditions::endIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
