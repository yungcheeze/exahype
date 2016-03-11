#include "EulerFlow/mappings/FaceDataExchange.h"

#include "EulerFlow/dg/Constants.h"

#include "EulerFlow/geometry/Mapping.h"
#include "EulerFlow/quad/GaussLegendre.h"

#include "EulerFlow/multiscalelinkedcell/HangingVertexBookkeeper.h"
#include "EulerFlow/VertexOperations.h"

#include "EulerFlow/problem/Problem.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::FaceDataExchange::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::
          SendDataAndStateAfterLastTouchVertexLastTime,
      false);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::FaceDataExchange::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::FaceDataExchange::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::FaceDataExchange::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::OnlyLeaves,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::FaceDataExchange::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::FaceDataExchange::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::FaceDataExchange::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::FaceDataExchange::_log(
    "exahype::mappings::FaceDataExchange");

exahype::mappings::FaceDataExchange::FaceDataExchange() {
  // do nothing
}

exahype::mappings::FaceDataExchange::~FaceDataExchange() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::FaceDataExchange::FaceDataExchange(
    const FaceDataExchange& masterThread) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::mergeWithWorkerThread(
    const FaceDataExchange& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::FaceDataExchange::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::FaceDataExchange::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::FaceDataExchange::prepareSendToWorker(
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

void exahype::mappings::FaceDataExchange::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::mergeWithMaster(
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

void exahype::mappings::FaceDataExchange::receiveDataFromMaster(
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

void exahype::mappings::FaceDataExchange::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::FaceDataExchange::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::setBoundaryGhostValues(
    records::CellDescription& dataSelf, const int patchIndexGhostSelf,
    const int patchIndexRealSelf, const int face, const int numberOfDofOnFace) {
  for (int dof = 0; dof < numberOfDofOnFace; dof++) {
    // copy extrapolated predictor
    DataHeap::getInstance()
        .getData(dataSelf.getExtrapolatedPredictor(
            patchIndexGhostSelf))[face * numberOfDofOnFace + dof]
        ._persistentRecords._u =
        DataHeap::getInstance()
            .getData(dataSelf.getExtrapolatedPredictor(
                patchIndexRealSelf))[face * numberOfDofOnFace + dof]
            ._persistentRecords._u;
    // copy extrapolated fluctuation/normal flux
    DataHeap::getInstance()
        .getData(dataSelf.getFluctuation(
            patchIndexGhostSelf))[face * numberOfDofOnFace + dof]
        ._persistentRecords._u =
        DataHeap::getInstance()
            .getData(dataSelf.getFluctuation(
                patchIndexRealSelf))[face * numberOfDofOnFace + dof]
            ._persistentRecords._u;
  }
}

void exahype::mappings::FaceDataExchange::copyGhostValues(
    records::CellDescription& dataSelf,
    const records::CellDescription& dataNeighbour,
    const int patchIndexGhostSelf, const int patchIndexNeighbour,
    const int face, const int numberOfDofOnFace) {
  for (int dof = 0; dof < numberOfDofOnFace; dof++) {
    double* QSelfGhost = &(DataHeap::getInstance()
                               .getData(dataSelf.getExtrapolatedPredictor(
                                   patchIndexGhostSelf))[0]
                               ._persistentRecords._u);
    double* QRealNeighbour =
        &(DataHeap::getInstance()
              .getData(dataNeighbour.getExtrapolatedPredictor(
                  patchIndexNeighbour))[0]
              ._persistentRecords._u);

    // copy extrapolated predictor
    DataHeap::getInstance()
        .getData(dataSelf.getExtrapolatedPredictor(
            patchIndexGhostSelf))[face * numberOfDofOnFace + dof]
        ._persistentRecords._u =
        DataHeap::getInstance()
            .getData(dataNeighbour.getExtrapolatedPredictor(
                patchIndexNeighbour))[face * numberOfDofOnFace + dof]
            ._persistentRecords._u;
    // copy extrapolated fluctuation/normal flux
    DataHeap::getInstance()
        .getData(dataSelf.getFluctuation(
            patchIndexGhostSelf))[face * numberOfDofOnFace + dof]
        ._persistentRecords._u =
        DataHeap::getInstance()
            .getData(dataNeighbour.getFluctuation(
                patchIndexNeighbour))[face * numberOfDofOnFace + dof]
            ._persistentRecords._u;
  }
}

// Begin of code for ADERDG method
void exahype::mappings::FaceDataExchange::initialiseGhostLayerOfPatches(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  int basisSize = EXAHYPE_ORDER + 1;
  int numberOfDofOnFace =
      EXAHYPE_NVARS * tarch::la::aPowI(DIMENSIONS - 1, basisSize);

  records::CellDescription& cellDescriptionSelf =
      CellDescriptionHeap::getInstance().getData(
          fineGridCell.getCellDescriptionsIndex())[0];

  // Read in neighbour information
  const tarch::la::Vector<THREE_POWER_D, int> cellDescriptionsOfAllNeighbours =
      multiscalelinkedcell::getIndicesAroundCell(
          exahype::VertexOperations::readCellDescriptionsIndex(
              fineGridVerticesEnumerator, fineGridVertices));

  int indexGhostSelf = 0;
  int indexRealSelf = 0;
  int indexRealNeighbor = 0;

  for (int i = 1; i < EXAHYPE_PATCH_SIZE_X + 1; i++) {
    ///////////////////////////////////////
    // FRONT
    ///////////////////////////////////////
    indexGhostSelf = i + (EXAHYPE_PATCH_SIZE_X + 2) * 0;
    indexRealSelf = i + (EXAHYPE_PATCH_SIZE_X + 2) * 1;
    indexRealNeighbor = i + (EXAHYPE_PATCH_SIZE_X + 2) * EXAHYPE_PATCH_SIZE_Y;

    if (cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_FRONT] ==
        multiscalelinkedcell::HangingVertexBookkeeper::
            DomainBoundaryAdjacencyIndex) {
      setBoundaryGhostValues(cellDescriptionSelf, indexGhostSelf, indexRealSelf,
                             EXAHYPE_FACE_BACK, numberOfDofOnFace);
    } else if (cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_FRONT] >
               multiscalelinkedcell::HangingVertexBookkeeper::
                   InvalidAdjacencyIndex) {
      records::CellDescription& cellDescriptionNeighbourFront =
          CellDescriptionHeap::getInstance().getData(
              cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_FRONT])[0];

      copyGhostValues(cellDescriptionSelf, cellDescriptionNeighbourFront,
                      indexGhostSelf, indexRealNeighbor, EXAHYPE_FACE_BACK,
                      numberOfDofOnFace);
    }
    ///////////////////////////////////////
    // BACK
    ///////////////////////////////////////
    indexGhostSelf =
        i + (EXAHYPE_PATCH_SIZE_X + 2) * (EXAHYPE_PATCH_SIZE_Y + 1);  // ghost
    indexRealSelf =
        i + (EXAHYPE_PATCH_SIZE_X + 2) * (EXAHYPE_PATCH_SIZE_Y);  // non-ghost
    indexRealNeighbor = i + (EXAHYPE_PATCH_SIZE_X + 2) * 1;       // non-ghost

    if (cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_BACK] ==
        multiscalelinkedcell::HangingVertexBookkeeper::
            DomainBoundaryAdjacencyIndex) {
      setBoundaryGhostValues(cellDescriptionSelf, indexGhostSelf, indexRealSelf,
                             EXAHYPE_FACE_FRONT, numberOfDofOnFace);
    } else if (cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_BACK] >
               multiscalelinkedcell::HangingVertexBookkeeper::
                   InvalidAdjacencyIndex) {
      records::CellDescription& cellDescriptionNeighbourBack =
          CellDescriptionHeap::getInstance().getData(
              cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_BACK])[0];

      copyGhostValues(cellDescriptionSelf, cellDescriptionNeighbourBack,
                      indexGhostSelf, indexRealNeighbor, EXAHYPE_FACE_FRONT,
                      numberOfDofOnFace);
    }

    for (int j = 1; j < EXAHYPE_PATCH_SIZE_Y + 1; j++) {
      ///////////////////////////////////////
      // LEFT
      ///////////////////////////////////////
      indexGhostSelf = 0 + (EXAHYPE_PATCH_SIZE_X + 2) * j;  // ghost
      indexRealSelf = 1 + (EXAHYPE_PATCH_SIZE_X + 2) * j;   // non-ghost
      indexRealNeighbor =
          EXAHYPE_PATCH_SIZE_X + (EXAHYPE_PATCH_SIZE_X + 2) * j;  // non-ghost

      if (cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_LEFT] ==
          multiscalelinkedcell::HangingVertexBookkeeper::
              DomainBoundaryAdjacencyIndex) {
        setBoundaryGhostValues(cellDescriptionSelf, indexGhostSelf,
                               indexRealSelf, EXAHYPE_FACE_RIGHT,
                               numberOfDofOnFace);
      } else if (cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_LEFT] >
                 multiscalelinkedcell::HangingVertexBookkeeper::
                     InvalidAdjacencyIndex) {
        records::CellDescription& cellDescriptionNeighbourLeft =
            CellDescriptionHeap::getInstance().getData(
                cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_LEFT])[0];
        copyGhostValues(cellDescriptionSelf, cellDescriptionNeighbourLeft,
                        indexGhostSelf, indexRealNeighbor, EXAHYPE_FACE_RIGHT,
                        numberOfDofOnFace);
      }

      ///////////////////////////////////////
      // RIGHT
      ///////////////////////////////////////
      indexGhostSelf =
          (EXAHYPE_PATCH_SIZE_X + 1) + (EXAHYPE_PATCH_SIZE_X + 2) * j;  // ghost
      indexRealSelf =
          EXAHYPE_PATCH_SIZE_X + (EXAHYPE_PATCH_SIZE_X + 2) * j;  // non-ghost
      indexRealNeighbor = 1 + (EXAHYPE_PATCH_SIZE_X + 2) * j;     // non-ghost

      if (cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_RIGHT] ==
          multiscalelinkedcell::HangingVertexBookkeeper::
              DomainBoundaryAdjacencyIndex) {
        setBoundaryGhostValues(cellDescriptionSelf, indexGhostSelf,
                               indexRealSelf, EXAHYPE_FACE_LEFT,
                               numberOfDofOnFace);
      } else if (cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_RIGHT] >
                 multiscalelinkedcell::HangingVertexBookkeeper::
                     InvalidAdjacencyIndex) {
        records::CellDescription& cellDescriptionNeighbourRight =
            CellDescriptionHeap::getInstance().getData(
                cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_RIGHT])[0];
        copyGhostValues(cellDescriptionSelf, cellDescriptionNeighbourRight,
                        indexGhostSelf, indexRealNeighbor, EXAHYPE_FACE_LEFT,
                        numberOfDofOnFace);
      }
    }
  }
}

void exahype::mappings::FaceDataExchange::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (!fineGridCell.isRefined()) {
    initialiseGhostLayerOfPatches(fineGridCell, fineGridVertices,
                                  fineGridVerticesEnumerator);
  }

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::FaceDataExchange::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::beginIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::endIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::FaceDataExchange::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
// End of code for ADERDG method
