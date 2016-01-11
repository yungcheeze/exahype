#include "EulerFlow/mappings/VolumeIntegral.h"

#include "EulerFlow/Constants.h"

#include "EulerFlow/quad/GaussLegendre.h"

#include "EulerFlow/geometry/Mapping.h"

#include "EulerFlow/problem/Problem.h"

#include "EulerFlow/dg/Constants.h"
#include "EulerFlow/dg/ADERDG.h"
#include "EulerFlow/dg/DGMatrices.h"

#include "stdlib.h"

#include "string.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification   exahype::mappings::VolumeIntegral::communicationSpecification() {
  return peano::CommunicationSpecification(peano::CommunicationSpecification::SendDataAndStateBeforeFirstTouchVertexFirstTime,peano::CommunicationSpecification::SendDataAndStateAfterLastTouchVertexLastTime,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::VolumeIntegral::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::VolumeIntegral::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::VolumeIntegral::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::OnlyLeaves,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::VolumeIntegral::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::VolumeIntegral::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::VolumeIntegral::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


tarch::logging::Log                exahype::mappings::VolumeIntegral::_log( "exahype::mappings::VolumeIntegral" ); 


exahype::mappings::VolumeIntegral::VolumeIntegral() {
  // do nothing
}


exahype::mappings::VolumeIntegral::~VolumeIntegral() {
  // do nothing
}


#if defined(SharedMemoryParallelisation)
exahype::mappings::VolumeIntegral::VolumeIntegral(const VolumeIntegral&  masterThread) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::mergeWithWorkerThread(const VolumeIntegral& workerThread) {
  // do nothing
}
#endif


void exahype::mappings::VolumeIntegral::createHangingVertex(
    exahype::Vertex&     fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
    exahype::Vertex * const   coarseGridVertices,
    const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
    exahype::Cell&       coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::destroyHangingVertex(
    const exahype::Vertex&   fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::createInnerVertex(
    exahype::Vertex&               fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::createBoundaryVertex(
    exahype::Vertex&               fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::destroyVertex(
    const exahype::Vertex&   fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::createCell(
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::destroyCell(
    const exahype::Cell&           fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::VolumeIntegral::mergeWithNeighbour(
    exahype::Vertex&  vertex,
    const exahype::Vertex&  neighbour,
    int                                           fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::VolumeIntegral::prepareSendToNeighbour(
    exahype::Vertex&  vertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::VolumeIntegral::prepareCopyToRemoteNode(
    exahype::Vertex&  localVertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::VolumeIntegral::prepareCopyToRemoteNode(
    exahype::Cell&  localCell,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::VolumeIntegral::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex&  localVertex,
    const exahype::Vertex&  masterOrWorkerVertex,
    int                                       fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&  x,
    const tarch::la::Vector<DIMENSIONS,double>&  h,
    int                                       level
) {
  // do nothing
}

void exahype::mappings::VolumeIntegral::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell&  localCell,
    const exahype::Cell&  masterOrWorkerCell,
    int                                       fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                       level
) {
  // do nothing
}

bool exahype::mappings::VolumeIntegral::prepareSendToWorker(
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
    int                                                                  worker
) {
  // do nothing
  return true;
}

void exahype::mappings::VolumeIntegral::prepareSendToMaster(
    exahype::Cell&                       localCell,
    exahype::Vertex *                    vertices,
    const peano::grid::VertexEnumerator&       verticesEnumerator,
    const exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
    const exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::mergeWithMaster(
    const exahype::Cell&           workerGridCell,
    exahype::Vertex * const        workerGridVertices,
    const peano::grid::VertexEnumerator& workerEnumerator,
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
    int                                                                  worker,
    const exahype::State&          workerState,
    exahype::State&                masterState
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::receiveDataFromMaster(
    exahype::Cell&                        receivedCell,
    exahype::Vertex *                     receivedVertices,
    const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
    exahype::Vertex * const               receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
    exahype::Cell&                        receivedCoarseGridCell,
    exahype::Vertex * const               workersCoarseGridVertices,
    const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
    exahype::Cell&                        workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::mergeWithWorker(
    exahype::Cell&           localCell,
    const exahype::Cell&     receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                          level
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::mergeWithWorker(
    exahype::Vertex&        localVertex,
    const exahype::Vertex&  receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}
#endif

void exahype::mappings::VolumeIntegral::touchVertexFirstTime(
    exahype::Vertex&               fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  // do nothing
}

extern const double exahype::quad::gaussLegendreNodes[EXAHYPE_ORDER+1];
extern const double exahype::quad::gaussLegendreWeights[EXAHYPE_ORDER+1];

void exahype::mappings::VolumeIntegral::touchVertexLastTime(
    exahype::Vertex&         fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  // do nothing
}

void exahype::mappings::VolumeIntegral::enterCell(
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "enterCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );

  if (!fineGridCell.isRefined()) {
    records::CellDescription& cellDescription =
        CellDescriptionHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[0];

    const tarch::la::Vector<DIMENSIONS,double> center = fineGridVerticesEnumerator.getCellCenter();  // the center of the cell
    const double dx = fineGridVerticesEnumerator.getCellSize()(0);
    const double dy = fineGridVerticesEnumerator.getCellSize()(1);

    const double dxPatch[2] = { dx/ (double) EXAHYPE_PATCH_SIZE_X,
                                dy/ (double) EXAHYPE_PATCH_SIZE_Y };

    for (int i=1; i<EXAHYPE_PATCH_SIZE_X+1; i++) { // loop over patches
      for (int j=1; j<EXAHYPE_PATCH_SIZE_Y+1; j++) {
        const int patchIndex = i + (EXAHYPE_PATCH_SIZE_X+2) * j;

        double * lFhi = &(DataHeap::getInstance().getData(cellDescription.getVolumeFlux(patchIndex))[0]._persistentRecords._u);
        double * lduh = &(DataHeap::getInstance().getData(cellDescription.getUpdate(patchIndex))    [0]._persistentRecords._u);

        dg::volumeIntegral<DIMENSIONS>(
            lduh,
            lFhi,
            dxPatch);
      }
    }
  }

  logTraceOutWith1Argument( "enterCell(...)", fineGridCell );
}


void exahype::mappings::VolumeIntegral::leaveCell(
    exahype::Cell&           fineGridCell,
    exahype::Vertex * const  fineGridVertices,
    const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::beginIteration(
    exahype::State&  solverState
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::endIteration(
    exahype::State&  solverState
) {
  // do nothing
}



void exahype::mappings::VolumeIntegral::descend(
    exahype::Cell * const          fineGridCells,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell
) {
  // do nothing
}


void exahype::mappings::VolumeIntegral::ascend(
    exahype::Cell * const    fineGridCells,
    exahype::Vertex * const  fineGridVertices,
    const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell
) {
  // do nothing
}
