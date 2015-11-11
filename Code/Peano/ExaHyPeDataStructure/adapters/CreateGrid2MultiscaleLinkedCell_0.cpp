#include "ExaHyPeDataStructure/adapters/CreateGrid2MultiscaleLinkedCell_0.h"

#include <sstream>

#include "peano/utils/Loop.h"
#include "peano/grid/CellFlags.h"


#include "ExaHyPeDataStructure/multiscalelinkedcell/HangingVertexBookkeeper.h"
#include "ExaHyPeDataStructure/VertexOperations.h"


peano::CommunicationSpecification   exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::communicationSpecification() {
  return peano::CommunicationSpecification::getMinimalSpecification();
}


peano::MappingSpecification   exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification   exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification   exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification   exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification   exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification   exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::CreateGrid2MultiscaleLinkedCell_0() {
}


exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::~CreateGrid2MultiscaleLinkedCell_0() {
}


#if defined(SharedMemoryParallelisation)
exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::CreateGrid2MultiscaleLinkedCell_0(const CreateGrid2MultiscaleLinkedCell_0&  masterThread) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::mergeWithWorkerThread(const CreateGrid2MultiscaleLinkedCell_0& workerThread) {
}
#endif


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::createHangingVertex(
  exahype::Vertex&     fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
  exahype::Vertex * const   coarseGridVertices,
  const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
  exahype::Cell&       coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  const int level = coarseGridVerticesEnumerator.getLevel()+1;
  
  fineGridVertex.getPatchIndex() = 
    multiscalelinkedcell::HangingVertexBookkeeper::getInstance().createHangingVertex(
      fineGridX,level,
      fineGridPositionOfVertex,
      VertexOperations::readPatchIndex(coarseGridVerticesEnumerator,coarseGridVertices)
    );
}



void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::destroyHangingVertex(
  const exahype::Vertex&   fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::createInnerVertex(
  exahype::Vertex&               fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  fineGridVertex.getPatchIndex() = multiscalelinkedcell::HangingVertexBookkeeper::getInstance().createVertexLinkMapForNewVertex();
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::createBoundaryVertex(
  exahype::Vertex&               fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  fineGridVertex.getPatchIndex() = multiscalelinkedcell::HangingVertexBookkeeper::getInstance().createVertexLinkMapForBoundaryVertex();
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::destroyVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::createCell(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::destroyCell(
  const exahype::Cell&           fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  multiscalelinkedcell::HangingVertexBookkeeper::getInstance().destroyCell(fineGridCell.getPatchIndex());
}


#ifdef Parallel
void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::mergeWithNeighbour(
  exahype::Vertex&  vertex,
  const exahype::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
  multiscalelinkedcell::HangingVertexBookkeeper::getInstance().updateCellIndicesInMergeWithNeighbour(vertex);
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::prepareSendToNeighbour(
      exahype::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::prepareCopyToRemoteNode(
      exahype::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::prepareCopyToRemoteNode(
      exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Vertex&  localVertex,
  const exahype::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Cell&  localCell,
  const exahype::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


bool exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::prepareSendToWorker(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  return false;
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::prepareSendToMaster(
      exahype::Cell&                       localCell,
      exahype::Vertex *                    vertices,
      const peano::grid::VertexEnumerator&       verticesEnumerator, 
      const exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      const exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::mergeWithMaster(
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
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::receiveDataFromMaster(
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
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::mergeWithWorker(
      exahype::Cell&           localCell, 
      const exahype::Cell&     receivedMasterCell,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::mergeWithWorker(
      exahype::Vertex&        localVertex,
      const exahype::Vertex&  receivedMasterVertex,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}
#endif


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::enterCell(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  dfor2(k)
    if (fineGridVertices[fineGridVerticesEnumerator(k)].isHangingNode()) {
      multiscalelinkedcell::HangingVertexBookkeeper::getInstance().getAdjacencyEntriesOfVertex( 
        fineGridVerticesEnumerator.getVertexPosition(k),
        fineGridVerticesEnumerator.getLevel()
      )(TWO_POWER_D-kScalar-1) = fineGridCell.getPatchIndex();
    }
    else {
      fineGridVertices[fineGridVerticesEnumerator(k)].getPatchIndex()(TWO_POWER_D-kScalar-1) = fineGridCell.getPatchIndex();
    }
  enddforx
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::leaveCell(
      exahype::Cell&           fineGridCell,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::beginIteration(
  exahype::State&  solverState
) {
  multiscalelinkedcell::HangingVertexBookkeeper::getInstance().beginIteration();
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::endIteration(
  exahype::State&  solverState
) {
  multiscalelinkedcell::HangingVertexBookkeeper::getInstance().endIteration();
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::descend(
  exahype::Cell * const          fineGridCells,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell
) {
}


void exahype::adapters::CreateGrid2MultiscaleLinkedCell_0::ascend(
  exahype::Cell * const    fineGridCells,
  exahype::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell
) {
}
