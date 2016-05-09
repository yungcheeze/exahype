#include "exahype/adapters/PlotAugmentedAMRGrid.h"


peano::CommunicationSpecification   exahype::adapters::PlotAugmentedAMRGrid::communicationSpecification() {
  return peano::CommunicationSpecification::getMinimalSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::communicationSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2MultiscaleLinkedCell_1::communicationSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::PlotAugmentedAMRGrid::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::touchVertexLastTimeSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2MultiscaleLinkedCell_1::touchVertexLastTimeSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::PlotAugmentedAMRGrid::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::touchVertexFirstTimeSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2MultiscaleLinkedCell_1::touchVertexFirstTimeSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::PlotAugmentedAMRGrid::enterCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::enterCellSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2MultiscaleLinkedCell_1::enterCellSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::PlotAugmentedAMRGrid::leaveCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::leaveCellSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2MultiscaleLinkedCell_1::leaveCellSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::PlotAugmentedAMRGrid::ascendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::ascendSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2MultiscaleLinkedCell_1::ascendSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::PlotAugmentedAMRGrid::descendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::descendSpecification()
   & exahype::adapters::PlotAugmentedAMRGrid2MultiscaleLinkedCell_1::descendSpecification()

  ;
}


exahype::adapters::PlotAugmentedAMRGrid::PlotAugmentedAMRGrid() {
}


exahype::adapters::PlotAugmentedAMRGrid::~PlotAugmentedAMRGrid() {
}


#if defined(SharedMemoryParallelisation)
exahype::adapters::PlotAugmentedAMRGrid::PlotAugmentedAMRGrid(const PlotAugmentedAMRGrid&  masterThread):
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0(masterThread._map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0) , 
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1(masterThread._map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1) 

{
}


void exahype::adapters::PlotAugmentedAMRGrid::mergeWithWorkerThread(const PlotAugmentedAMRGrid& workerThread) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.mergeWithWorkerThread(workerThread._map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0);
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.mergeWithWorkerThread(workerThread._map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1);

}
#endif


void exahype::adapters::PlotAugmentedAMRGrid::createHangingVertex(
      exahype::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      exahype::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      exahype::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void exahype::adapters::PlotAugmentedAMRGrid::destroyHangingVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::PlotAugmentedAMRGrid::createInnerVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::PlotAugmentedAMRGrid::createBoundaryVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::PlotAugmentedAMRGrid::destroyVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::PlotAugmentedAMRGrid::createCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::PlotAugmentedAMRGrid::destroyCell(
      const exahype::Cell&           fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


#ifdef Parallel
void exahype::adapters::PlotAugmentedAMRGrid::mergeWithNeighbour(
  exahype::Vertex&  vertex,
  const exahype::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
   _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );

}


void exahype::adapters::PlotAugmentedAMRGrid::prepareSendToNeighbour(
  exahype::Vertex&  vertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.prepareSendToNeighbour( vertex, toRank, x, h, level );

}


void exahype::adapters::PlotAugmentedAMRGrid::prepareCopyToRemoteNode(
  exahype::Vertex&  localVertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );

}


void exahype::adapters::PlotAugmentedAMRGrid::prepareCopyToRemoteNode(
  exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
   _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.prepareCopyToRemoteNode( localCell, toRank, x, h, level );

}


void exahype::adapters::PlotAugmentedAMRGrid::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Vertex&  localVertex,
  const exahype::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );

}


void exahype::adapters::PlotAugmentedAMRGrid::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Cell&  localCell,
  const exahype::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );

}


bool exahype::adapters::PlotAugmentedAMRGrid::prepareSendToWorker(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  bool result = false;
   result |= _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );

  return result;
}


void exahype::adapters::PlotAugmentedAMRGrid::prepareSendToMaster(
  exahype::Cell&                       localCell,
  exahype::Vertex *                    vertices,
  const peano::grid::VertexEnumerator&       verticesEnumerator, 
  const exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
  const exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
   _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::PlotAugmentedAMRGrid::mergeWithMaster(
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
   _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );

}


void exahype::adapters::PlotAugmentedAMRGrid::receiveDataFromMaster(
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
   _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::PlotAugmentedAMRGrid::mergeWithWorker(
  exahype::Cell&           localCell, 
  const exahype::Cell&     receivedMasterCell,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level
) {
   _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );

}


void exahype::adapters::PlotAugmentedAMRGrid::mergeWithWorker(
  exahype::Vertex&        localVertex,
  const exahype::Vertex&  receivedMasterVertex,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );

}
#endif


void exahype::adapters::PlotAugmentedAMRGrid::touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::PlotAugmentedAMRGrid::touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::PlotAugmentedAMRGrid::enterCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::PlotAugmentedAMRGrid::leaveCell(
      exahype::Cell&           fineGridCell,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::PlotAugmentedAMRGrid::beginIteration(
  exahype::State&  solverState
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.beginIteration( solverState );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.beginIteration( solverState );

}


void exahype::adapters::PlotAugmentedAMRGrid::endIteration(
  exahype::State&  solverState
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.endIteration( solverState );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.endIteration( solverState );

}




void exahype::adapters::PlotAugmentedAMRGrid::descend(
  exahype::Cell * const          fineGridCells,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}


void exahype::adapters::PlotAugmentedAMRGrid::ascend(
  exahype::Cell * const    fineGridCells,
  exahype::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell
) {
  _map2PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2PlotAugmentedAMRGrid2MultiscaleLinkedCell_1.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}
