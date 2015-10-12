#include "ExaHyPE/adapters/CreateGrid.h"


peano::CommunicationSpecification   ExaHyPE::adapters::CreateGrid::communicationSpecification() {
  return peano::CommunicationSpecification::getMinimalSpecification()
   & ExaHyPE::mappings::CreateGrid::communicationSpecification()

   & ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::communicationSpecification()

;
}


peano::MappingSpecification   ExaHyPE::adapters::CreateGrid::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & ExaHyPE::mappings::CreateGrid::touchVertexLastTimeSpecification()

   & ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::touchVertexLastTimeSpecification()

;
}


peano::MappingSpecification   ExaHyPE::adapters::CreateGrid::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification::getMinimalSpecification()
   & ExaHyPE::mappings::CreateGrid::touchVertexFirstTimeSpecification()

   & ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::touchVertexFirstTimeSpecification()

;
}


peano::MappingSpecification   ExaHyPE::adapters::CreateGrid::enterCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & ExaHyPE::mappings::CreateGrid::enterCellSpecification()

   & ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::enterCellSpecification()

;
}


peano::MappingSpecification   ExaHyPE::adapters::CreateGrid::leaveCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & ExaHyPE::mappings::CreateGrid::leaveCellSpecification()

   & ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::leaveCellSpecification()

;
}


peano::MappingSpecification   ExaHyPE::adapters::CreateGrid::ascendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & ExaHyPE::mappings::CreateGrid::ascendSpecification()

   & ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::ascendSpecification()

;
}


peano::MappingSpecification   ExaHyPE::adapters::CreateGrid::descendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & ExaHyPE::mappings::CreateGrid::descendSpecification()

   & ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::descendSpecification()

;
}


ExaHyPE::adapters::CreateGrid::CreateGrid() {
}


ExaHyPE::adapters::CreateGrid::~CreateGrid() {
}


#if defined(SharedMemoryParallelisation)
ExaHyPE::adapters::CreateGrid::CreateGrid(const CreateGrid&  masterThread):
  _map2CreateGrid(masterThread._map2CreateGrid) 

 ,
  _map2CreateGrid2VTKGridVisualiser_0(masterThread._map2CreateGrid2VTKGridVisualiser_0) 

{
}


void ExaHyPE::adapters::CreateGrid::mergeWithWorkerThread(const CreateGrid& workerThread) {

  _map2CreateGrid.mergeWithWorkerThread(workerThread._map2CreateGrid);

  _map2CreateGrid2VTKGridVisualiser_0.mergeWithWorkerThread(workerThread._map2CreateGrid2VTKGridVisualiser_0);

}
#endif


void ExaHyPE::adapters::CreateGrid::createHangingVertex(
      ExaHyPE::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      ExaHyPE::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      ExaHyPE::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  _map2CreateGrid2VTKGridVisualiser_0.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2CreateGrid.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void ExaHyPE::adapters::CreateGrid::destroyHangingVertex(
      const ExaHyPE::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2CreateGrid.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2CreateGrid2VTKGridVisualiser_0.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void ExaHyPE::adapters::CreateGrid::createInnerVertex(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2CreateGrid2VTKGridVisualiser_0.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2CreateGrid.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void ExaHyPE::adapters::CreateGrid::createBoundaryVertex(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2CreateGrid2VTKGridVisualiser_0.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2CreateGrid.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void ExaHyPE::adapters::CreateGrid::destroyVertex(
      const ExaHyPE::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {

  _map2CreateGrid.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

  _map2CreateGrid2VTKGridVisualiser_0.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void ExaHyPE::adapters::CreateGrid::createCell(
      ExaHyPE::Cell&                 fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2CreateGrid2VTKGridVisualiser_0.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );


  _map2CreateGrid.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void ExaHyPE::adapters::CreateGrid::destroyCell(
      const ExaHyPE::Cell&           fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {

  _map2CreateGrid.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

  _map2CreateGrid2VTKGridVisualiser_0.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}

#ifdef Parallel
void ExaHyPE::adapters::CreateGrid::mergeWithNeighbour(
  ExaHyPE::Vertex&  vertex,
  const ExaHyPE::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
   _map2CreateGrid2VTKGridVisualiser_0.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );


   _map2CreateGrid.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );

}

void ExaHyPE::adapters::CreateGrid::prepareSendToNeighbour(
  ExaHyPE::Vertex&  vertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2CreateGrid.prepareSendToNeighbour( vertex, toRank, x, h, level );

   _map2CreateGrid2VTKGridVisualiser_0.prepareSendToNeighbour( vertex, toRank, x, h, level );

}

void ExaHyPE::adapters::CreateGrid::prepareCopyToRemoteNode(
  ExaHyPE::Vertex&  localVertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2CreateGrid.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );

   _map2CreateGrid2VTKGridVisualiser_0.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );

}

void ExaHyPE::adapters::CreateGrid::prepareCopyToRemoteNode(
  ExaHyPE::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
   _map2CreateGrid.prepareCopyToRemoteNode( localCell, toRank, x, h, level );

   _map2CreateGrid2VTKGridVisualiser_0.prepareCopyToRemoteNode( localCell, toRank, x, h, level );

}

void ExaHyPE::adapters::CreateGrid::mergeWithRemoteDataDueToForkOrJoin(
  ExaHyPE::Vertex&  localVertex,
  const ExaHyPE::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2CreateGrid2VTKGridVisualiser_0.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );


   _map2CreateGrid.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );

}

void ExaHyPE::adapters::CreateGrid::mergeWithRemoteDataDueToForkOrJoin(
  ExaHyPE::Cell&  localCell,
  const ExaHyPE::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2CreateGrid2VTKGridVisualiser_0.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );


   _map2CreateGrid.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );

}

bool ExaHyPE::adapters::CreateGrid::prepareSendToWorker(
  ExaHyPE::Cell&                 fineGridCell,
  ExaHyPE::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  ExaHyPE::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  ExaHyPE::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  bool result = false;
   result |= _map2CreateGrid.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );

   result |= _map2CreateGrid2VTKGridVisualiser_0.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );

  return result;
}

void ExaHyPE::adapters::CreateGrid::prepareSendToMaster(
  ExaHyPE::Cell&                       localCell,
  ExaHyPE::Vertex *                    vertices,
  const peano::grid::VertexEnumerator&       verticesEnumerator, 
  const ExaHyPE::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
  const ExaHyPE::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {

   _map2CreateGrid.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

   _map2CreateGrid2VTKGridVisualiser_0.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}

void ExaHyPE::adapters::CreateGrid::mergeWithMaster(
  const ExaHyPE::Cell&           workerGridCell,
  ExaHyPE::Vertex * const        workerGridVertices,
  const peano::grid::VertexEnumerator& workerEnumerator,
  ExaHyPE::Cell&                 fineGridCell,
  ExaHyPE::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  ExaHyPE::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  ExaHyPE::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker,
    const ExaHyPE::State&          workerState,
  ExaHyPE::State&                masterState
) {
   _map2CreateGrid2VTKGridVisualiser_0.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );


   _map2CreateGrid.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );

}

void ExaHyPE::adapters::CreateGrid::receiveDataFromMaster(
      ExaHyPE::Cell&                        receivedCell, 
      ExaHyPE::Vertex *                     receivedVertices,
      const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
      ExaHyPE::Vertex * const               receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
      ExaHyPE::Cell&                        receivedCoarseGridCell,
      ExaHyPE::Vertex * const               workersCoarseGridVertices,
      const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
      ExaHyPE::Cell&                        workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
) {
   _map2CreateGrid2VTKGridVisualiser_0.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );


   _map2CreateGrid.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );

}


void ExaHyPE::adapters::CreateGrid::mergeWithWorker(
  ExaHyPE::Cell&           localCell, 
  const ExaHyPE::Cell&     receivedMasterCell,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level
) {
   _map2CreateGrid2VTKGridVisualiser_0.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );


   _map2CreateGrid.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );

}

void ExaHyPE::adapters::CreateGrid::mergeWithWorker(
  ExaHyPE::Vertex&        localVertex,
  const ExaHyPE::Vertex&  receivedMasterVertex,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {

   _map2CreateGrid.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );

   _map2CreateGrid2VTKGridVisualiser_0.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );

}
#endif

void ExaHyPE::adapters::CreateGrid::touchVertexFirstTime(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2CreateGrid2VTKGridVisualiser_0.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2CreateGrid.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void ExaHyPE::adapters::CreateGrid::touchVertexLastTime(
      ExaHyPE::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {

  _map2CreateGrid.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

  _map2CreateGrid2VTKGridVisualiser_0.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void ExaHyPE::adapters::CreateGrid::enterCell(
      ExaHyPE::Cell&                 fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2CreateGrid2VTKGridVisualiser_0.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );


  _map2CreateGrid.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void ExaHyPE::adapters::CreateGrid::leaveCell(
      ExaHyPE::Cell&           fineGridCell,
      ExaHyPE::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {

  _map2CreateGrid.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

  _map2CreateGrid2VTKGridVisualiser_0.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void ExaHyPE::adapters::CreateGrid::beginIteration(
  ExaHyPE::State&  solverState
) {
  _map2CreateGrid2VTKGridVisualiser_0.beginIteration( solverState );


  _map2CreateGrid.beginIteration( solverState );

}


void ExaHyPE::adapters::CreateGrid::endIteration(
  ExaHyPE::State&  solverState
) {

  _map2CreateGrid.endIteration( solverState );

  _map2CreateGrid2VTKGridVisualiser_0.endIteration( solverState );

}




void ExaHyPE::adapters::CreateGrid::descend(
  ExaHyPE::Cell * const          fineGridCells,
  ExaHyPE::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  ExaHyPE::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  ExaHyPE::Cell&                 coarseGridCell
) {
  _map2CreateGrid2VTKGridVisualiser_0.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );


  _map2CreateGrid.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}


void ExaHyPE::adapters::CreateGrid::ascend(
  ExaHyPE::Cell * const    fineGridCells,
  ExaHyPE::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  ExaHyPE::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  ExaHyPE::Cell&           coarseGridCell
) {

  _map2CreateGrid.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

  _map2CreateGrid2VTKGridVisualiser_0.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}
