#include "ExplicitEulerForHeatEquation/adapters/CreateGridAndPlot.h"


peano::CommunicationSpecification   myproject::adapters::CreateGridAndPlot::communicationSpecification() {
  return peano::CommunicationSpecification::getMinimalSpecification()
   & myproject::mappings::CreateGrid::communicationSpecification()

   & myproject::adapters::CreateGridAndPlot2VTKPlotCellValue_0::communicationSpecification()
   & myproject::adapters::CreateGridAndPlot2VTKPlotVertexValue_1::communicationSpecification()

;
}


peano::MappingSpecification   myproject::adapters::CreateGridAndPlot::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & myproject::mappings::CreateGrid::touchVertexLastTimeSpecification()

   & myproject::adapters::CreateGridAndPlot2VTKPlotCellValue_0::touchVertexLastTimeSpecification()
   & myproject::adapters::CreateGridAndPlot2VTKPlotVertexValue_1::touchVertexLastTimeSpecification()

;
}


peano::MappingSpecification   myproject::adapters::CreateGridAndPlot::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification::getMinimalSpecification()
   & myproject::mappings::CreateGrid::touchVertexFirstTimeSpecification()

   & myproject::adapters::CreateGridAndPlot2VTKPlotCellValue_0::touchVertexFirstTimeSpecification()
   & myproject::adapters::CreateGridAndPlot2VTKPlotVertexValue_1::touchVertexFirstTimeSpecification()

;
}


peano::MappingSpecification   myproject::adapters::CreateGridAndPlot::enterCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & myproject::mappings::CreateGrid::enterCellSpecification()

   & myproject::adapters::CreateGridAndPlot2VTKPlotCellValue_0::enterCellSpecification()
   & myproject::adapters::CreateGridAndPlot2VTKPlotVertexValue_1::enterCellSpecification()

;
}


peano::MappingSpecification   myproject::adapters::CreateGridAndPlot::leaveCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & myproject::mappings::CreateGrid::leaveCellSpecification()

   & myproject::adapters::CreateGridAndPlot2VTKPlotCellValue_0::leaveCellSpecification()
   & myproject::adapters::CreateGridAndPlot2VTKPlotVertexValue_1::leaveCellSpecification()

;
}


peano::MappingSpecification   myproject::adapters::CreateGridAndPlot::ascendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & myproject::mappings::CreateGrid::ascendSpecification()

   & myproject::adapters::CreateGridAndPlot2VTKPlotCellValue_0::ascendSpecification()
   & myproject::adapters::CreateGridAndPlot2VTKPlotVertexValue_1::ascendSpecification()

;
}


peano::MappingSpecification   myproject::adapters::CreateGridAndPlot::descendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & myproject::mappings::CreateGrid::descendSpecification()

   & myproject::adapters::CreateGridAndPlot2VTKPlotCellValue_0::descendSpecification()
   & myproject::adapters::CreateGridAndPlot2VTKPlotVertexValue_1::descendSpecification()

;
}


myproject::adapters::CreateGridAndPlot::CreateGridAndPlot() {
}


myproject::adapters::CreateGridAndPlot::~CreateGridAndPlot() {
}


#if defined(SharedMemoryParallelisation)
myproject::adapters::CreateGridAndPlot::CreateGridAndPlot(const CreateGridAndPlot&  masterThread):
  _map2CreateGrid(masterThread._map2CreateGrid) 

 ,
  _map2CreateGridAndPlot2VTKPlotCellValue_0(masterThread._map2CreateGridAndPlot2VTKPlotCellValue_0) , 
  _map2CreateGridAndPlot2VTKPlotVertexValue_1(masterThread._map2CreateGridAndPlot2VTKPlotVertexValue_1) 

{
}


void myproject::adapters::CreateGridAndPlot::mergeWithWorkerThread(const CreateGridAndPlot& workerThread) {

  _map2CreateGrid.mergeWithWorkerThread(workerThread._map2CreateGrid);

  _map2CreateGridAndPlot2VTKPlotCellValue_0.mergeWithWorkerThread(workerThread._map2CreateGridAndPlot2VTKPlotCellValue_0);
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.mergeWithWorkerThread(workerThread._map2CreateGridAndPlot2VTKPlotVertexValue_1);

}
#endif


void myproject::adapters::CreateGridAndPlot::createHangingVertex(
      myproject::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      myproject::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      myproject::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  _map2CreateGridAndPlot2VTKPlotCellValue_0.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2CreateGrid.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void myproject::adapters::CreateGridAndPlot::destroyHangingVertex(
      const myproject::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2CreateGrid.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2CreateGridAndPlot2VTKPlotCellValue_0.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void myproject::adapters::CreateGridAndPlot::createInnerVertex(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2CreateGridAndPlot2VTKPlotCellValue_0.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2CreateGrid.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void myproject::adapters::CreateGridAndPlot::createBoundaryVertex(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2CreateGridAndPlot2VTKPlotCellValue_0.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2CreateGrid.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void myproject::adapters::CreateGridAndPlot::destroyVertex(
      const myproject::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {

  _map2CreateGrid.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

  _map2CreateGridAndPlot2VTKPlotCellValue_0.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void myproject::adapters::CreateGridAndPlot::createCell(
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2CreateGridAndPlot2VTKPlotCellValue_0.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );


  _map2CreateGrid.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void myproject::adapters::CreateGridAndPlot::destroyCell(
      const myproject::Cell&           fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {

  _map2CreateGrid.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

  _map2CreateGridAndPlot2VTKPlotCellValue_0.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}

#ifdef Parallel
void myproject::adapters::CreateGridAndPlot::mergeWithNeighbour(
  myproject::Vertex&  vertex,
  const myproject::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
   _map2CreateGridAndPlot2VTKPlotCellValue_0.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2CreateGridAndPlot2VTKPlotVertexValue_1.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );


   _map2CreateGrid.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );

}

void myproject::adapters::CreateGridAndPlot::prepareSendToNeighbour(
  myproject::Vertex&  vertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2CreateGrid.prepareSendToNeighbour( vertex, toRank, x, h, level );

   _map2CreateGridAndPlot2VTKPlotCellValue_0.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2CreateGridAndPlot2VTKPlotVertexValue_1.prepareSendToNeighbour( vertex, toRank, x, h, level );

}

void myproject::adapters::CreateGridAndPlot::prepareCopyToRemoteNode(
  myproject::Vertex&  localVertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2CreateGrid.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );

   _map2CreateGridAndPlot2VTKPlotCellValue_0.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2CreateGridAndPlot2VTKPlotVertexValue_1.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );

}

void myproject::adapters::CreateGridAndPlot::prepareCopyToRemoteNode(
  myproject::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
   _map2CreateGrid.prepareCopyToRemoteNode( localCell, toRank, x, h, level );

   _map2CreateGridAndPlot2VTKPlotCellValue_0.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2CreateGridAndPlot2VTKPlotVertexValue_1.prepareCopyToRemoteNode( localCell, toRank, x, h, level );

}

void myproject::adapters::CreateGridAndPlot::mergeWithRemoteDataDueToForkOrJoin(
  myproject::Vertex&  localVertex,
  const myproject::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2CreateGridAndPlot2VTKPlotCellValue_0.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2CreateGridAndPlot2VTKPlotVertexValue_1.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );


   _map2CreateGrid.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );

}

void myproject::adapters::CreateGridAndPlot::mergeWithRemoteDataDueToForkOrJoin(
  myproject::Cell&  localCell,
  const myproject::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2CreateGridAndPlot2VTKPlotCellValue_0.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2CreateGridAndPlot2VTKPlotVertexValue_1.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );


   _map2CreateGrid.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );

}

bool myproject::adapters::CreateGridAndPlot::prepareSendToWorker(
  myproject::Cell&                 fineGridCell,
  myproject::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  myproject::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  myproject::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  bool result = false;
   result |= _map2CreateGrid.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );

   result |= _map2CreateGridAndPlot2VTKPlotCellValue_0.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2CreateGridAndPlot2VTKPlotVertexValue_1.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );

  return result;
}

void myproject::adapters::CreateGridAndPlot::prepareSendToMaster(
  myproject::Cell&                       localCell,
  myproject::Vertex *                    vertices,
  const peano::grid::VertexEnumerator&       verticesEnumerator, 
  const myproject::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
  const myproject::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {

   _map2CreateGrid.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

   _map2CreateGridAndPlot2VTKPlotCellValue_0.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2CreateGridAndPlot2VTKPlotVertexValue_1.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}

void myproject::adapters::CreateGridAndPlot::mergeWithMaster(
  const myproject::Cell&           workerGridCell,
  myproject::Vertex * const        workerGridVertices,
  const peano::grid::VertexEnumerator& workerEnumerator,
  myproject::Cell&                 fineGridCell,
  myproject::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  myproject::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  myproject::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker,
    const myproject::State&          workerState,
  myproject::State&                masterState
) {
   _map2CreateGridAndPlot2VTKPlotCellValue_0.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2CreateGridAndPlot2VTKPlotVertexValue_1.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );


   _map2CreateGrid.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );

}

void myproject::adapters::CreateGridAndPlot::receiveDataFromMaster(
      myproject::Cell&                        receivedCell, 
      myproject::Vertex *                     receivedVertices,
      const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
      myproject::Vertex * const               receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
      myproject::Cell&                        receivedCoarseGridCell,
      myproject::Vertex * const               workersCoarseGridVertices,
      const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
      myproject::Cell&                        workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
) {
   _map2CreateGridAndPlot2VTKPlotCellValue_0.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2CreateGridAndPlot2VTKPlotVertexValue_1.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );


   _map2CreateGrid.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );

}


void myproject::adapters::CreateGridAndPlot::mergeWithWorker(
  myproject::Cell&           localCell, 
  const myproject::Cell&     receivedMasterCell,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level
) {
   _map2CreateGridAndPlot2VTKPlotCellValue_0.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2CreateGridAndPlot2VTKPlotVertexValue_1.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );


   _map2CreateGrid.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );

}

void myproject::adapters::CreateGridAndPlot::mergeWithWorker(
  myproject::Vertex&        localVertex,
  const myproject::Vertex&  receivedMasterVertex,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {

   _map2CreateGrid.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );

   _map2CreateGridAndPlot2VTKPlotCellValue_0.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2CreateGridAndPlot2VTKPlotVertexValue_1.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );

}
#endif

void myproject::adapters::CreateGridAndPlot::touchVertexFirstTime(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2CreateGridAndPlot2VTKPlotCellValue_0.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2CreateGrid.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void myproject::adapters::CreateGridAndPlot::touchVertexLastTime(
      myproject::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {

  _map2CreateGrid.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

  _map2CreateGridAndPlot2VTKPlotCellValue_0.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void myproject::adapters::CreateGridAndPlot::enterCell(
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2CreateGridAndPlot2VTKPlotCellValue_0.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );


  _map2CreateGrid.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void myproject::adapters::CreateGridAndPlot::leaveCell(
      myproject::Cell&           fineGridCell,
      myproject::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {

  _map2CreateGrid.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

  _map2CreateGridAndPlot2VTKPlotCellValue_0.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void myproject::adapters::CreateGridAndPlot::beginIteration(
  myproject::State&  solverState
) {
  _map2CreateGridAndPlot2VTKPlotCellValue_0.beginIteration( solverState );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.beginIteration( solverState );


  _map2CreateGrid.beginIteration( solverState );

}


void myproject::adapters::CreateGridAndPlot::endIteration(
  myproject::State&  solverState
) {

  _map2CreateGrid.endIteration( solverState );

  _map2CreateGridAndPlot2VTKPlotCellValue_0.endIteration( solverState );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.endIteration( solverState );

}




void myproject::adapters::CreateGridAndPlot::descend(
  myproject::Cell * const          fineGridCells,
  myproject::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  myproject::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  myproject::Cell&                 coarseGridCell
) {
  _map2CreateGridAndPlot2VTKPlotCellValue_0.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );


  _map2CreateGrid.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}


void myproject::adapters::CreateGridAndPlot::ascend(
  myproject::Cell * const    fineGridCells,
  myproject::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  myproject::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  myproject::Cell&           coarseGridCell
) {

  _map2CreateGrid.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

  _map2CreateGridAndPlot2VTKPlotCellValue_0.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2CreateGridAndPlot2VTKPlotVertexValue_1.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}
