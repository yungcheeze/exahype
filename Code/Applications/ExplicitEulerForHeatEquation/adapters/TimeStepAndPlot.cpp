#include "ExplicitEulerForHeatEquation/adapters/TimeStepAndPlot.h"


peano::CommunicationSpecification   myproject::adapters::TimeStepAndPlot::communicationSpecification() {
  return peano::CommunicationSpecification::getMinimalSpecification()
   & myproject::mappings::TimeStep::communicationSpecification()
   & myproject::mappings::InterpolateHangingNodes::communicationSpecification()
   & myproject::mappings::Inject::communicationSpecification()

   & myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::communicationSpecification()

;
}


peano::MappingSpecification   myproject::adapters::TimeStepAndPlot::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & myproject::mappings::TimeStep::touchVertexLastTimeSpecification()
   & myproject::mappings::InterpolateHangingNodes::touchVertexLastTimeSpecification()
   & myproject::mappings::Inject::touchVertexLastTimeSpecification()

   & myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::touchVertexLastTimeSpecification()

;
}


peano::MappingSpecification   myproject::adapters::TimeStepAndPlot::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification::getMinimalSpecification()
   & myproject::mappings::TimeStep::touchVertexFirstTimeSpecification()
   & myproject::mappings::InterpolateHangingNodes::touchVertexFirstTimeSpecification()
   & myproject::mappings::Inject::touchVertexFirstTimeSpecification()

   & myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::touchVertexFirstTimeSpecification()

;
}


peano::MappingSpecification   myproject::adapters::TimeStepAndPlot::enterCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & myproject::mappings::TimeStep::enterCellSpecification()
   & myproject::mappings::InterpolateHangingNodes::enterCellSpecification()
   & myproject::mappings::Inject::enterCellSpecification()

   & myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::enterCellSpecification()

;
}


peano::MappingSpecification   myproject::adapters::TimeStepAndPlot::leaveCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & myproject::mappings::TimeStep::leaveCellSpecification()
   & myproject::mappings::InterpolateHangingNodes::leaveCellSpecification()
   & myproject::mappings::Inject::leaveCellSpecification()

   & myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::leaveCellSpecification()

;
}


peano::MappingSpecification   myproject::adapters::TimeStepAndPlot::ascendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & myproject::mappings::TimeStep::ascendSpecification()
   & myproject::mappings::InterpolateHangingNodes::ascendSpecification()
   & myproject::mappings::Inject::ascendSpecification()

   & myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::ascendSpecification()

;
}


peano::MappingSpecification   myproject::adapters::TimeStepAndPlot::descendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & myproject::mappings::TimeStep::descendSpecification()
   & myproject::mappings::InterpolateHangingNodes::descendSpecification()
   & myproject::mappings::Inject::descendSpecification()

   & myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::descendSpecification()

;
}


myproject::adapters::TimeStepAndPlot::TimeStepAndPlot() {
}


myproject::adapters::TimeStepAndPlot::~TimeStepAndPlot() {
}


#if defined(SharedMemoryParallelisation)
myproject::adapters::TimeStepAndPlot::TimeStepAndPlot(const TimeStepAndPlot&  masterThread):
  _map2TimeStep(masterThread._map2TimeStep) , 
  _map2InterpolateHangingNodes(masterThread._map2InterpolateHangingNodes) , 
  _map2Inject(masterThread._map2Inject) 

 ,
  _map2TimeStepAndPlot2VTKPlotVertexValue_0(masterThread._map2TimeStepAndPlot2VTKPlotVertexValue_0) 

{
}


void myproject::adapters::TimeStepAndPlot::mergeWithWorkerThread(const TimeStepAndPlot& workerThread) {

  _map2TimeStep.mergeWithWorkerThread(workerThread._map2TimeStep);
  _map2InterpolateHangingNodes.mergeWithWorkerThread(workerThread._map2InterpolateHangingNodes);
  _map2Inject.mergeWithWorkerThread(workerThread._map2Inject);

  _map2TimeStepAndPlot2VTKPlotVertexValue_0.mergeWithWorkerThread(workerThread._map2TimeStepAndPlot2VTKPlotVertexValue_0);

}
#endif


void myproject::adapters::TimeStepAndPlot::createHangingVertex(
      myproject::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      myproject::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      myproject::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  _map2TimeStepAndPlot2VTKPlotVertexValue_0.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2TimeStep.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2InterpolateHangingNodes.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Inject.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void myproject::adapters::TimeStepAndPlot::destroyHangingVertex(
      const myproject::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2TimeStep.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2InterpolateHangingNodes.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Inject.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2TimeStepAndPlot2VTKPlotVertexValue_0.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void myproject::adapters::TimeStepAndPlot::createInnerVertex(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2TimeStepAndPlot2VTKPlotVertexValue_0.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2TimeStep.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2InterpolateHangingNodes.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Inject.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void myproject::adapters::TimeStepAndPlot::createBoundaryVertex(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2TimeStepAndPlot2VTKPlotVertexValue_0.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2TimeStep.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2InterpolateHangingNodes.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Inject.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void myproject::adapters::TimeStepAndPlot::destroyVertex(
      const myproject::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {

  _map2TimeStep.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2InterpolateHangingNodes.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Inject.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

  _map2TimeStepAndPlot2VTKPlotVertexValue_0.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void myproject::adapters::TimeStepAndPlot::createCell(
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2TimeStepAndPlot2VTKPlotVertexValue_0.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );


  _map2TimeStep.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2InterpolateHangingNodes.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2Inject.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void myproject::adapters::TimeStepAndPlot::destroyCell(
      const myproject::Cell&           fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {

  _map2TimeStep.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2InterpolateHangingNodes.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2Inject.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

  _map2TimeStepAndPlot2VTKPlotVertexValue_0.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}

#ifdef Parallel
void myproject::adapters::TimeStepAndPlot::mergeWithNeighbour(
  myproject::Vertex&  vertex,
  const myproject::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
   _map2TimeStepAndPlot2VTKPlotVertexValue_0.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );


   _map2TimeStep.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2InterpolateHangingNodes.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2Inject.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );

}

void myproject::adapters::TimeStepAndPlot::prepareSendToNeighbour(
  myproject::Vertex&  vertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2TimeStep.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2InterpolateHangingNodes.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2Inject.prepareSendToNeighbour( vertex, toRank, x, h, level );

   _map2TimeStepAndPlot2VTKPlotVertexValue_0.prepareSendToNeighbour( vertex, toRank, x, h, level );

}

void myproject::adapters::TimeStepAndPlot::prepareCopyToRemoteNode(
  myproject::Vertex&  localVertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2TimeStep.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2InterpolateHangingNodes.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2Inject.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );

   _map2TimeStepAndPlot2VTKPlotVertexValue_0.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );

}

void myproject::adapters::TimeStepAndPlot::prepareCopyToRemoteNode(
  myproject::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
   _map2TimeStep.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2InterpolateHangingNodes.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2Inject.prepareCopyToRemoteNode( localCell, toRank, x, h, level );

   _map2TimeStepAndPlot2VTKPlotVertexValue_0.prepareCopyToRemoteNode( localCell, toRank, x, h, level );

}

void myproject::adapters::TimeStepAndPlot::mergeWithRemoteDataDueToForkOrJoin(
  myproject::Vertex&  localVertex,
  const myproject::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2TimeStepAndPlot2VTKPlotVertexValue_0.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );


   _map2TimeStep.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2InterpolateHangingNodes.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2Inject.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );

}

void myproject::adapters::TimeStepAndPlot::mergeWithRemoteDataDueToForkOrJoin(
  myproject::Cell&  localCell,
  const myproject::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2TimeStepAndPlot2VTKPlotVertexValue_0.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );


   _map2TimeStep.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2InterpolateHangingNodes.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2Inject.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );

}

bool myproject::adapters::TimeStepAndPlot::prepareSendToWorker(
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
   result |= _map2TimeStep.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2InterpolateHangingNodes.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2Inject.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );

   result |= _map2TimeStepAndPlot2VTKPlotVertexValue_0.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );

  return result;
}

void myproject::adapters::TimeStepAndPlot::prepareSendToMaster(
  myproject::Cell&                       localCell,
  myproject::Vertex *                    vertices,
  const peano::grid::VertexEnumerator&       verticesEnumerator, 
  const myproject::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
  const myproject::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {

   _map2TimeStep.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2InterpolateHangingNodes.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2Inject.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

   _map2TimeStepAndPlot2VTKPlotVertexValue_0.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}

void myproject::adapters::TimeStepAndPlot::mergeWithMaster(
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
   _map2TimeStepAndPlot2VTKPlotVertexValue_0.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );


   _map2TimeStep.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2InterpolateHangingNodes.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2Inject.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );

}

void myproject::adapters::TimeStepAndPlot::receiveDataFromMaster(
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
   _map2TimeStepAndPlot2VTKPlotVertexValue_0.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );


   _map2TimeStep.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2InterpolateHangingNodes.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2Inject.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );

}


void myproject::adapters::TimeStepAndPlot::mergeWithWorker(
  myproject::Cell&           localCell, 
  const myproject::Cell&     receivedMasterCell,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level
) {
   _map2TimeStepAndPlot2VTKPlotVertexValue_0.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );


   _map2TimeStep.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2InterpolateHangingNodes.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2Inject.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );

}

void myproject::adapters::TimeStepAndPlot::mergeWithWorker(
  myproject::Vertex&        localVertex,
  const myproject::Vertex&  receivedMasterVertex,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {

   _map2TimeStep.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2InterpolateHangingNodes.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2Inject.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );

   _map2TimeStepAndPlot2VTKPlotVertexValue_0.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );

}
#endif

void myproject::adapters::TimeStepAndPlot::touchVertexFirstTime(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2TimeStepAndPlot2VTKPlotVertexValue_0.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


  _map2TimeStep.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2InterpolateHangingNodes.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Inject.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void myproject::adapters::TimeStepAndPlot::touchVertexLastTime(
      myproject::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {

  _map2TimeStep.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2InterpolateHangingNodes.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Inject.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

  _map2TimeStepAndPlot2VTKPlotVertexValue_0.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void myproject::adapters::TimeStepAndPlot::enterCell(
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2TimeStepAndPlot2VTKPlotVertexValue_0.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );


  _map2TimeStep.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2InterpolateHangingNodes.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2Inject.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void myproject::adapters::TimeStepAndPlot::leaveCell(
      myproject::Cell&           fineGridCell,
      myproject::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {

  _map2TimeStep.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2InterpolateHangingNodes.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2Inject.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

  _map2TimeStepAndPlot2VTKPlotVertexValue_0.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void myproject::adapters::TimeStepAndPlot::beginIteration(
  myproject::State&  solverState
) {
  _map2TimeStepAndPlot2VTKPlotVertexValue_0.beginIteration( solverState );


  _map2TimeStep.beginIteration( solverState );
  _map2InterpolateHangingNodes.beginIteration( solverState );
  _map2Inject.beginIteration( solverState );

}


void myproject::adapters::TimeStepAndPlot::endIteration(
  myproject::State&  solverState
) {

  _map2TimeStep.endIteration( solverState );
  _map2InterpolateHangingNodes.endIteration( solverState );
  _map2Inject.endIteration( solverState );

  _map2TimeStepAndPlot2VTKPlotVertexValue_0.endIteration( solverState );

}




void myproject::adapters::TimeStepAndPlot::descend(
  myproject::Cell * const          fineGridCells,
  myproject::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  myproject::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  myproject::Cell&                 coarseGridCell
) {
  _map2TimeStepAndPlot2VTKPlotVertexValue_0.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );


  _map2TimeStep.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2InterpolateHangingNodes.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2Inject.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}


void myproject::adapters::TimeStepAndPlot::ascend(
  myproject::Cell * const    fineGridCells,
  myproject::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  myproject::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  myproject::Cell&           coarseGridCell
) {

  _map2TimeStep.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2InterpolateHangingNodes.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2Inject.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

  _map2TimeStepAndPlot2VTKPlotVertexValue_0.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}
