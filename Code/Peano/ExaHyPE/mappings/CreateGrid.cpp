#include "ExaHyPE/mappings/CreateGrid.h"



/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification   ExaHyPE::mappings::CreateGrid::communicationSpecification() {
  return peano::CommunicationSpecification(peano::CommunicationSpecification::SendDataAndStateBeforeFirstTouchVertexFirstTime,peano::CommunicationSpecification::SendDataAndStateAfterLastTouchVertexLastTime,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   ExaHyPE::mappings::CreateGrid::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   ExaHyPE::mappings::CreateGrid::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   ExaHyPE::mappings::CreateGrid::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   ExaHyPE::mappings::CreateGrid::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   ExaHyPE::mappings::CreateGrid::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidCoarseGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   ExaHyPE::mappings::CreateGrid::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidCoarseGridRaces);
}


tarch::logging::Log                ExaHyPE::mappings::CreateGrid::_log( "ExaHyPE::mappings::CreateGrid" ); 


ExaHyPE::mappings::CreateGrid::CreateGrid() {
  logTraceIn( "CreateGrid()" );
  // @todo Insert your code here
  logTraceOut( "CreateGrid()" );
}


ExaHyPE::mappings::CreateGrid::~CreateGrid() {
  logTraceIn( "~CreateGrid()" );
  // @todo Insert your code here
  logTraceOut( "~CreateGrid()" );
}


#if defined(SharedMemoryParallelisation)
ExaHyPE::mappings::CreateGrid::CreateGrid(const CreateGrid&  masterThread) {
  logTraceIn( "CreateGrid(CreateGrid)" );
  // @todo Insert your code here
  logTraceOut( "CreateGrid(CreateGrid)" );
}


void ExaHyPE::mappings::CreateGrid::mergeWithWorkerThread(const CreateGrid& workerThread) {
  logTraceIn( "mergeWithWorkerThread(CreateGrid)" );
  // @todo Insert your code here
  logTraceOut( "mergeWithWorkerThread(CreateGrid)" );
}
#endif


void ExaHyPE::mappings::CreateGrid::createHangingVertex(
      ExaHyPE::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      ExaHyPE::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      ExaHyPE::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "createHangingVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "createHangingVertex(...)", fineGridVertex );
}


void ExaHyPE::mappings::CreateGrid::destroyHangingVertex(
      const ExaHyPE::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "destroyHangingVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "destroyHangingVertex(...)", fineGridVertex );
}


void ExaHyPE::mappings::CreateGrid::createInnerVertex(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "createInnerVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "createInnerVertex(...)", fineGridVertex );
}


void ExaHyPE::mappings::CreateGrid::createBoundaryVertex(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "createBoundaryVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "createBoundaryVertex(...)", fineGridVertex );
}


void ExaHyPE::mappings::CreateGrid::destroyVertex(
      const ExaHyPE::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "destroyVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "destroyVertex(...)", fineGridVertex );
}


void ExaHyPE::mappings::CreateGrid::createCell(
      ExaHyPE::Cell&                 fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "createCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
  // @todo Insert your code here
  logTraceOutWith1Argument( "createCell(...)", fineGridCell );
}


void ExaHyPE::mappings::CreateGrid::destroyCell(
      const ExaHyPE::Cell&           fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "destroyCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
  // @todo Insert your code here
  logTraceOutWith1Argument( "destroyCell(...)", fineGridCell );
}

#ifdef Parallel
void ExaHyPE::mappings::CreateGrid::mergeWithNeighbour(
  ExaHyPE::Vertex&  vertex,
  const ExaHyPE::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
  logTraceInWith6Arguments( "mergeWithNeighbour(...)", vertex, neighbour, fromRank, fineGridX, fineGridH, level );
  // @todo Insert your code here
  logTraceOut( "mergeWithNeighbour(...)" );
}

void ExaHyPE::mappings::CreateGrid::prepareSendToNeighbour(
  ExaHyPE::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
  logTraceInWith3Arguments( "prepareSendToNeighbour(...)", vertex, toRank, level );
  // @todo Insert your code here
  logTraceOut( "prepareSendToNeighbour(...)" );
}

void ExaHyPE::mappings::CreateGrid::prepareCopyToRemoteNode(
  ExaHyPE::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
  logTraceInWith5Arguments( "prepareCopyToRemoteNode(...)", localVertex, toRank, x, h, level );
  // @todo Insert your code here
  logTraceOut( "prepareCopyToRemoteNode(...)" );
}

void ExaHyPE::mappings::CreateGrid::prepareCopyToRemoteNode(
  ExaHyPE::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
) {
  logTraceInWith5Arguments( "prepareCopyToRemoteNode(...)", localCell, toRank, cellCentre, cellSize, level );
  // @todo Insert your code here
  logTraceOut( "prepareCopyToRemoteNode(...)" );
}

void ExaHyPE::mappings::CreateGrid::mergeWithRemoteDataDueToForkOrJoin(
  ExaHyPE::Vertex&  localVertex,
  const ExaHyPE::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
  logTraceInWith6Arguments( "mergeWithRemoteDataDueToForkOrJoin(...)", localVertex, masterOrWorkerVertex, fromRank, x, h, level );
  // @todo Insert your code here
  logTraceOut( "mergeWithRemoteDataDueToForkOrJoin(...)" );
}

void ExaHyPE::mappings::CreateGrid::mergeWithRemoteDataDueToForkOrJoin(
  ExaHyPE::Cell&  localCell,
  const ExaHyPE::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                       level
) {
  logTraceInWith3Arguments( "mergeWithRemoteDataDueToForkOrJoin(...)", localCell, masterOrWorkerCell, fromRank );
  // @todo Insert your code here
  logTraceOut( "mergeWithRemoteDataDueToForkOrJoin(...)" );
}

bool ExaHyPE::mappings::CreateGrid::prepareSendToWorker(
  ExaHyPE::Cell&                 fineGridCell,
  ExaHyPE::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  ExaHyPE::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  ExaHyPE::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  logTraceIn( "prepareSendToWorker(...)" );
  // @todo Insert your code here
  logTraceOutWith1Argument( "prepareSendToWorker(...)", true );
  return true;
}

void ExaHyPE::mappings::CreateGrid::prepareSendToMaster(
  ExaHyPE::Cell&                       localCell,
  ExaHyPE::Vertex *                    vertices,
  const peano::grid::VertexEnumerator&       verticesEnumerator, 
  const ExaHyPE::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
  const ExaHyPE::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
  logTraceInWith2Arguments( "prepareSendToMaster(...)", localCell, verticesEnumerator.toString() );
  // @todo Insert your code here
  logTraceOut( "prepareSendToMaster(...)" );
}


void ExaHyPE::mappings::CreateGrid::mergeWithMaster(
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
  logTraceIn( "mergeWithMaster(...)" );
  // @todo Insert your code here
  logTraceOut( "mergeWithMaster(...)" );
}


void ExaHyPE::mappings::CreateGrid::receiveDataFromMaster(
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
  logTraceIn( "receiveDataFromMaster(...)" );
  // @todo Insert your code here
  logTraceOut( "receiveDataFromMaster(...)" );
}


void ExaHyPE::mappings::CreateGrid::mergeWithWorker(
  ExaHyPE::Cell&           localCell, 
  const ExaHyPE::Cell&     receivedMasterCell,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level
) {
  logTraceInWith2Arguments( "mergeWithWorker(...)", localCell.toString(), receivedMasterCell.toString() );
  // @todo Insert your code here
  logTraceOutWith1Argument( "mergeWithWorker(...)", localCell.toString() );
}


void ExaHyPE::mappings::CreateGrid::mergeWithWorker(
  ExaHyPE::Vertex&        localVertex,
  const ExaHyPE::Vertex&  receivedMasterVertex,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
  logTraceInWith2Arguments( "mergeWithWorker(...)", localVertex.toString(), receivedMasterVertex.toString() );
  // @todo Insert your code here
  logTraceOutWith1Argument( "mergeWithWorker(...)", localVertex.toString() );
}
#endif

void ExaHyPE::mappings::CreateGrid::touchVertexFirstTime(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "touchVertexFirstTime(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "touchVertexFirstTime(...)", fineGridVertex );
}


void ExaHyPE::mappings::CreateGrid::touchVertexLastTime(
      ExaHyPE::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "touchVertexLastTime(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // @todo Insert your code here
  logTraceOutWith1Argument( "touchVertexLastTime(...)", fineGridVertex );
}


void ExaHyPE::mappings::CreateGrid::enterCell(
      ExaHyPE::Cell&                 fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "enterCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
  // @todo Insert your code here
  logTraceOutWith1Argument( "enterCell(...)", fineGridCell );
}


void ExaHyPE::mappings::CreateGrid::leaveCell(
      ExaHyPE::Cell&           fineGridCell,
      ExaHyPE::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "leaveCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
  // @todo Insert your code here
  logTraceOutWith1Argument( "leaveCell(...)", fineGridCell );
}


void ExaHyPE::mappings::CreateGrid::beginIteration(
  ExaHyPE::State&  solverState
) {
  logTraceInWith1Argument( "beginIteration(State)", solverState );
  // @todo Insert your code here
  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


void ExaHyPE::mappings::CreateGrid::endIteration(
  ExaHyPE::State&  solverState
) {
  logTraceInWith1Argument( "endIteration(State)", solverState );
  // @todo Insert your code here
  logTraceOutWith1Argument( "endIteration(State)", solverState);
}



void ExaHyPE::mappings::CreateGrid::descend(
  ExaHyPE::Cell * const          fineGridCells,
  ExaHyPE::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  ExaHyPE::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  ExaHyPE::Cell&                 coarseGridCell
) {
  logTraceInWith2Arguments( "descend(...)", coarseGridCell.toString(), coarseGridVerticesEnumerator.toString() );
  // @todo Insert your code here
  logTraceOut( "descend(...)" );
}


void ExaHyPE::mappings::CreateGrid::ascend(
  ExaHyPE::Cell * const    fineGridCells,
  ExaHyPE::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  ExaHyPE::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  ExaHyPE::Cell&           coarseGridCell
) {
  logTraceInWith2Arguments( "ascend(...)", coarseGridCell.toString(), coarseGridVerticesEnumerator.toString() );
  // @todo Insert your code here
  logTraceOut( "ascend(...)" );
}
