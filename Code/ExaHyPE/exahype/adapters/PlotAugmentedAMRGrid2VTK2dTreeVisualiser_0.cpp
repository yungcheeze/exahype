#include "exahype/adapters/PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0.h"

#include <sstream>

#include "peano/utils/Loop.h"
#include "peano/grid/CellFlags.h"

#ifdef Parallel
#include "tarch/parallel/Node.h"
#endif


int exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::_snapshotCounter      = 0;
double exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::SqueezeZAxis = 4.0;
double exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::TreeConnectionsValue = -100.0;


peano::CommunicationSpecification   exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::communicationSpecification() {
  return peano::CommunicationSpecification::getPessimisticSpecification();
}


peano::MappingSpecification   exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification   exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial);
}


peano::MappingSpecification   exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial);
}


peano::MappingSpecification   exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification   exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification   exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


std::map<tarch::la::Vector<DIMENSIONS+1,double> , int, tarch::la::VectorCompare<DIMENSIONS+1> >  exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::_vertex2IndexMap;
std::map<tarch::la::Vector<DIMENSIONS+1,double> , int, tarch::la::VectorCompare<DIMENSIONS+1> >  exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::_cellCenter2IndexMap;



exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0():
  _vtkWriter(0),
  _vertexWriter(0),
  _cellWriter(0),
  _cellNumberWriter(0) {
}


exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::~PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0() {
}


#if defined(SharedMemoryParallelisation)
exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0(const PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0&  masterThread):
  _vtkWriter(masterThread._vtkWriter),
  _vertexWriter(masterThread._vertexWriter),
  _cellWriter(masterThread._cellWriter),
  _cellNumberWriter(masterThread._cellNumberWriter) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::mergeWithWorkerThread(const PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0& workerThread) {
}
#endif





void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::plotVertex(
  const exahype::Vertex&                 fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  int                                          level
) {
  tarch::la::Vector<DIMENSIONS+1,double> y;
  for( int i=0; i<DIMENSIONS; i++) y(i) = x(i);
  y(DIMENSIONS)=level;

  tarch::la::Vector<DIMENSIONS+1,double> plotY = y;
  plotY(DIMENSIONS)=level / SqueezeZAxis;
  
  if ( _vertex2IndexMap.find(y) == _vertex2IndexMap.end() ) {  
    _vertex2IndexMap[y] = _vertexWriter->plotVertex(plotY);
  }
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::createHangingVertex(
      exahype::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      exahype::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      exahype::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  plotVertex( fineGridVertex, fineGridX, coarseGridVerticesEnumerator.getLevel()+1 ); 
}



void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::destroyHangingVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::createInnerVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::createBoundaryVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::destroyVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::createCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::destroyCell(
      const exahype::Cell&           fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


#ifdef Parallel
void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::mergeWithNeighbour(
  exahype::Vertex&  vertex,
  const exahype::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::prepareSendToNeighbour(
      exahype::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::prepareCopyToRemoteNode(
      exahype::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::prepareCopyToRemoteNode(
      exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Vertex&  localVertex,
  const exahype::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Cell&  localCell,
  const exahype::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


bool exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::prepareSendToWorker(
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


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::prepareSendToMaster(
      exahype::Cell&                       localCell,
      exahype::Vertex *                    vertices,
      const peano::grid::VertexEnumerator&       verticesEnumerator, 
      const exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      const exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::mergeWithMaster(
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


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::receiveDataFromMaster(
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


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::mergeWithWorker(
      exahype::Cell&           localCell, 
      const exahype::Cell&     receivedMasterCell,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::mergeWithWorker(
      exahype::Vertex&        localVertex,
      const exahype::Vertex&  receivedMasterVertex,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}
#endif


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  plotVertex( fineGridVertex, fineGridX, coarseGridVerticesEnumerator.getLevel()+1 ); 
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::enterCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  int vertexIndex[TWO_POWER_D];
  tarch::la::Vector<DIMENSIONS+1,double> currentVertexPosition;
  currentVertexPosition(DIMENSIONS) = fineGridVerticesEnumerator.getLevel();
  
   dfor2(i)
    for (int d=0; d<DIMENSIONS; d++) {
      currentVertexPosition(d) = fineGridVerticesEnumerator.getVertexPosition(i)(d);
    }
    assertion2 ( _vertex2IndexMap.find(currentVertexPosition) != _vertex2IndexMap.end(), currentVertexPosition, fineGridVertices[ fineGridVerticesEnumerator(i) ].toString() );
    vertexIndex[iScalar] = _vertex2IndexMap[currentVertexPosition];
  enddforx
  
  _cellWriter->plotQuadrangle(vertexIndex);
  
  const int cellIndex = _cellWriter->plotQuadrangle(vertexIndex);
    
  _cellNumberWriter->plotCell(cellIndex,_cellCounter);
  _cellCounter++;

  tarch::la::Vector<DIMENSIONS+1,double> fineCellCenter;
  tarch::la::Vector<DIMENSIONS+1,double> coarseCellCenter;
  for( int i=0; i<DIMENSIONS; i++) fineCellCenter(i)   = fineGridVerticesEnumerator.getCellCenter()(i);
  for( int i=0; i<DIMENSIONS; i++) coarseCellCenter(i) = coarseGridVerticesEnumerator.getCellCenter()(i);
  fineCellCenter(DIMENSIONS)=fineGridVerticesEnumerator.getLevel();
  coarseCellCenter(DIMENSIONS)=fineGridVerticesEnumerator.getLevel()-1;

  tarch::la::Vector<DIMENSIONS+1,double> plotY = fineCellCenter;
  plotY(DIMENSIONS)=(fineCellCenter(DIMENSIONS)-fineGridVerticesEnumerator.getCellSize()(0)/2.0) / SqueezeZAxis;

  _cellCenter2IndexMap[fineCellCenter] = _vertexWriter->plotVertex(plotY);
  const int cellCenterIndex = _cellWriter->plotPoint(_cellCenter2IndexMap[fineCellCenter]);
  _cellNumberWriter->plotCell(cellCenterIndex,TreeConnectionsValue);
  
  if ( _cellCenter2IndexMap.count(coarseCellCenter)>0 ) {
    int indices[2];
    indices[0] = _cellCenter2IndexMap[fineCellCenter];
    indices[1] = _cellCenter2IndexMap[coarseCellCenter];
    const int treeLinkIndex = _cellWriter->plotLine(indices);
    _cellNumberWriter->plotCell(treeLinkIndex,TreeConnectionsValue);
  }
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::leaveCell(
      exahype::Cell&           fineGridCell,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::beginIteration(
  exahype::State&  solverState
) {
  assertion( _vtkWriter==0 );
  
  _cellCounter = 0;
  
  _vtkWriter = new UsedWriter();
  
  _vertexWriter     = _vtkWriter->createVertexWriter();
  _cellWriter       = _vtkWriter->createCellWriter();
  
  _cellNumberWriter = _vtkWriter->createCellDataWriter( "cell-number" ,1);
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::endIteration(
  exahype::State&  solverState
) {
  _vertexWriter->close();
  _cellWriter->close();
  
  _cellNumberWriter->close();
  
  delete _vertexWriter;
  delete _cellWriter;

  delete _cellNumberWriter;
  
  _vertexWriter                  = 0;
  _cellWriter                    = 0;

  _cellNumberWriter              = 0;
  
  std::ostringstream snapshotFileName;
  snapshotFileName << "tree"
                   #ifdef Parallel
                   << "-rank-" << tarch::parallel::Node::getInstance().getRank()
                   #endif
                   << "-" << _snapshotCounter
                   << ".vtk";
  _vtkWriter->writeToFile( snapshotFileName.str() );
  
  _snapshotCounter++;                  
  
  _vertex2IndexMap.clear();
  
  delete _vtkWriter;
  _vtkWriter = 0;
}




void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::descend(
  exahype::Cell * const          fineGridCells,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell
) {
}


void exahype::adapters::PlotAugmentedAMRGrid2VTK2dTreeVisualiser_0::ascend(
  exahype::Cell * const    fineGridCells,
  exahype::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell
) {
}
