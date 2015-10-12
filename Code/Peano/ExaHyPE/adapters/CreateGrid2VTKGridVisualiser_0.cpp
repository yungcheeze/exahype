#include "ExaHyPE/adapters/CreateGrid2VTKGridVisualiser_0.h"

#include <sstream>

#include "peano/utils/Loop.h"
#include "peano/grid/CellFlags.h"

#ifdef Parallel
#include "tarch/parallel/Node.h"
#endif


int ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::_snapshotCounter = 0;



peano::CommunicationSpecification   ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::communicationSpecification() {
  return peano::CommunicationSpecification::getPessimisticSpecification();
}


peano::MappingSpecification   ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification   ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial);
}


peano::MappingSpecification   ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification   ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial);
}


peano::MappingSpecification   ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification   ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


std::map<tarch::la::Vector<DIMENSIONS,double> , int, tarch::la::VectorCompare<DIMENSIONS> >  ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::_vertex2IndexMap;


ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::CreateGrid2VTKGridVisualiser_0():
  _vtkWriter(0),
  _vertexWriter(0),
  _cellWriter(0),
  _vertexTypeWriter(0),
  _vertexRefinementControlWriter(0),
  _vertexAdjacentCellsHeight(0),
  _cellStateWriter(0) {
}


ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::~CreateGrid2VTKGridVisualiser_0() {
}


#if defined(SharedMemoryParallelisation)
ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::CreateGrid2VTKGridVisualiser_0(const CreateGrid2VTKGridVisualiser_0&  masterThread):
  _vtkWriter(masterThread._vtkWriter),
  _vertexWriter(masterThread._vertexWriter),
  _cellWriter(masterThread._cellWriter),
  _vertexTypeWriter(masterThread._vertexTypeWriter),
  _vertexRefinementControlWriter(masterThread._vertexRefinementControlWriter),
  _vertexAdjacentCellsHeight(masterThread._vertexAdjacentCellsHeight),
  _cellStateWriter(masterThread._cellStateWriter)
{
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::mergeWithWorkerThread(const CreateGrid2VTKGridVisualiser_0& workerThread) {
}
#endif





void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::plotVertex(
  const ExaHyPE::Vertex&                 fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&  fineGridX
) {
  if ( _vertex2IndexMap.find(fineGridX) == _vertex2IndexMap.end() ) {
    assertion( _vertexWriter                  != nullptr );
    assertion( _vertexTypeWriter              != nullptr );
    assertion( _vertexRefinementControlWriter != nullptr );
    assertion( _vertexAdjacentCellsHeight     != nullptr );
    
    #if defined(Dim2) || defined(Dim3)
    _vertex2IndexMap[fineGridX] = _vertexWriter->plotVertex(fineGridX);
    #else
    _vertex2IndexMap[fineGridX] = _vertexWriter->plotVertex(tarch::la::Vector<3,double>(fineGridX.data()));
    #endif
    const int boundaryFlag = fineGridVertex.isHangingNode() ? -1 : fineGridVertex.isBoundary() ? 1 : 0;
    _vertexTypeWriter->plotVertex             (_vertex2IndexMap[fineGridX], boundaryFlag);
    _vertexRefinementControlWriter->plotVertex(_vertex2IndexMap[fineGridX],fineGridVertex.getRefinementControl() );
    _vertexAdjacentCellsHeight->plotVertex    (_vertex2IndexMap[fineGridX],fineGridVertex.getAdjacentCellsHeightOfPreviousIteration() );
  }
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::createHangingVertex(
      ExaHyPE::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      ExaHyPE::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      ExaHyPE::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  plotVertex( fineGridVertex, fineGridX ); 
}



void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::destroyHangingVertex(
      const ExaHyPE::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::createInnerVertex(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::createBoundaryVertex(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::destroyVertex(
      const ExaHyPE::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::createCell(
      ExaHyPE::Cell&                 fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::destroyCell(
      const ExaHyPE::Cell&           fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


#ifdef Parallel
void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::mergeWithNeighbour(
  ExaHyPE::Vertex&  vertex,
  const ExaHyPE::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::prepareSendToNeighbour(
      ExaHyPE::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::prepareCopyToRemoteNode(
      ExaHyPE::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::prepareCopyToRemoteNode(
      ExaHyPE::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::mergeWithRemoteDataDueToForkOrJoin(
  ExaHyPE::Vertex&  localVertex,
  const ExaHyPE::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::mergeWithRemoteDataDueToForkOrJoin(
  ExaHyPE::Cell&  localCell,
  const ExaHyPE::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


bool ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::prepareSendToWorker(
  ExaHyPE::Cell&                 fineGridCell,
  ExaHyPE::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  ExaHyPE::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  ExaHyPE::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  return false;
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::prepareSendToMaster(
      ExaHyPE::Cell&                       localCell,
      ExaHyPE::Vertex *                    vertices,
      const peano::grid::VertexEnumerator&       verticesEnumerator, 
      const ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      const ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::mergeWithMaster(
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
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::receiveDataFromMaster(
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
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::mergeWithWorker(
      ExaHyPE::Cell&           localCell, 
      const ExaHyPE::Cell&     receivedMasterCell,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::mergeWithWorker(
      ExaHyPE::Vertex&        localVertex,
      const ExaHyPE::Vertex&  receivedMasterVertex,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}
#endif


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::touchVertexFirstTime(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  if (
    fineGridVertex.getRefinementControl()==ExaHyPE::Vertex::Records::Unrefined ||
    fineGridVertex.getRefinementControl()==ExaHyPE::Vertex::Records::RefinementTriggered ||
    fineGridVertex.getRefinementControl()==ExaHyPE::Vertex::Records::Erasing
  ) {
    plotVertex( fineGridVertex, fineGridX ); 
  }
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::touchVertexLastTime(
      ExaHyPE::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::enterCell(
      ExaHyPE::Cell&                 fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::leaveCell(
      ExaHyPE::Cell&           fineGridCell,
      ExaHyPE::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  #ifdef Parallel
  if (fineGridCell.isLeaf() && !fineGridCell.isAssignedToRemoteRank()) {
  #else
  if (fineGridCell.isLeaf()) {
  #endif
    assertion( DIMENSIONS==2 || DIMENSIONS==3 );
    int vertexIndex[TWO_POWER_D];
     dfor2(i)
      tarch::la::Vector<DIMENSIONS,double> currentVertexPosition = fineGridVerticesEnumerator.getVertexPosition(i);
      assertion2 ( _vertex2IndexMap.find(currentVertexPosition) != _vertex2IndexMap.end(), currentVertexPosition, fineGridVertices[ fineGridVerticesEnumerator(i) ].toString() );
      vertexIndex[iScalar] = _vertex2IndexMap[currentVertexPosition];
    enddforx
  
    int cellIndex;
    if (DIMENSIONS==2) {
      cellIndex = _cellWriter->plotQuadrangle(vertexIndex);
    }
    if (DIMENSIONS==3) {
      cellIndex = _cellWriter->plotHexahedron(vertexIndex);
    }
    
    _cellStateWriter->plotCell(cellIndex,fineGridVerticesEnumerator.getCellFlags());
  }
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::beginIteration(
  ExaHyPE::State&  solverState
) {
  assertion( _vtkWriter==0 );
  
  _vtkWriter = new tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter();
  
  _vertexWriter     = _vtkWriter->createVertexWriter();
  _cellWriter       = _vtkWriter->createCellWriter();
  
  _vertexTypeWriter               = _vtkWriter->createVertexDataWriter(ExaHyPE::Vertex::Records::getInsideOutsideDomainMapping()+"/Hanging=-1" ,1);
  _vertexRefinementControlWriter  = _vtkWriter->createVertexDataWriter(ExaHyPE::Vertex::Records::getRefinementControlMapping() ,1);
  _vertexAdjacentCellsHeight      = _vtkWriter->createVertexDataWriter( peano::grid::getCellFlagsLegend(),1);

  _cellStateWriter                = _vtkWriter->createCellDataWriter( "cell-flag(>=-1=stationary,-1=parallel-boundary,<=-2=not-stationary" ,1);
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::endIteration(
  ExaHyPE::State&  solverState
) {
  _vertexWriter->close();
  _cellWriter->close();
  
  _vertexTypeWriter->close();
  _vertexRefinementControlWriter->close();
  _vertexAdjacentCellsHeight->close();
  _cellStateWriter->close();
  
  delete _vertexWriter;
  delete _cellWriter;
  delete _vertexTypeWriter;
  delete _vertexRefinementControlWriter;
  delete _vertexAdjacentCellsHeight;
  delete _cellStateWriter;
  
  _vertexWriter                  = nullptr;
  _cellWriter                    = nullptr;
  _vertexTypeWriter              = nullptr;
  _vertexRefinementControlWriter = nullptr;
  _vertexAdjacentCellsHeight     = nullptr;
  _cellStateWriter               = nullptr;
  
  std::ostringstream snapshotFileName;
  snapshotFileName << "grid"
                   #ifdef Parallel
                   << "-rank-" << tarch::parallel::Node::getInstance().getRank()
                   #endif
                   << "-" << _snapshotCounter
                   << ".vtk";
  _vtkWriter->writeToFile( snapshotFileName.str() );
  
  _snapshotCounter++;                  
  
  _vertex2IndexMap.clear();
  
  delete _vtkWriter;
  _vtkWriter = nullptr;
}




void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::descend(
  ExaHyPE::Cell * const          fineGridCells,
  ExaHyPE::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  ExaHyPE::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  ExaHyPE::Cell&                 coarseGridCell
) {
}


void ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0::ascend(
  ExaHyPE::Cell * const    fineGridCells,
  ExaHyPE::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  ExaHyPE::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  ExaHyPE::Cell&           coarseGridCell
) {
}
