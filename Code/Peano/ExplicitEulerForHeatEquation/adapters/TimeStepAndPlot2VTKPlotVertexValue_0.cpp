#include "ExplicitEulerForHeatEquation/adapters/TimeStepAndPlot2VTKPlotVertexValue_0.h"

#include <sstream>

#include "peano/utils/Loop.h"

#ifdef Parallel
#include "tarch/parallel/Node.h"
#endif


int myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::_snapshotCounter = 0;


peano::CommunicationSpecification   myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::communicationSpecification() {
  return peano::CommunicationSpecification::getPessimisticSpecification();
}


peano::MappingSpecification   myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification   myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial);
}


peano::MappingSpecification   myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial);
}


peano::MappingSpecification   myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification   myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification   myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


std::map<tarch::la::Vector<DIMENSIONS,double> , int, tarch::la::VectorCompare<DIMENSIONS> >  myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::_vertex2IndexMap;


myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::TimeStepAndPlot2VTKPlotVertexValue_0():
  _vtkWriter(0),
  _vertexWriter(0),
  _cellWriter(0),
  _vertexValueWriter(0) {
}


myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::~TimeStepAndPlot2VTKPlotVertexValue_0() {
}


#if defined(SharedMemoryParallelisation)
myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::TimeStepAndPlot2VTKPlotVertexValue_0(const TimeStepAndPlot2VTKPlotVertexValue_0&  masterThread):
  _vtkWriter(masterThread._vtkWriter),
  _vertexWriter(masterThread._vertexWriter),
  _cellWriter(masterThread._cellWriter),
  _vertexValueWriter(masterThread._vertexValueWriter) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::mergeWithWorkerThread(const TimeStepAndPlot2VTKPlotVertexValue_0& workerThread) {
}
#endif





void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::plotVertex(
  const myproject::Vertex&                 fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&  fineGridX
) {
  if ( 
    fineGridVertex.getRefinementControl() != myproject::Vertex::Records::Refined &&
    fineGridVertex.getRefinementControl() != myproject::Vertex::Records::Refining &&
    _vertex2IndexMap.find(fineGridX) == _vertex2IndexMap.end() 
  ) {  
    #if defined(Dim2) || defined(Dim3)
    _vertex2IndexMap[fineGridX] = _vertexWriter->plotVertex(fineGridX);
    #else
    _vertex2IndexMap[fineGridX] = _vertexWriter->plotVertex(tarch::la::Vector<3,double>(fineGridX.data()));
    #endif
    
    #ifdef Parallel
    const bool plotValue =
      !fineGridVertex.isOutside() &&
      fineGridVertex.isAdjacentToDomainOf( tarch::parallel::Node::getInstance().getRank() );
    #else
    const bool plotValue = true;
    #endif

    const double value = plotValue ? fineGridVertex.getU() : 0.0;

    _vertexValueWriter->plotVertex             (_vertex2IndexMap[fineGridX],value );
  }
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::createHangingVertex(
      myproject::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      myproject::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      myproject::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  plotVertex( fineGridVertex, fineGridX ); 
}



void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::destroyHangingVertex(
      const myproject::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::createInnerVertex(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::createBoundaryVertex(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::destroyVertex(
      const myproject::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::createCell(
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::destroyCell(
      const myproject::Cell&           fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


#ifdef Parallel
void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::mergeWithNeighbour(
  myproject::Vertex&  vertex,
  const myproject::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::prepareSendToNeighbour(
      myproject::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::prepareCopyToRemoteNode(
      myproject::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::prepareCopyToRemoteNode(
      myproject::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::mergeWithRemoteDataDueToForkOrJoin(
  myproject::Vertex&  localVertex,
  const myproject::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::mergeWithRemoteDataDueToForkOrJoin(
  myproject::Cell&  localCell,
  const myproject::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


bool myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::prepareSendToWorker(
  myproject::Cell&                 fineGridCell,
  myproject::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  myproject::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  myproject::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  return false;
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::prepareSendToMaster(
      myproject::Cell&                       localCell,
      myproject::Vertex *                    vertices,
      const peano::grid::VertexEnumerator&       verticesEnumerator, 
      const myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      const myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::mergeWithMaster(
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
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::receiveDataFromMaster(
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
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::mergeWithWorker(
      myproject::Cell&           localCell, 
      const myproject::Cell&     receivedMasterCell,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::mergeWithWorker(
      myproject::Vertex&        localVertex,
      const myproject::Vertex&  receivedMasterVertex,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}
#endif


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::touchVertexFirstTime(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  plotVertex( fineGridVertex, fineGridX ); 
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::touchVertexLastTime(
      myproject::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::enterCell(
  myproject::Cell&                 fineGridCell,
  myproject::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  myproject::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  myproject::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&            fineGridPositionOfCell
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
      assertion4 ( 
        _vertex2IndexMap.find(currentVertexPosition) != _vertex2IndexMap.end(), 
        currentVertexPosition, 
        fineGridCell.toString(),
        fineGridVerticesEnumerator.toString(),
        fineGridVertices[ fineGridVerticesEnumerator(i) ].toString()
      );
      vertexIndex[iScalar] = _vertex2IndexMap[currentVertexPosition];
    enddforx
  
    if (DIMENSIONS==2) {
      _cellWriter->plotQuadrangle(vertexIndex);
    }
    if (DIMENSIONS==3) {
      _cellWriter->plotHexahedron(vertexIndex);
    }
  }
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::leaveCell(
      myproject::Cell&           fineGridCell,
      myproject::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::beginIteration(
  myproject::State&  solverState
) {
  assertion( _vtkWriter==0 );
  
  _vtkWriter = new tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter();
  
  _vertexWriter     = _vtkWriter->createVertexWriter();
  _cellWriter       = _vtkWriter->createCellWriter();
  
  _vertexValueWriter               = _vtkWriter->createVertexDataWriter("u",1);
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::endIteration(
  myproject::State&  solverState
) {
  _vertexWriter->close();
  _cellWriter->close();
  _vertexValueWriter->close();
  
  delete _vertexWriter;
  delete _cellWriter;
  delete _vertexValueWriter;
  
  _vertexWriter                  = 0;
  _cellWriter                    = 0;
  _vertexValueWriter             = 0;
  
  std::ostringstream snapshotFileName;
  snapshotFileName << "result"
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




void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::descend(
  myproject::Cell * const          fineGridCells,
  myproject::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  myproject::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  myproject::Cell&                 coarseGridCell
) {
}


void myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0::ascend(
  myproject::Cell * const    fineGridCells,
  myproject::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  myproject::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  myproject::Cell&           coarseGridCell
) {
}
