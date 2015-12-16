#include "EulerFlow3d/mappings/VTKExport.h"

#include "EulerFlow3d/Constants.h"

#include "EulerFlow3d/quad/GaussLegendre.h"

#include "EulerFlow3d/geometry/Mapping.h"

#include "EulerFlow3d/problem/Problem.h"

#include "EulerFlow3d/dg/DGMatrices.h"

#include "string.h"

int exahype::mappings::VTKExport::_snapshotCounter = 0;

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification   exahype::mappings::VTKExport::communicationSpecification() {
  return peano::CommunicationSpecification(peano::CommunicationSpecification::SendDataAndStateBeforeFirstTouchVertexFirstTime,peano::CommunicationSpecification::SendDataAndStateAfterLastTouchVertexLastTime,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::VTKExport::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::VTKExport::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::VTKExport::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::OnlyLeaves,peano::MappingSpecification::Serial);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::VTKExport::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::VTKExport::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::VTKExport::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


tarch::logging::Log                exahype::mappings::VTKExport::_log( "exahype::mappings::VTKExport" ); 

exahype::mappings::VTKExport::VTKExport():
          _vtkWriter(0),
          _vertexWriter(0),
          _cellWriter(0),
          _vertexValueWriter(0)
{
  // do nothing
}


exahype::mappings::VTKExport::~VTKExport() {
  // do nothing
}


#if defined(SharedMemoryParallelisation)
exahype::mappings::VTKExport::VTKExport(const VTKExport&  masterThread):
          _vtkWriter(masterThread._vtkWriter),
          _vertexWriter(masterThread._vertexWriter),
          _cellWriter(masterThread._cellWriter),
          _vertexValueWriter(masterThread._vertexValueWriter)
//      _cellValueWriter(masterThread._cellValueWriter)
  {
}


void exahype::mappings::VTKExport::mergeWithWorkerThread(const VTKExport& workerThread) {
}
#endif

void exahype::mappings::VTKExport::createHangingVertex(
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


void exahype::mappings::VTKExport::destroyHangingVertex(
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


void exahype::mappings::VTKExport::createInnerVertex(
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


void exahype::mappings::VTKExport::createBoundaryVertex(
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


void exahype::mappings::VTKExport::destroyVertex(
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


void exahype::mappings::VTKExport::createCell(
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


void exahype::mappings::VTKExport::destroyCell(
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
void exahype::mappings::VTKExport::mergeWithNeighbour(
    exahype::Vertex&  vertex,
    const exahype::Vertex&  neighbour,
    int                                           fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::VTKExport::prepareSendToNeighbour(
    exahype::Vertex&  vertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::VTKExport::prepareCopyToRemoteNode(
    exahype::Vertex&  localVertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::VTKExport::prepareCopyToRemoteNode(
    exahype::Cell&  localCell,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::VTKExport::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex&  localVertex,
    const exahype::Vertex&  masterOrWorkerVertex,
    int                                       fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&  x,
    const tarch::la::Vector<DIMENSIONS,double>&  h,
    int                                       level
) {
  // do nothing
}

void exahype::mappings::VTKExport::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell&  localCell,
    const exahype::Cell&  masterOrWorkerCell,
    int                                       fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                       level
) {
  // do nothing
}

bool exahype::mappings::VTKExport::prepareSendToWorker(
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
    int                                                                  worker
) {
  logTraceIn( "prepareSendToWorker(...)" );
  return false;
  logTraceOutWith1Argument( "prepareSendToWorker(...)", true );
  return true;
}

void exahype::mappings::VTKExport::prepareSendToMaster(
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


void exahype::mappings::VTKExport::mergeWithMaster(
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


void exahype::mappings::VTKExport::receiveDataFromMaster(
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


void exahype::mappings::VTKExport::mergeWithWorker(
    exahype::Cell&           localCell,
    const exahype::Cell&     receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                          level
) {
  // do nothing
}


void exahype::mappings::VTKExport::mergeWithWorker(
    exahype::Vertex&        localVertex,
    const exahype::Vertex&  receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}
#endif

void exahype::mappings::VTKExport::touchVertexFirstTime(
    exahype::Vertex&               fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "touchVertexFirstTime(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );

  if (fineGridVertex.getRefinementControl() != exahype::Vertex::Records::Refined &&
      fineGridVertex.getRefinementControl() != exahype::Vertex::Records::Refining) {
    //    plotVertex( fineGridVertex, fineGridX );
  }

  logTraceOutWith1Argument( "touchVertexFirstTime(...)", fineGridVertex );
}


void exahype::mappings::VTKExport::touchVertexLastTime(
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

void exahype::mappings::VTKExport::enterCell(
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "enterCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );

#ifdef Parallel
  if (!fineGridCell.isRefined() && !fineGridCell.isAssignedToRemoteRank()) {
#else
    if (!fineGridCell.isRefined()) {
#endif
      assertion( DIMENSIONS==2);

      records::CellDescription& cellDescription =
          CellDescriptionHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[0];

      const tarch::la::Vector<DIMENSIONS,double> center = fineGridVerticesEnumerator.getCellCenter();  // the center of the cell
      const double dx = fineGridVerticesEnumerator.getCellSize()(0);
      const double dy = fineGridVerticesEnumerator.getCellSize()(1);

      const double dxPatch = dx/EXAHYPE_PATCH_SIZE_X;
      const double dyPatch = dy/EXAHYPE_PATCH_SIZE_Y;

      const int basisSize = EXAHYPE_ORDER+1;
      const int nvar      = EXAHYPE_NVARS;

      // helper variables
      double x = 0;
      double y = 0;

      // BEGIN: move into static fields todo
      int    indexMapping      [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];
      double uniformCoordinates[DIMENSIONS     ][EXAHYPE_NBASIS_POWER_DIMENSIONS];
      double uniformPartition  [TWO_POWER_D    ][EXAHYPE_NBASIS_POWER_DIMENSIONS];
      double uniformDoF        [EXAHYPE_NVARS][EXAHYPE_NBASIS_POWER_DIMENSIONS];

      // define sub nodes
      for (int ii=0; ii < basisSize; ii++) {
        for (int jj=0; jj < basisSize; jj++) {
          const int uniformNodeIndex = ii + basisSize * jj;

          indexMapping[ii][jj] = uniformNodeIndex;
          uniformCoordinates[0][uniformNodeIndex] = (double) ii / (double) EXAHYPE_ORDER;
          uniformCoordinates[1][uniformNodeIndex] = (double) jj / (double) EXAHYPE_ORDER;
        }
      }

      // define sub quadrangles/hexahedrons
      for (int ii=0; ii < basisSize; ii++) {
        for (int jj=0; jj < basisSize; jj++) {
          const int uniformNodeIndex = ii + basisSize * jj;

          uniformPartition[0][uniformNodeIndex] = indexMapping[ii  ][jj  ];
          uniformPartition[1][uniformNodeIndex] = indexMapping[ii+1][jj  ];
          uniformPartition[2][uniformNodeIndex] = indexMapping[ii+1][jj+1];
          uniformPartition[3][uniformNodeIndex] = indexMapping[ii  ][jj+1];
        }
      }

      // END: move into static fields todo

      for (int i=1; i<EXAHYPE_PATCH_SIZE_X+1; i++) { // loop over patches
        for (int j=1; j<EXAHYPE_PATCH_SIZE_Y+1; j++) {
          const int patchIndex = i + (EXAHYPE_PATCH_SIZE_X+2) * j;

          // Map Gauss-Legendre nodes to equidistant subgrid coordinates
//          std::memset((double *) &subData[0],0,sizeof(double) * EXAHYPE_NVARS * EXAHYPE_NBASIS_POWER_DIMENSIONS);
          for (int ii=0; ii<basisSize; ii++) { // mem zero
             for (int jj=0; jj<basisSize; jj++) {
               const int uniformNodeIndex = ii + basisSize * jj;

               for (int ivar=0; ivar < nvar; ivar++) {
                 uniformDoF[ivar][uniformNodeIndex] = 0;
               }
             }
          }

          double* luh = &(DataHeap::getInstance().getData(cellDescription.getSolution(patchIndex))[0]._persistentRecords._u);

          for (int ii=0; ii<basisSize; ii++) { // project on subgrid coordinates
            for (int jj=0; jj<basisSize; jj++) {
              const int uniformNodeIndex = ii + basisSize * jj;

              for (int mm=0; mm<basisSize; mm++) { // project on subgrid coordinates
                for (int nn=0; nn<basisSize; nn++) {
                  const int nodeIndex     = mm + basisSize * nn;
                  const int dofStartIndex = nodeIndex * nvar;

                  for (int ivar=0; ivar < nvar; ivar++) {
                    uniformDoF[ivar][uniformNodeIndex] += luh[dofStartIndex+ivar] * dg::subOutputMatrix[nodeIndex][uniformNodeIndex];
                  }
                }
              }
            }
          }

          for (int ii=0; ii<basisSize; ii++) { // loop over dof
            for (int jj=0; jj<basisSize; jj++) {
              // location and index of nodal degrees of freedom
              const int nodeIndex     = ii + basisSize * jj;
              const int dofStartIndex = nodeIndex * nvar;

              const double r = uniformCoordinates[0][nodeIndex];
              const double s = uniformCoordinates[1][nodeIndex];

//              const double r = quad::gaussLegendreNodes[basisSize-1][ii];
//              const double s = quad::gaussLegendreNodes[basisSize-1][jj];

              geometry::mapping2d(center(0),center(1),dx,dy,dxPatch,dyPatch,i,j,r,s,&x,&y);
              tarch::la::Vector<DIMENSIONS,double> currentVertexPosition(x,y);

              const int vtkNodeIndex = _vertexWriter->plotVertex(currentVertexPosition);
              _cellWriter->plotPoint(vtkNodeIndex);

              for (int ivar=4; ivar < 5; ivar++) {
//                const double dofValue = u[dofStartIndex+ivar];
                const double dofValue = uniformDoF[ivar][nodeIndex];
                _vertexValueWriter->plotVertex(vtkNodeIndex,dofValue);
              }
            }
          }
        }
      }
    }
    logTraceOutWith1Argument( "enterCell(...)", fineGridCell );
  }


  void exahype::mappings::VTKExport::leaveCell(
      exahype::Cell&           fineGridCell,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
  ) {
    logTraceInWith4Arguments( "leaveCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
    // do nothing
    logTraceOutWith1Argument( "leaveCell(...)", fineGridCell );
  }


  void exahype::mappings::VTKExport::beginIteration(
      exahype::State&  solverState
  ) {
    logTraceInWith1Argument( "beginIteration(State)", solverState );

    assertion( _vtkWriter==0 );
    _vtkWriter         = new tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter();
    _vertexWriter      = _vtkWriter->createVertexWriter();
    _cellWriter        = _vtkWriter->createCellWriter();

    _vertexValueWriter = _vtkWriter->createVertexDataWriter("Qh",1);
//    _cellValueWriter   = _vtkWriter->createCellDataWriter("concentration",1);

    logTraceOutWith1Argument( "beginIteration(State)", solverState);
  }


  void exahype::mappings::VTKExport::endIteration(
      exahype::State&  solverState
  ) {
    logTraceInWith1Argument( "endIteration(State)", solverState );

    _vertexWriter->close();
    _cellWriter->close();
    _vertexValueWriter->close();
    //  _cellValueWriter->close();

    delete _vertexWriter;
    delete _vertexValueWriter;
    delete _cellWriter;
    //  delete _cellValueWriter;

    _vertexWriter      = nullptr;
    _cellWriter        = nullptr;
    _vertexValueWriter = nullptr;
    //  _cellValueWriter = nullptr;

    std::ostringstream snapshotFileName;
    snapshotFileName << "solution"
#ifdef Parallel
        << "-rank-" << tarch::parallel::Node::getInstance().getRank()
#endif
        << "-" << _snapshotCounter
        << ".vtk";
    _vtkWriter->writeToFile( snapshotFileName.str() );

    _snapshotCounter++;

    delete _vtkWriter;
    _vtkWriter = 0;

    logTraceOutWith1Argument( "endIteration(State)", solverState);
  }

  void exahype::mappings::VTKExport::descend(
      exahype::Cell * const          fineGridCells,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell
  ) {
    logTraceInWith2Arguments( "descend(...)", coarseGridCell.toString(), coarseGridVerticesEnumerator.toString() );
    // do nothing
    logTraceOut( "descend(...)" );
  }


  void exahype::mappings::VTKExport::ascend(
      exahype::Cell * const    fineGridCells,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell
  ) {
    logTraceInWith2Arguments( "ascend(...)", coarseGridCell.toString(), coarseGridVerticesEnumerator.toString() );
    // do nothing
    logTraceOut( "ascend(...)" );
  }
