#include "EulerFlow3d/mappings/InitialCondition.h"

#include "EulerFlow3d/Constants.h"

#include "EulerFlow3d/math/quad/Gausslegendre.h"

#include "EulerFlow3d/geometry/Mapping.h"

#include "EulerFlow3d/problem/Problem.h"

#include "stdlib.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification   exahype::mappings::InitialCondition::communicationSpecification() {
  return peano::CommunicationSpecification(peano::CommunicationSpecification::SendDataAndStateBeforeFirstTouchVertexFirstTime,peano::CommunicationSpecification::SendDataAndStateAfterLastTouchVertexLastTime,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::InitialCondition::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::InitialCondition::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::InitialCondition::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::InitialCondition::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::InitialCondition::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidCoarseGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::InitialCondition::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidCoarseGridRaces);
}


tarch::logging::Log                exahype::mappings::InitialCondition::_log( "exahype::mappings::InitialCondition" ); 


exahype::mappings::InitialCondition::InitialCondition() {
  logTraceIn( "InitialCondition()" );
  // do nothing
  logTraceOut( "InitialCondition()" );
}


exahype::mappings::InitialCondition::~InitialCondition() {
  logTraceIn( "~InitialCondition()" );
  // do nothing
  logTraceOut( "~InitialCondition()" );
}


#if defined(SharedMemoryParallelisation)
exahype::mappings::InitialCondition::InitialCondition(const InitialCondition&  masterThread) {
  logTraceIn( "InitialCondition(InitialCondition)" );
  // do nothing
  logTraceOut( "InitialCondition(InitialCondition)" );
}


void exahype::mappings::InitialCondition::mergeWithWorkerThread(const InitialCondition& workerThread) {
  logTraceIn( "mergeWithWorkerThread(InitialCondition)" );
  // do nothing
  logTraceOut( "mergeWithWorkerThread(InitialCondition)" );
}
#endif


void exahype::mappings::InitialCondition::createHangingVertex(
    exahype::Vertex&     fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
    exahype::Vertex * const   coarseGridVertices,
    const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
    exahype::Cell&       coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "createHangingVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // do nothing
  logTraceOutWith1Argument( "createHangingVertex(...)", fineGridVertex );
}


void exahype::mappings::InitialCondition::destroyHangingVertex(
    const exahype::Vertex&   fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "destroyHangingVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // do nothing
  logTraceOutWith1Argument( "destroyHangingVertex(...)", fineGridVertex );
}


void exahype::mappings::InitialCondition::createInnerVertex(
    exahype::Vertex&               fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "createInnerVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // do nothing
  logTraceOutWith1Argument( "createInnerVertex(...)", fineGridVertex );
}


void exahype::mappings::InitialCondition::createBoundaryVertex(
    exahype::Vertex&               fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "createBoundaryVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // do nothing
  logTraceOutWith1Argument( "createBoundaryVertex(...)", fineGridVertex );
}


void exahype::mappings::InitialCondition::destroyVertex(
    const exahype::Vertex&   fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "destroyVertex(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // do nothing
  logTraceOutWith1Argument( "destroyVertex(...)", fineGridVertex );
}


void exahype::mappings::InitialCondition::createCell(
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "createCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
  // do nothing
  logTraceOutWith1Argument( "createCell(...)", fineGridCell );
}


void exahype::mappings::InitialCondition::destroyCell(
    const exahype::Cell&           fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "destroyCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );
  // do nothing
  logTraceOutWith1Argument( "destroyCell(...)", fineGridCell );
}

#ifdef Parallel
void exahype::mappings::InitialCondition::mergeWithNeighbour(
    exahype::Vertex&  vertex,
    const exahype::Vertex&  neighbour,
    int                                           fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
    int                                           level
) {
  logTraceInWith6Arguments( "mergeWithNeighbour(...)", vertex, neighbour, fromRank, fineGridX, fineGridH, level );
  // do nothing
  logTraceOut( "mergeWithNeighbour(...)" );
}

void exahype::mappings::InitialCondition::prepareSendToNeighbour(
    exahype::Vertex&  vertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  logTraceInWith3Arguments( "prepareSendToNeighbour(...)", vertex, toRank, level );
  // do nothing
  logTraceOut( "prepareSendToNeighbour(...)" );
}

void exahype::mappings::InitialCondition::prepareCopyToRemoteNode(
    exahype::Vertex&  localVertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  logTraceInWith5Arguments( "prepareCopyToRemoteNode(...)", localVertex, toRank, x, h, level );
  // do nothing
  logTraceOut( "prepareCopyToRemoteNode(...)" );
}

void exahype::mappings::InitialCondition::prepareCopyToRemoteNode(
    exahype::Cell&  localCell,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
    int                                           level
) {
  logTraceInWith5Arguments( "prepareCopyToRemoteNode(...)", localCell, toRank, cellCentre, cellSize, level );
  // do nothing
  logTraceOut( "prepareCopyToRemoteNode(...)" );
}

void exahype::mappings::InitialCondition::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex&  localVertex,
    const exahype::Vertex&  masterOrWorkerVertex,
    int                                       fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&  x,
    const tarch::la::Vector<DIMENSIONS,double>&  h,
    int                                       level
) {
  logTraceInWith6Arguments( "mergeWithRemoteDataDueToForkOrJoin(...)", localVertex, masterOrWorkerVertex, fromRank, x, h, level );
  // do nothing
  logTraceOut( "mergeWithRemoteDataDueToForkOrJoin(...)" );
}

void exahype::mappings::InitialCondition::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell&  localCell,
    const exahype::Cell&  masterOrWorkerCell,
    int                                       fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                       level
) {
  logTraceInWith3Arguments( "mergeWithRemoteDataDueToForkOrJoin(...)", localCell, masterOrWorkerCell, fromRank );
  // do nothing
  logTraceOut( "mergeWithRemoteDataDueToForkOrJoin(...)" );
}

bool exahype::mappings::InitialCondition::prepareSendToWorker(
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
  // do nothing
  logTraceOutWith1Argument( "prepareSendToWorker(...)", true );
  return true;
}

void exahype::mappings::InitialCondition::prepareSendToMaster(
    exahype::Cell&                       localCell,
    exahype::Vertex *                    vertices,
    const peano::grid::VertexEnumerator&       verticesEnumerator,
    const exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
    const exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
  logTraceInWith2Arguments( "prepareSendToMaster(...)", localCell, verticesEnumerator.toString() );
  // do nothing
  logTraceOut( "prepareSendToMaster(...)" );
}


void exahype::mappings::InitialCondition::mergeWithMaster(
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
  logTraceIn( "mergeWithMaster(...)" );
  // do nothing
  logTraceOut( "mergeWithMaster(...)" );
}


void exahype::mappings::InitialCondition::receiveDataFromMaster(
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
  logTraceIn( "receiveDataFromMaster(...)" );
  // do nothing
  logTraceOut( "receiveDataFromMaster(...)" );
}


void exahype::mappings::InitialCondition::mergeWithWorker(
    exahype::Cell&           localCell,
    const exahype::Cell&     receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                          level
) {
  logTraceInWith2Arguments( "mergeWithWorker(...)", localCell.toString(), receivedMasterCell.toString() );
  // do nothing
  logTraceOutWith1Argument( "mergeWithWorker(...)", localCell.toString() );
}


void exahype::mappings::InitialCondition::mergeWithWorker(
    exahype::Vertex&        localVertex,
    const exahype::Vertex&  receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  logTraceInWith2Arguments( "mergeWithWorker(...)", localVertex.toString(), receivedMasterVertex.toString() );
  // do nothing
  logTraceOutWith1Argument( "mergeWithWorker(...)", localVertex.toString() );
}
#endif

void exahype::mappings::InitialCondition::touchVertexFirstTime(
    exahype::Vertex&               fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "touchVertexFirstTime(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // do nothing
  logTraceOutWith1Argument( "touchVertexFirstTime(...)", fineGridVertex );
}


void exahype::mappings::InitialCondition::touchVertexLastTime(
    exahype::Vertex&         fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "touchVertexLastTime(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );
  // do nothing
  logTraceOutWith1Argument( "touchVertexLastTime(...)", fineGridVertex );
}

void exahype::mappings::InitialCondition::enterCell(
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "enterCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );

  // ! Begin of code for the DG method.
  if (!fineGridCell.isRefined()) {
    records::CellDescription& cellDescription =
             CellDescriptionHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[0];

    const tarch::la::Vector<DIMENSIONS,double> center = fineGridVerticesEnumerator.getCellCenter();  // the center of the cell
    const double dx = fineGridVerticesEnumerator.getCellSize()(0);
    const double dy = fineGridVerticesEnumerator.getCellSize()(1);

    const double dxPatch = dx/ (double) EXAHYPE_PATCH_SIZE_X;
    const double dyPatch = dy/ (double) EXAHYPE_PATCH_SIZE_Y;

    const int basisSize = EXAHYPE_ORDER+1;
    const int nvar      = EXAHYPE_NVARS;

    // helper variables
    double x,y;
    double* value = (double*) std::malloc(nvar * sizeof(double));

    for (int i=1; i<EXAHYPE_PATCH_SIZE_X+1; i++) { // loop over patches
      for (int j=1; j<EXAHYPE_PATCH_SIZE_Y+1; j++) {
        const int patchIndex = i + (EXAHYPE_PATCH_SIZE_X+2) * j;

        double* luh  = &(DataHeap::getInstance().getData(cellDescription.getSolution(patchIndex))[0]._persistentRecords._u);
        double* lduh = &(DataHeap::getInstance().getData(cellDescription.getUpdate(patchIndex))  [0]._persistentRecords._u);

        for (int ii=0; ii<basisSize; ii++) { // loop over dof
          for (int jj=0; jj<basisSize; jj++) {
            // location and index of nodal degrees of freedom
            const int nodeIndex = ii + basisSize * jj;

            const double qr = exahype::quad::gaussLegendreNodes[basisSize-1][ii];
            const double qs = exahype::quad::gaussLegendreNodes[basisSize-1][jj];
            exahype::geometry::mapping2d(center(0),center(1),dx,dy,dxPatch,dyPatch,i,j,qr,qs,&x,&y);

            // read initial condition
            exahype::problem::PDEInitialValue2d(x,y,nvar,value);

            // set the DoF
            const int dofStartIndex  = nodeIndex * nvar;

            for (int ivar=0; ivar < nvar; ivar++) {
              luh[dofStartIndex+ivar] = value[ivar];

              lduh[dofStartIndex+ivar] = 0.;
            }
          }
        }
      }
    }

    // clean up
    std::free(value);
  }

  // ! End of code for the DG method.

  logTraceOutWith1Argument( "enterCell(...)", fineGridCell );
}


void exahype::mappings::InitialCondition::leaveCell(
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


void exahype::mappings::InitialCondition::beginIteration(
    exahype::State&  solverState
) {
  logTraceInWith1Argument( "beginIteration(State)", solverState );
  // do nothing
  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


void exahype::mappings::InitialCondition::endIteration(
    exahype::State&  solverState
) {
  logTraceInWith1Argument( "endIteration(State)", solverState );
  // do nothing
  logTraceOutWith1Argument( "endIteration(State)", solverState);
}



void exahype::mappings::InitialCondition::descend(
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


void exahype::mappings::InitialCondition::ascend(
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
