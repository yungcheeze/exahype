#include "EulerFlow3d/mappings/ExchangeFaceData.h"

#include "EulerFlow3d/dg/Constants.h"

#include "EulerFlow3d/geometry/Mapping.h"
#include "EulerFlow3d/dg/DGHelpers.h"
#include "EulerFlow3d/math/quad/Gausslegendre.h"

#include "EulerFlow3d/multiscalelinkedcell/HangingVertexBookkeeper.h"
#include "EulerFlow3d/VertexOperations.h"

#include "EulerFlow3d/problem/Problem.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification   exahype::mappings::ExchangeFaceData::communicationSpecification() {
  return peano::CommunicationSpecification(peano::CommunicationSpecification::SendDataAndStateBeforeFirstTouchVertexFirstTime,peano::CommunicationSpecification::SendDataAndStateAfterLastTouchVertexLastTime,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::ExchangeFaceData::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::ExchangeFaceData::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::ExchangeFaceData::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::ExchangeFaceData::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::ExchangeFaceData::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidCoarseGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::ExchangeFaceData::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidCoarseGridRaces);
}


tarch::logging::Log                exahype::mappings::ExchangeFaceData::_log( "exahype::mappings::ExchangeFaceData" ); 


exahype::mappings::ExchangeFaceData::ExchangeFaceData() {
  logTraceIn( "ExchangeFaceData()" );
  // do nothing
  logTraceOut( "ExchangeFaceData()" );
}


exahype::mappings::ExchangeFaceData::~ExchangeFaceData() {
  logTraceIn( "~ExchangeFaceData()" );
  // do nothing
  logTraceOut( "~ExchangeFaceData()" );
}


#if defined(SharedMemoryParallelisation)
exahype::mappings::ExchangeFaceData::ExchangeFaceData(const ExchangeFaceData&  masterThread) {
  logTraceIn( "ExchangeFaceData(ExchangeFaceData)" );
  // do nothing
  logTraceOut( "ExchangeFaceData(ExchangeFaceData)" );
}


void exahype::mappings::ExchangeFaceData::mergeWithWorkerThread(const ExchangeFaceData& workerThread) {
  logTraceIn( "mergeWithWorkerThread(ExchangeFaceData)" );
  // do nothing
  logTraceOut( "mergeWithWorkerThread(ExchangeFaceData)" );
}
#endif


void exahype::mappings::ExchangeFaceData::createHangingVertex(
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


void exahype::mappings::ExchangeFaceData::destroyHangingVertex(
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


void exahype::mappings::ExchangeFaceData::createInnerVertex(
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


void exahype::mappings::ExchangeFaceData::createBoundaryVertex(
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


void exahype::mappings::ExchangeFaceData::destroyVertex(
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



void exahype::mappings::ExchangeFaceData::createCell(
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


void exahype::mappings::ExchangeFaceData::destroyCell(
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
void exahype::mappings::ExchangeFaceData::mergeWithNeighbour(
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

void exahype::mappings::ExchangeFaceData::prepareSendToNeighbour(
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

void exahype::mappings::ExchangeFaceData::prepareCopyToRemoteNode(
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

void exahype::mappings::ExchangeFaceData::prepareCopyToRemoteNode(
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

void exahype::mappings::ExchangeFaceData::mergeWithRemoteDataDueToForkOrJoin(
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

void exahype::mappings::ExchangeFaceData::mergeWithRemoteDataDueToForkOrJoin(
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

bool exahype::mappings::ExchangeFaceData::prepareSendToWorker(
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

void exahype::mappings::ExchangeFaceData::prepareSendToMaster(
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


void exahype::mappings::ExchangeFaceData::mergeWithMaster(
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


void exahype::mappings::ExchangeFaceData::receiveDataFromMaster(
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


void exahype::mappings::ExchangeFaceData::mergeWithWorker(
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


void exahype::mappings::ExchangeFaceData::mergeWithWorker(
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

void exahype::mappings::ExchangeFaceData::touchVertexFirstTime(
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


void exahype::mappings::ExchangeFaceData::touchVertexLastTime(
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

void addSelfContribution(exahype::records::CellDescription& cellDescription,const double selfFlux) {
  exahype::DataHeap::getInstance().getData(cellDescription.getUpdate())[0]._persistentRecords._u +=
      selfFlux * exahype::DataHeap::getInstance().getData(cellDescription.getSolution())[0]._persistentRecords._u;
}

void setNeighbourContribution(exahype::records::CellDescription& cellDescription_c,
                              exahype::records::CellDescription& cellDescription_b,
                              const int face_c,const double flux_cb) {
  exahype::DataHeap::getInstance().getData(cellDescription_c.getFluctuation(face_c))[0]._persistentRecords._u =
      flux_cb * exahype::DataHeap::getInstance().getData(cellDescription_b.getSolution())[0]._persistentRecords._u;
}

void exahype::mappings::ExchangeFaceData::enterCell(
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "enterCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );

  if (!fineGridCell.isRefined()) {
    initialiseGhostLayerOfPatch(fineGridCell,fineGridVertices,fineGridVerticesEnumerator);


    const int n_gauss_points=2;

    // cell geometry
    const double h                                    = fineGridVerticesEnumerator.getCellSize()(0);
    const double ds                                   = 0.5*h;                                       // the length of the path element; constant on the face
    const tarch::la::Vector<DIMENSIONS,double> center = fineGridVerticesEnumerator.getCellCenter();  // the center of the cell
    double x                                          = 0.;                                          // physical x coordinate
    double y                                          = 0.;                                          // physical y coordinate

    // helper variables
    int    face_b  =0;
    double fluxL   =0.;
    double fluxR   =0.;
    double flux_cc =0.;
    double flux_cb =0.;
    double flux_bb =0.;
    double flux_bc =0.;

    records::CellDescription& cellDescriptionForPde_c =
        CellDescriptionHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[0];

    tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>  cellDescriptionsOfNeighbors;
    const tarch::la::Vector<THREE_POWER_D,int> cellDescriptionsOfAllNeighbours =
        multiscalelinkedcell::getIndicesAroundCell(
            exahype::VertexOperations::readCellDescriptionsIndex(fineGridVerticesEnumerator,fineGridVertices));

    // Store neighbor flux in LEFT (0,r=-1),RIGHT (1,r=+1),FRONT (2,s=-1),BACK (3,s=+1),BOTTOM (4,t=-1),TOP (5,t=+1)
    // manner, where r,s,t refer to the reference coordinates.
    cellDescriptionsOfNeighbors[EXAHYPE_FACE_LEFT  ] = cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_LEFT  ];
    cellDescriptionsOfNeighbors[EXAHYPE_FACE_RIGHT ] = cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_RIGHT ];
    cellDescriptionsOfNeighbors[EXAHYPE_FACE_FRONT ] = cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_FRONT ];
    cellDescriptionsOfNeighbors[EXAHYPE_FACE_BACK  ] = cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_BACK  ];


    for (int face_c=0; face_c < DIMENSIONS_TIMES_TWO; face_c++) {
      if (cellDescriptionsOfNeighbors[face_c] > multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) { // todo iterate over right, back, top
        records::CellDescription& cellDescriptionForPde_b =
            CellDescriptionHeap::getInstance().getData(cellDescriptionsOfNeighbors[face_c])[0];
        exahype::dg::GetNeighborFace(face_c,&face_b);

        // non-optimized quadrature loop;
        double qr;
        double qs;
        for (int iq=0; iq < n_gauss_points; iq++) {
          exahype::dg::GetFaceQr(n_gauss_points,iq,face_c,&qr);
          exahype::dg::GetFaceQs(n_gauss_points,iq,face_c,&qs);
          exahype::geometry::mapping2d(center(0),center(1),h,qr,qs,&x,&y);

          // solve the Riemann problem
          exahype::problem::DGRiemannSolver(x,y,exahype::dg::normal[0][face_c],exahype::dg::normal[1][face_c],&fluxL,&fluxR);
          // flux_cc -= exahype::quad::gaussLegendreWeights[n_gauss_points-1][iq] * 0 * ds;
          // flux_cb += exahype::quad::gaussLegendreWeights[n_gauss_points-1][iq] * 0 * ds;
          // flux_bb -= exahype::quad::gaussLegendreWeights[n_gauss_points-1][iq] * 0 * ds;
          // flux_bc += exahype::quad::gaussLegendreWeights[n_gauss_points-1][iq] * 0 * ds;
        }
        addSelfContribution(cellDescriptionForPde_c,flux_cc);
        setNeighbourContribution(cellDescriptionForPde_b,cellDescriptionForPde_c,face_b,flux_bc);
        addSelfContribution(cellDescriptionForPde_b,flux_bb);
        setNeighbourContribution(cellDescriptionForPde_c,cellDescriptionForPde_b,face_c,flux_cb);
      }
    }
  }
  logTraceOutWith1Argument( "enterCell(...)", fineGridCell );
}


void exahype::mappings::ExchangeFaceData::leaveCell(
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


void exahype::mappings::ExchangeFaceData::beginIteration(
    exahype::State&  solverState
) {
  logTraceInWith1Argument( "beginIteration(State)", solverState );
  // do nothing
  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


void exahype::mappings::ExchangeFaceData::endIteration(
    exahype::State&  solverState
) {
  logTraceInWith1Argument( "endIteration(State)", solverState );
  // do nothing
  logTraceOutWith1Argument( "endIteration(State)", solverState);
}



void exahype::mappings::ExchangeFaceData::descend(
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


void exahype::mappings::ExchangeFaceData::ascend(
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

// Begin of code for ADERDG method
void exahype::mappings::ExchangeFaceData::initialiseGhostLayerOfPatch(
    exahype::Cell& fineGridCell,
    exahype::Vertex * const  fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator
) {
  int basisSize         = exahype::order[0]+1;
  int numberOfDofOnFace = exahype::numberOfVariables[0] * tarch::la::aPowI(DIMENSIONS-1,basisSize);

  records::CellDescription& dataSelf =
      CellDescriptionHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[0];

  // Read in neighbor information
  tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>  cellDescriptionsOfNeighbours;

  const tarch::la::Vector<THREE_POWER_D,int> cellDescriptionsOfAllNeighbours =
      multiscalelinkedcell::getIndicesAroundCell(
          exahype::VertexOperations::readCellDescriptionsIndex(fineGridVerticesEnumerator,fineGridVertices));

  int indexGhostSelf = 0;
  int indexNeighbor  = 0;

  for (int i=0; i<EXAHYPE_PATCH_SIZE_X+2; i++) {
    // front
    indexGhostSelf = i + (EXAHYPE_PATCH_SIZE_X+2) * 0;
    indexNeighbor  = i + (EXAHYPE_PATCH_SIZE_X+2) * EXAHYPE_PATCH_SIZE_Y;

    records::CellDescription& dataNeighbourFront =
        CellDescriptionHeap::getInstance().getData(cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_FRONT])[0];
    for (int dof=0; dof<numberOfDofOnFace; dof++) {
      DataHeap::getInstance().getData(dataSelf.getExtrapolatedPredictor(indexGhostSelf))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_BACK)*numberOfDofOnFace + dof]._persistentRecords._u =
          DataHeap::getInstance().getData(dataNeighbourFront.getExtrapolatedPredictor(indexNeighbor))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_BACK)*numberOfDofOnFace + dof]._persistentRecords._u;
      DataHeap::getInstance().getData(dataSelf.getFluctuation(indexGhostSelf))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_BACK)*numberOfDofOnFace + dof]._persistentRecords._u =
                DataHeap::getInstance().getData(dataNeighbourFront.getFluctuation(indexNeighbor))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_BACK) + dof]._persistentRecords._u;
    }


    // back
    indexGhostSelf = i + (EXAHYPE_PATCH_SIZE_X+2) * (EXAHYPE_PATCH_SIZE_Y+1); // ghost
    indexNeighbor  = i + (EXAHYPE_PATCH_SIZE_X+2) * 1;                        // non-ghost

    records::CellDescription& dataNeighbourBack =
        CellDescriptionHeap::getInstance().getData(cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_BACK])[0];
    for (int dof=0; dof<numberOfDofOnFace; dof++) {
      DataHeap::getInstance().getData(dataSelf.getExtrapolatedPredictor(indexGhostSelf))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_FRONT)*numberOfDofOnFace + dof]._persistentRecords._u =
          DataHeap::getInstance().getData(dataNeighbourBack.getExtrapolatedPredictor(indexNeighbor))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_FRONT) + dof]._persistentRecords._u;
      DataHeap::getInstance().getData(dataSelf.getFluctuation(indexGhostSelf))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_FRONT)*numberOfDofOnFace + dof]._persistentRecords._u =
          DataHeap::getInstance().getData(dataNeighbourBack.getFluctuation(indexNeighbor))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_FRONT) + dof]._persistentRecords._u;
    }
  }

  for (int j=0; j<EXAHYPE_PATCH_SIZE_Y+2; j++) {
    // left
    indexGhostSelf = 0 +                    (EXAHYPE_PATCH_SIZE_X+2) * j;
    indexNeighbor  = EXAHYPE_PATCH_SIZE_X + (EXAHYPE_PATCH_SIZE_X+2) * j;

    records::CellDescription& dataNeighbourLeft =
        CellDescriptionHeap::getInstance().getData(cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_LEFT])[0];
    for (int dof=0; dof<numberOfDofOnFace; dof++) {
      DataHeap::getInstance().getData(dataSelf.getExtrapolatedPredictor(indexGhostSelf))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_RIGHT)*numberOfDofOnFace + dof]._persistentRecords._u =
          DataHeap::getInstance().getData(dataNeighbourLeft.getExtrapolatedPredictor(indexNeighbor))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_RIGHT) + dof]._persistentRecords._u;
      DataHeap::getInstance().getData(dataSelf.getFluctuation(indexGhostSelf))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_RIGHT)*numberOfDofOnFace + dof]._persistentRecords._u =
          DataHeap::getInstance().getData(dataNeighbourLeft.getFluctuation(indexNeighbor))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_RIGHT) + dof]._persistentRecords._u;
    }

    // right
    indexGhostSelf = (EXAHYPE_PATCH_SIZE_X+1) + (EXAHYPE_PATCH_SIZE_X+2) * j;
    indexNeighbor  = 1                        + (EXAHYPE_PATCH_SIZE_X+2) * j;

    records::CellDescription& dataNeighbourRight =
        CellDescriptionHeap::getInstance().getData(cellDescriptionsOfAllNeighbours[PEANO_2D_NEIGHBOUR_RIGHT])[0];
    for (int dof=0; dof<numberOfDofOnFace; dof++) {
      DataHeap::getInstance().getData(dataSelf.getExtrapolatedPredictor(indexGhostSelf))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_LEFT)*numberOfDofOnFace + dof]._persistentRecords._u =
          DataHeap::getInstance().getData(dataNeighbourRight.getExtrapolatedPredictor(indexNeighbor))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_LEFT) + dof]._persistentRecords._u;
      DataHeap::getInstance().getData(dataSelf.getFluctuation(indexGhostSelf))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_LEFT)*numberOfDofOnFace + dof]._persistentRecords._u =
          DataHeap::getInstance().getData(dataNeighbourRight.getFluctuation(indexNeighbor))[(DIMENSIONS_TIMES_TWO+EXAHYPE_FACE_LEFT) + dof]._persistentRecords._u;
    }
  }
}
// End of code for ADERDG method

