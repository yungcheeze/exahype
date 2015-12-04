#include "EulerFlow3d/mappings/SpaceTimePredictor.h"

#include "EulerFlow3d/Constants.h"

#include "EulerFlow3d/math/quad/Gausslegendre.h"

#include "EulerFlow3d/geometry/Mapping.h"

#include "EulerFlow3d/problem/Problem.h"

#include "EulerFlow3d/dg/Constants.h"
#include "EulerFlow3d/dg/DGMatrices.h"

#include "stdlib.h"

#include "string.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification   exahype::mappings::SpaceTimePredictor::communicationSpecification() {
  return peano::CommunicationSpecification(peano::CommunicationSpecification::SendDataAndStateBeforeFirstTouchVertexFirstTime,peano::CommunicationSpecification::SendDataAndStateAfterLastTouchVertexLastTime,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::SpaceTimePredictor::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::SpaceTimePredictor::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::SpaceTimePredictor::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::SpaceTimePredictor::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::SpaceTimePredictor::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidCoarseGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::SpaceTimePredictor::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidCoarseGridRaces);
}


tarch::logging::Log                exahype::mappings::SpaceTimePredictor::_log( "exahype::mappings::SpaceTimePredictor" ); 


exahype::mappings::SpaceTimePredictor::SpaceTimePredictor() {
  logTraceIn( "SpaceTimePredictor()" );
  // do nothing
  logTraceOut( "SpaceTimePredictor()" );
}


exahype::mappings::SpaceTimePredictor::~SpaceTimePredictor() {
  logTraceIn( "~SpaceTimePredictor()" );
  // do nothing
  logTraceOut( "~SpaceTimePredictor()" );
}


#if defined(SharedMemoryParallelisation)
exahype::mappings::SpaceTimePredictor::SpaceTimePredictor(const SpaceTimePredictor&  masterThread) {
  logTraceIn( "SpaceTimePredictor(SpaceTimePredictor)" );
  // do nothing
  logTraceOut( "SpaceTimePredictor(SpaceTimePredictor)" );
}


void exahype::mappings::SpaceTimePredictor::mergeWithWorkerThread(const SpaceTimePredictor& workerThread) {
  logTraceIn( "mergeWithWorkerThread(SpaceTimePredictor)" );
  // do nothing
  logTraceOut( "mergeWithWorkerThread(SpaceTimePredictor)" );
}
#endif


void exahype::mappings::SpaceTimePredictor::createHangingVertex(
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


void exahype::mappings::SpaceTimePredictor::destroyHangingVertex(
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


void exahype::mappings::SpaceTimePredictor::createInnerVertex(
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


void exahype::mappings::SpaceTimePredictor::createBoundaryVertex(
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


void exahype::mappings::SpaceTimePredictor::destroyVertex(
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


void exahype::mappings::SpaceTimePredictor::createCell(
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


void exahype::mappings::SpaceTimePredictor::destroyCell(
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
void exahype::mappings::SpaceTimePredictor::mergeWithNeighbour(
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

void exahype::mappings::SpaceTimePredictor::prepareSendToNeighbour(
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

void exahype::mappings::SpaceTimePredictor::prepareCopyToRemoteNode(
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

void exahype::mappings::SpaceTimePredictor::prepareCopyToRemoteNode(
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

void exahype::mappings::SpaceTimePredictor::mergeWithRemoteDataDueToForkOrJoin(
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

void exahype::mappings::SpaceTimePredictor::mergeWithRemoteDataDueToForkOrJoin(
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

bool exahype::mappings::SpaceTimePredictor::prepareSendToWorker(
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

void exahype::mappings::SpaceTimePredictor::prepareSendToMaster(
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


void exahype::mappings::SpaceTimePredictor::mergeWithMaster(
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


void exahype::mappings::SpaceTimePredictor::receiveDataFromMaster(
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


void exahype::mappings::SpaceTimePredictor::mergeWithWorker(
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


void exahype::mappings::SpaceTimePredictor::mergeWithWorker(
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

void exahype::mappings::SpaceTimePredictor::touchVertexFirstTime(
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


void exahype::mappings::SpaceTimePredictor::touchVertexLastTime(
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


void exahype::mappings::SpaceTimePredictor::computePredictor(
    double *lQi,
    double* lFi,
    double* luh,
    double* lQhi,
    double* lFhi,
    double* lQhbnd,
    double* lFhbnd,
    const tarch::la::Vector<DIMENSIONS,double> center,
    const double dxPatch,const double dyPatch,
    const double dt,
    const int patchIndex,
    const int nvar,
    const int order) {
  // helper variables
  const int basisSize = order+1;
  int numberOfSpaceTimeDof  = nvar * tarch::la::aPowI(DIMENSIONS+1,basisSize);

  double* rhs0 = (double*) std::malloc(numberOfSpaceTimeDof * sizeof(double)); // todo remove all mallocs, nvar is known at compile time and the same for optimization problems
  double* rhs  = (double*) std::malloc(numberOfSpaceTimeDof * sizeof(double)); // todo remove all mallocs, nvar is known at compile time and the same for optimization problems

  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    for (int jj=0; jj<basisSize; jj++) {
      for (int ll=0; ll<basisSize; ll++) { // loop over dof
        // location and index of nodal degrees of freedom
        const int nodeIndex          = ii + basisSize * jj;
        const int spaceTimeNodeIndex = ii + basisSize * jj  + basisSize * basisSize * ll;

        const int dofStartIndex           = nodeIndex * nvar;
        const int spaceTimeDofStartIndex  = spaceTimeNodeIndex * nvar;

        for (int ivar=0; ivar < nvar; ivar++) {
          // Trivial initial guess (can be significantly improved)
          lQi[spaceTimeDofStartIndex+ivar] = luh[dofStartIndex+ivar];

          // Compute the contribution of the initial condition uh to the time update. I prefer to compute it once
          // and store it in rhs0, but if you think it is faster, you can also recompute this contribution
          // inside the Picard loop (DO iter = 1, N+1)
          rhs0[spaceTimeDofStartIndex+ivar] =
              quad::gaussLegendreWeights[basisSize-1][ii] *
              quad::gaussLegendreWeights[basisSize-1][jj] *
              dg::F0[ll] *
              luh[dofStartIndex+ivar];
        }
      }
    }
  }
  // Above seems to work!

  double*    Q;
  double*    f;
  double*    g;

  double* tmp  = (double*) std::malloc(nvar * basisSize * sizeof(double));
  double* dqdt = (double*) std::malloc(nvar * basisSize * sizeof(double)); // todo remove all mallocs, nvar is known at compile time and the same for optimization problems

  // Discrete Picard iterations. This set of nested loops should (theoretically) be a dream for vectorization, since they are rather independent...
  for (int iter=1; iter < basisSize+1; iter++) {

    // Compute the fluxes (once these fluxes are available, the subsequent operations are independent from each other)
    for (int ll=0; ll<basisSize; ll++) { // loop over dof (time)
      for (int ii=0; ii<basisSize; ii++) { // loop over dof
        for (int jj=0; jj<basisSize; jj++) {
          const int nodeIndex          = ii + basisSize * jj;
          const int spaceTimeNodeIndex = nodeIndex  + basisSize * basisSize * ll;

          const int spaceTimeDofStartIndex     = spaceTimeNodeIndex * nvar;
          const int spaceTimeFluxDofStartIndex = spaceTimeDofStartIndex * DIMENSIONS;

          Q = &lQi [spaceTimeDofStartIndex        ];
          f = &lFi[spaceTimeFluxDofStartIndex     ];
          g = &lFi[spaceTimeFluxDofStartIndex+nvar];
          exahype::problem::PDEFlux(Q,nvar,f,g);
        }
      }
      // Above seems okay!

      // x direction (independent from the y and z derivatives)
      // Kxi : basisSize * basisSize
      // lFh : nvar * basisSize

      // Compute the "derivatives" (contributions of the stiffness matrix)
      // x direction (independent from the y and z derivatives)
      for (int ii=0; ii<basisSize; ii++) { // loop over dof
        for (int jj=0; jj<basisSize; jj++) {
          const int nodeIndex              = ii + basisSize * jj;
          const int spaceTimeNodeIndex     = nodeIndex  + basisSize * basisSize * ll;
          const int spaceTimeDofStartIndex = spaceTimeNodeIndex * nvar;

          double weight = quad::gaussLegendreWeights[basisSize-1][ll] *
              quad::gaussLegendreWeights[basisSize-1][jj];

          // COMPUTE SPATIAL DERIVATIVES FOR TESTING PURPOSES
          //          for (int ll=0; ll < basisSize; ll++) { // set tmp = 0
          //            for(int ivar=0; ivar < nvar; ivar++) {
          //              tmp[ivar + nvar*ll] = 0.;
          //            }
          //          }
          //
          //          for(int mm=0; mm < basisSize; mm++) {
          //            const int mmNodeIndex                  = mm + basisSize * jj;
          //            const int mmSpaceTimeNodeIndex         = mmNodeIndex  + basisSize * basisSize * ll;
          //            const int mmSpaceTimeDofStartIndex     = mmSpaceTimeNodeIndex * nvar;
          //            const int mmSpaceTimeFluxDofStartIndex = mmSpaceTimeDofStartIndex * DIMENSIONS;
          //
          //            Q = &(lQi [mmSpaceTimeDofStartIndex]);
          //            for(int ivar=0; ivar < nvar; ivar++) {
          //              tmp[ivar] += 1./dxPatch * dg::dudx[ii][mm] * Q[ivar];
          //            }
          //          }
          for(int ivar=0; ivar < nvar; ivar++) {
            rhs[spaceTimeDofStartIndex+ivar] = rhs0[spaceTimeDofStartIndex+ivar];
          }

          for(int mm=0; mm < basisSize; mm++) {
            const int mmNodeIndex                  = mm + basisSize * jj;
            const int mmSpaceTimeNodeIndex         = mmNodeIndex  + basisSize * basisSize * ll;
            const int mmSpaceTimeDofStartIndex     = mmSpaceTimeNodeIndex * nvar;
            const int mmSpaceTimeFluxDofStartIndex = mmSpaceTimeDofStartIndex * DIMENSIONS;

            f = &lFi[mmSpaceTimeFluxDofStartIndex];

            for(int ivar=0; ivar < nvar; ivar++) {
              rhs[spaceTimeDofStartIndex+ivar]
                  -= weight * dt/dxPatch * dg::Kxi[mm][ii] * f[ivar];
            }
          }
        }
      }
      // Above seems okay!

      // Compute the "derivatives" (contributions of the stiffness matrix)
      // y direction (independent from the x and z derivatives)
      for (int ii=0; ii<basisSize; ii++) { // loop over dof
        for (int jj=0; jj<basisSize; jj++) {
          const int nodeIndex              = ii + basisSize * jj;
          const int spaceTimeNodeIndex     = nodeIndex  + basisSize * basisSize * ll;
          const int spaceTimeDofStartIndex = spaceTimeNodeIndex * nvar;

          double weight = quad::gaussLegendreWeights[basisSize-1][ll] *
              quad::gaussLegendreWeights[basisSize-1][ii];

          for(int mm=0; mm < basisSize; mm++) {
            const int mmNodeIndex                  = ii + basisSize * mm;
            const int mmSpaceTimeNodeIndex         = mmNodeIndex  + basisSize * basisSize * ll;
            const int mmSpaceTimeDofStartIndex     = mmSpaceTimeNodeIndex * nvar;
            const int mmSpaceTimeFluxDofStartIndex = mmSpaceTimeDofStartIndex * DIMENSIONS;

            g = &lFi[mmSpaceTimeFluxDofStartIndex+nvar];

            for(int ivar=0; ivar < nvar; ivar++) {
              rhs[spaceTimeDofStartIndex+ivar]
                  -= weight * dt/dyPatch * dg::Kxi[mm][jj] * g[ivar];
            }
          }
        }
      }
    } // end of time dof loop

    // Above seems okay!

    for (int ii=0; ii<basisSize; ii++) {  // loop over dof
      for (int jj=0; jj<basisSize; jj++) {
        const int nodeIndex = ii + basisSize * jj;

        double iWeight = 1./(quad::gaussLegendreWeights[basisSize-1][ii] * quad::gaussLegendreWeights[basisSize-1][jj]);

        for (int ll=0; ll < basisSize; ll++) { // set tmp = 0
          for(int ivar=0; ivar < nvar; ivar++) {
            tmp[ivar + nvar*ll] = 0.;
          }
        }

        for (int ll=0; ll<basisSize; ll++) { // loop over dof

          for(int nn=0; nn < basisSize; nn++) {
            const int nnSpaceTimeNodeIndex     = nodeIndex  + basisSize * basisSize * nn;
            const int nnSpaceTimeDofStartIndex = nnSpaceTimeNodeIndex * nvar;

            for(int ivar=0; ivar < nvar; ivar++) {
              tmp[ivar + nvar*ll] += iWeight * dg::iK1[ll][nn] * rhs[nnSpaceTimeDofStartIndex+ivar];
            }
          }
        }

        for (int ll=0; ll<basisSize; ll++) { // loop over dof
          const int spaceTimeNodeIndex     = nodeIndex  + basisSize * basisSize * ll;
          const int spaceTimeDofStartIndex = spaceTimeNodeIndex * nvar;

          for(int ivar=0; ivar < nvar; ivar++) {
            lQi[spaceTimeDofStartIndex+ivar] = tmp[ivar + nvar*ll];
          }
        }

        // dqdt
        for (int ll=0; ll < basisSize; ll++) { // set tmp = 0
          for(int ivar=0; ivar < nvar; ivar++) {
            dqdt[ivar + nvar*ll] = 0.;
          }
        }

        for (int ll=0; ll<basisSize; ll++) { // loop over dof

          for(int ivar=0; ivar < nvar; ivar++) {

            for(int nn=0; nn < basisSize; nn++) {
              const int nnSpaceTimeNodeIndex         = nodeIndex  + basisSize * basisSize * nn;
              const int nnSpaceTimeDofStartIndex     = nnSpaceTimeNodeIndex * nvar;

              dqdt[ivar + nvar*ll] += 1./dt * dg::dudx[ll][nn] *
                  lQi[nnSpaceTimeDofStartIndex+ivar];
            }
          }
        }
      }
    }
  } // end of Picard iteration

  /////////////////////////////////////////////////
  // Post processing of the predictor:
  // Immediately compute the time-averaged space-time polynomials
  /////////////////////////////////////////////////
  int numberOfDof      = nvar * tarch::la::aPowI(DIMENSIONS,basisSize);
  int numberOfFluxDof  = numberOfDof * DIMENSIONS;

  memset((double *) lQhi,0,sizeof(double) * numberOfDof);
  memset((double *) lFhi,0,sizeof(double) * numberOfFluxDof);

  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    for (int jj=0; jj<basisSize; jj++) {
      const int nodeIndex     = ii + basisSize * jj;
      const int dofStartIndex = nodeIndex * nvar;
      const int fluxDofStartIndex = DIMENSIONS * dofStartIndex;

      for (int ll=0; ll<basisSize; ll++) { // loop over dof
        const int spaceTimeNodeIndex         = nodeIndex  + basisSize * basisSize * ll;
        const int spaceTimeDofStartIndex     = spaceTimeNodeIndex * nvar;
        const int spaceTimeFluxDofStartIndex = spaceTimeDofStartIndex * DIMENSIONS;

        Q = &lQi[spaceTimeDofStartIndex];

        f = &lFi[spaceTimeFluxDofStartIndex     ];
        g = &lFi[spaceTimeFluxDofStartIndex+nvar];

        double weight = quad::gaussLegendreWeights[basisSize-1][ll];

        double * temp = &(lQhi[dofStartIndex]);
        for(int ivar=0; ivar < nvar; ivar++) {
          lQhi[dofStartIndex+ivar] += weight * Q[ivar];

          lFhi[fluxDofStartIndex+ivar     ] += weight * f[ivar];
          lFhi[fluxDofStartIndex+nvar+ivar] += weight * g[ivar];
        }
      }
    }
  }

  /////////////////////////////////////////////////
  // Compute the bounday-extrapolated values for Q and F*n
  /////////////////////////////////////////////////
  int numberOfFaceDof = nvar * tarch::la::aPowI(DIMENSIONS-1,basisSize);

  memset((double *) &lQhbnd[0],0,sizeof(double) * numberOfFaceDof * DIMENSIONS_TIMES_TWO);
  memset((double *) &lFhbnd[0],0,sizeof(double) * numberOfFaceDof * DIMENSIONS_TIMES_TWO);

  // x-direction: face 0 (left) and face 1 (right)
  for (int jj=0; jj<basisSize; jj++) {
    const int nodeIndex      = jj;
    const int dofStartIndexL = EXAHYPE_FACE_LEFT  * numberOfFaceDof + nodeIndex * nvar;
    const int dofStartIndexR = EXAHYPE_FACE_RIGHT * numberOfFaceDof + nodeIndex * nvar;

    double * tempQL = &lQhbnd[dofStartIndexL];
    double * tempQR = &lQhbnd[dofStartIndexR];
    double * tempL  = &lFhbnd[dofStartIndexL];
    double * tempR  = &lFhbnd[dofStartIndexR];

    for (int mm=0; mm<basisSize; mm++) { // loop over dof
      const int mmNodeIndex         = mm  + basisSize * jj;
      const int mmDofStartIndex     = mmNodeIndex * nvar;
      const int mmFluxDofStartIndex = mmDofStartIndex * DIMENSIONS;

      Q = &lQhi[mmDofStartIndex    ];
      f = &lFhi[mmFluxDofStartIndex];

      for(int ivar=0; ivar < nvar; ivar++) {
        lQhbnd[dofStartIndexL+ivar] += dg::FLCoeff[mm] * Q[ivar];
        lQhbnd[dofStartIndexR+ivar] += dg::FRCoeff[mm] * Q[ivar];

        lFhbnd[dofStartIndexL+ivar] += dg::FLCoeff[mm] * f[ivar];
        lFhbnd[dofStartIndexR+ivar] += dg::FRCoeff[mm] * f[ivar];
      }
    }
    continue;
  }

  // y-direction: face 2 (left) and face 3 (right)
  for (int ii=0; ii<basisSize; ii++) {
    const int nodeIndex      = ii;
    const int dofStartIndexL = EXAHYPE_FACE_FRONT * numberOfFaceDof + nodeIndex * nvar;
    const int dofStartIndexR = EXAHYPE_FACE_BACK  * numberOfFaceDof + nodeIndex * nvar;

    for (int mm=0; mm<basisSize; mm++) {
      const int mmNodeIndex         = ii  + basisSize * mm;
      const int mmDofStartIndex     = mmNodeIndex * nvar;
      const int mmFluxDofStartIndex = mmDofStartIndex * DIMENSIONS;

      Q = &lQhi [mmDofStartIndex         ];
      g = &lFhi[mmFluxDofStartIndex+nvar];

      for(int ivar=0; ivar < nvar; ivar++) {
        lQhbnd[dofStartIndexL+ivar] += dg::FLCoeff[mm] * Q[ivar];
        lQhbnd[dofStartIndexR+ivar] += dg::FRCoeff[mm] * Q[ivar];

        lFhbnd[dofStartIndexL+ivar] += dg::FLCoeff[mm] * g[ivar];
        lFhbnd[dofStartIndexR+ivar] += dg::FRCoeff[mm] * g[ivar];
      }
    }
    continue;
  }

  // clean up
  std::free(tmp);
  std::free(rhs0);
  std::free(dqdt);
}

void exahype::mappings::SpaceTimePredictor::enterCell(
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

    for (int i=1; i<EXAHYPE_PATCH_SIZE_X+1; i++) { // loop over patches
      for (int j=1; j<EXAHYPE_PATCH_SIZE_Y+1; j++) {
        const int patchIndex = i + (EXAHYPE_PATCH_SIZE_X+2) * j;

        // space-time DoF (basisSize**(DIMENSIONS+1))
        double* lQi = &(DataHeap::getInstance().getData(cellDescription.getSpaceTimePredictor(patchIndex)) [0]._persistentRecords._u);
        double* lFi = &(DataHeap::getInstance().getData(cellDescription.getSpaceTimeVolumeFlux(patchIndex))[0]._persistentRecords._u);

        // volume DoF (basisSize**(DIMENSIONS))
        double* luh  = &(DataHeap::getInstance().getData(cellDescription.getSolution(patchIndex))  [0]._persistentRecords._u);
        double* lQhi = &(DataHeap::getInstance().getData(cellDescription.getPredictor(patchIndex)) [0]._persistentRecords._u);
        double* lFhi = &(DataHeap::getInstance().getData(cellDescription.getVolumeFlux(patchIndex))[0]._persistentRecords._u);

        // face DoF (basisSize**(DIMENSIONS-1))
        double* lQhbnd = &(DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor(patchIndex))[0]._persistentRecords._u);
        double* lFhbnd = &(DataHeap::getInstance().getData(cellDescription.getFluctuation(patchIndex))          [0]._persistentRecords._u);

        computePredictor(
            lQi,
            lFi,
            luh,
            lQhi,
            lFhi,
            lQhbnd,
            lFhbnd,
            center,
            dxPatch,dyPatch,
            this->_timeStepSize,
            patchIndex,
            nvar,
            basisSize-1);
      }
    }
  }


  logTraceOutWith1Argument( "enterCell(...)", fineGridCell );
}


void exahype::mappings::SpaceTimePredictor::leaveCell(
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


void exahype::mappings::SpaceTimePredictor::beginIteration(
    exahype::State&  solverState
) {
  logTraceInWith1Argument( "beginIteration(State)", solverState );

  // Begin of code for ADERDG method
  _timeStepSize = solverState.getTimeStepSize();
  // End of code for ADERDG method

  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


void exahype::mappings::SpaceTimePredictor::endIteration(
    exahype::State&  solverState
) {
  logTraceInWith1Argument( "endIteration(State)", solverState );
  // do nothing
  logTraceOutWith1Argument( "endIteration(State)", solverState);
}



void exahype::mappings::SpaceTimePredictor::descend(
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


void exahype::mappings::SpaceTimePredictor::ascend(
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
