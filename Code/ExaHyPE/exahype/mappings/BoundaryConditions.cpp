#include "exahype/mappings/BoundaryConditions.h"

#include "tarch/multicore/Lock.h"

#include "peano/utils/Globals.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/Constants.h"
#include "exahype/aderdg/ADERDG.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification   exahype::mappings::BoundaryConditions::communicationSpecification() {
  return peano::CommunicationSpecification(peano::CommunicationSpecification::SendDataAndStateBeforeFirstTouchVertexFirstTime,peano::CommunicationSpecification::SendDataAndStateAfterLastTouchVertexLastTime,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::BoundaryConditions::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::OnlyLeaves,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::BoundaryConditions::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::BoundaryConditions::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::BoundaryConditions::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::BoundaryConditions::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::BoundaryConditions::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


tarch::logging::Log                exahype::mappings::BoundaryConditions::_log( "exahype::mappings::BoundaryConditions" ); 


exahype::mappings::BoundaryConditions::BoundaryConditions() {
  // do nothing
}


exahype::mappings::BoundaryConditions::~BoundaryConditions() {
  // do nothing
}


#if defined(SharedMemoryParallelisation)
exahype::mappings::BoundaryConditions::BoundaryConditions(const BoundaryConditions&  masterThread) {
  logTraceIn( "RiemannSolver(RiemannSolver)" );

  _localState.setOldTimeStepSize(masterThread.getState().getOldTimeStepSize());

  logTraceOut( "RiemannSolver(RiemannSolver)" );
}


void exahype::mappings::BoundaryConditions::mergeWithWorkerThread(const BoundaryConditions& workerThread) {
  // do nothing
}
#endif


void exahype::mappings::BoundaryConditions::createHangingVertex(
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


void exahype::mappings::BoundaryConditions::destroyHangingVertex(
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


void exahype::mappings::BoundaryConditions::createInnerVertex(
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


void exahype::mappings::BoundaryConditions::createBoundaryVertex(
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


void exahype::mappings::BoundaryConditions::destroyVertex(
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


void exahype::mappings::BoundaryConditions::createCell(
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


void exahype::mappings::BoundaryConditions::destroyCell(
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
void exahype::mappings::BoundaryConditions::mergeWithNeighbour(
    exahype::Vertex&  vertex,
    const exahype::Vertex&  neighbour,
    int                                           fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::prepareSendToNeighbour(
    exahype::Vertex&  vertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::prepareCopyToRemoteNode(
    exahype::Vertex&  localVertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::prepareCopyToRemoteNode(
    exahype::Cell&  localCell,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex&  localVertex,
    const exahype::Vertex&  masterOrWorkerVertex,
    int                                       fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&  x,
    const tarch::la::Vector<DIMENSIONS,double>&  h,
    int                                       level
) {
  // do nothing
}

void exahype::mappings::BoundaryConditions::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell&  localCell,
    const exahype::Cell&  masterOrWorkerCell,
    int                                       fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                       level
) {
  // do nothing
}

bool exahype::mappings::BoundaryConditions::prepareSendToWorker(
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
    int                                                                  worker
) {
  // do nothing
  return true;
}

void exahype::mappings::BoundaryConditions::prepareSendToMaster(
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


void exahype::mappings::BoundaryConditions::mergeWithMaster(
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


void exahype::mappings::BoundaryConditions::receiveDataFromMaster(
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


void exahype::mappings::BoundaryConditions::mergeWithWorker(
    exahype::Cell&           localCell,
    const exahype::Cell&     receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                          level
) {
  // do nothing
}


void exahype::mappings::BoundaryConditions::mergeWithWorker(
    exahype::Vertex&        localVertex,
    const exahype::Vertex&  receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}
#endif

void exahype::mappings::BoundaryConditions::touchVertexFirstTime(
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

void exahype::mappings::BoundaryConditions::applyBoundaryConditions(
    tarch::la::Vector<TWO_POWER_D,int>& adjacentADERDGCellDescriptionsIndices,
    const int cellIndexL,
    const int cellIndexR,
    const int faceIndexL,
    const int faceIndexR,
    const int numberOfFaceDof,
    const double * const normal
) {
  bool noBoundaryFace       = true;
  bool bothAreBoundaryFaces = false;
  int  cellIndex            = multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex;
  int  faceIndex            = -1;

  if ( // If cellIndexL refers to a boundary index then cellIndexR refers to a cell.
      adjacentADERDGCellDescriptionsIndices[cellIndexL] == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex
  ) {
    noBoundaryFace       = false;
    bothAreBoundaryFaces = true;
    cellIndex            = cellIndexR;
    faceIndex            = faceIndexR;
  }

  if (
      adjacentADERDGCellDescriptionsIndices[cellIndexR] == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex
  ) {
    noBoundaryFace       = false;
    bothAreBoundaryFaces = bothAreBoundaryFaces & true;
    cellIndex            = cellIndexL;
    faceIndex            = faceIndexL;
  }

  // Only continue if this is a boundary face. See multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex.
  if (
      noBoundaryFace
      ||
      bothAreBoundaryFaces
  ) {
    return;
  }

  assertion1(
      adjacentADERDGCellDescriptionsIndices[cellIndexL] > multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex
      ||
      adjacentADERDGCellDescriptionsIndices[cellIndexR] > multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex
      ,
      adjacentADERDGCellDescriptionsIndices.toString()
  );


  records::ADERDGCellDescription&  cellDescription =
      ADERDGADERDGCellDescriptionHeap::getInstance().getData(adjacentADERDGCellDescriptionsIndices[cellIndex])[0];

  bool riemannSolveNotPerformed = false;
  {
    // Lock the critical multithreading area.
    // Note that two boundary vertices can operate on the same face at the same time.
    tarch::multicore::Lock lock(_semaphore);
    riemannSolveNotPerformed = !cellDescription.getRiemannSolvePerformed(faceIndex);

    if(riemannSolveNotPerformed) {
      cellDescription.setRiemannSolvePerformed(faceIndex,true);
    }
  } // Unlock the critical multithreading area by letting lock go out of scope.

  // Apply the boundary conditions.
  // todo Copy and periodic boundary conditions should be configured by additional mappings after
  // the initial grid refinement and the initialisation of the cell descriptions.
  // The configuration should only involve an edit of the index maps generated
  // by the multiscalelinkedcell toolbox. In case of consequent mesh-refinement,
  // these indices should be propagated down to the finer cells.
  // The resulting Riemann problems are then simply solved
  // by exahype::mappings::RiemannSolver.
  if (riemannSolveNotPerformed) {
    constexpr double superfluousArgument = 0;
    // work vectors
    double QavL   [EXAHYPE_NVARS]; // av: average
    double QavR   [EXAHYPE_NVARS]; // need Qav twice!
    double lambdaL[EXAHYPE_NVARS];
    double lambdaR[EXAHYPE_NVARS];

    double * Qhbnd = &(DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor())[faceIndex * numberOfFaceDof]._persistentRecords._u);
    double * Fhbnd = &(DataHeap::getInstance().getData(cellDescription.getFluctuation())          [faceIndex * numberOfFaceDof]._persistentRecords._u);

    // Invoke user defined boundary condition function todo
    // At the moment, we simply copy the cell solution to the boundary.
    aderdg::riemannSolver<DIMENSIONS>(
        Fhbnd,
        Fhbnd,
        Qhbnd,
        Qhbnd,
        QavL,
        QavR,
        lambdaL,
        lambdaR,
        _localState.getMaxTimeStepSize(),

        superfluousArgument,
        normal);
  }
}

void exahype::mappings::BoundaryConditions::touchVertexLastTime(
    exahype::Vertex&         fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "touchVertexLastTime(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );

  if (
      fineGridVertex.getRefinementControl()==Vertex::Records::Unrefined // todo Replace by something that works for multiple PDEs. Must possible move into applyBoundaryConditions.
  ) {
    assertion1WithExplanation(_localState.getOldTimeStepSize() < std::numeric_limits<double>::max(),_localState.getOldTimeStepSize(),"Old time step size was not initialised correctly!");

    tarch::la::Vector<TWO_POWER_D,int>& adjacentADERDGCellDescriptionsIndices = fineGridVertex.getADERDGCellDescriptionsIndex();
    // todo: DEC: Reverse engineered indices from
    // PatchInitialisation2MultiscaleLinkedCell_1::touchVertexLastTime(...)
    // Not sure what happens with hanging nodes.

    /*
     * Right cell-left cell   pair indices: 0,1; 2,3;   4,5; 6;7
     * Front cell-back cell   pair indices: 0,2; 1,3;   4,6; 5;7
     * Top   cell-bottom cell pair indices: 0,4; 1,5;   2,6; 3;7
     *
     * Note that from the viewpoint of a cell, the face
     * has always the "opposite" index, i.e., we solve a Riemann
     * problem on the left face of the right cell (which
     * is the right face of the left cell).
     */

    constexpr int basisSize   = EXAHYPE_ORDER+1;
    constexpr int nvar        = EXAHYPE_NVARS;
    const int numberOfFaceDof = nvar * tarch::la::aPowI(DIMENSIONS-1,basisSize);

    // index maps (
    constexpr int cellIndicesLeft   [4] = { 0, 2, 4, 6 };
    constexpr int cellIndicesRight  [4] = { 1, 3, 5, 7 };
    constexpr int cellIndicesFront  [4] = { 0, 1, 4, 5 };
    constexpr int cellIndicesBack   [4] = { 2, 3, 6, 7 };
#if DIMENSIONS==3
    constexpr int cellIndicesBottom [4] = { 0, 1, 2, 3 };
    constexpr int cellIndicesTop    [4] = { 4, 5, 6, 7 };
#endif

    // normal vectors
    const double nx[3]= { 1., 0., 0. };
    const double ny[3]= { 0., 1., 0. };
#if DIMENSIONS==3
    const double nz[3]= { 0., 0., 1. };
#endif

    // Left/right face
    for (int i=0; i<TWO_POWER_D_DIVIDED_BY_TWO; i++) {
      applyBoundaryConditions(
          adjacentADERDGCellDescriptionsIndices,
          cellIndicesLeft [i],
          cellIndicesRight[i],
          EXAHYPE_FACE_RIGHT,
          EXAHYPE_FACE_LEFT,
          numberOfFaceDof,
          nx);

      applyBoundaryConditions(
          adjacentADERDGCellDescriptionsIndices,
          cellIndicesFront[i],
          cellIndicesBack [i],
          EXAHYPE_FACE_BACK,
          EXAHYPE_FACE_FRONT,
          numberOfFaceDof,
          ny);

#if DIMENSIONS==3
      applyBoundaryConditions(
          adjacentADERDGCellDescriptionsIndices,
          cellIndicesBottom[i],
          cellIndicesTop   [i],
          EXAHYPE_FACE_TOP,
          EXAHYPE_FACE_BOTTOM,
          numberOfFaceDof,
          nz);
#endif
    }
  }
  logTraceOutWith1Argument( "touchVertexLastTime(...)", fineGridVertex );
}

void exahype::mappings::BoundaryConditions::enterCell(
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


void exahype::mappings::BoundaryConditions::leaveCell(
    exahype::Cell&           fineGridCell,
    exahype::Vertex * const  fineGridVertices,
    const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  // do nothing
}


void exahype::mappings::BoundaryConditions::beginIteration(
    exahype::State&  solverState
) {
  logTraceInWith1Argument( "beginIteration(State)", solverState );

  _localState = solverState;

  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


void exahype::mappings::BoundaryConditions::endIteration(
    exahype::State&  solverState
) {
  // do nothing
}



void exahype::mappings::BoundaryConditions::descend(
    exahype::Cell * const          fineGridCells,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell
) {
  // do nothing
}


void exahype::mappings::BoundaryConditions::ascend(
    exahype::Cell * const    fineGridCells,
    exahype::Vertex * const  fineGridVertices,
    const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell
) {
  // do nothing
}
