#include "exahype/mappings/BoundaryConditions.h"

#include "peano/utils/Globals.h"

#include "tarch/multicore/Lock.h"
#include "tarch/multicore/Loop.h"
#include "peano/datatraversal/autotuning/Oracle.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/solvers/Solve.h"
#include "exahype/solvers/Solver.h"

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
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::BoundaryConditions::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::OnlyLeaves,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
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
exahype::mappings::BoundaryConditions::BoundaryConditions(const BoundaryConditions&  masterThread):
      _localState( masterThread._localState ) {
  _localState.deepCopySolveRegistry ( masterThread._localState );
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
  logTraceInWith6Arguments( "touchVertextFirstTime(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );

  assertion1WithExplanation(_localState.getPreviousMinTimeStepSize() < std::numeric_limits<double>::max(),_localState.toString(),"Old time step size is not initialised correctly!");

  tarch::la::Vector<TWO_POWER_D,int>& adjacentADERDGCellDescriptionsIndices = fineGridVertex.getADERDGCellDescriptionsIndex();
  // todo: DEC: Reverse engineered indices from
  // PatchInitialisation2MultiscaleLinkedCell_1::touchVertextFirstTime(...)
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
  // index maps (
  constexpr int cellIndicesLeft   [4] = { 0, 2, 4, 6 };
  constexpr int cellIndicesRight  [4] = { 1, 3, 5, 7 };
  constexpr int cellIndicesFront  [4] = { 0, 1, 4, 5 };
  constexpr int cellIndicesBack   [4] = { 2, 3, 6, 7 };
#if DIMENSIONS==3
  constexpr int cellIndicesBottom [4] = { 0, 1, 2, 3 };
  constexpr int cellIndicesTop    [4] = { 4, 5, 6, 7 };
#endif

  // Left/right face
  for (int i=0; i<TWO_POWER_D_DIVIDED_BY_TWO; i++) {
    applyBoundaryConditions(
        adjacentADERDGCellDescriptionsIndices,
        cellIndicesLeft [i],
        cellIndicesRight[i],
        EXAHYPE_FACE_RIGHT,
        EXAHYPE_FACE_LEFT,
        0);

    applyBoundaryConditions(
        adjacentADERDGCellDescriptionsIndices,
        cellIndicesFront[i],
        cellIndicesBack [i],
        EXAHYPE_FACE_BACK,
        EXAHYPE_FACE_FRONT,
        1);

#if DIMENSIONS==3
    applyBoundaryConditions(
        adjacentADERDGCellDescriptionsIndices,
        cellIndicesBottom[i],
        cellIndicesTop   [i],
        EXAHYPE_FACE_TOP,
        EXAHYPE_FACE_BOTTOM,
        2);
#endif
  }
  logTraceOutWith1Argument( "touchVertextFirstTime(...)", fineGridVertex );
}

void exahype::mappings::BoundaryConditions::applyBoundaryConditions(
    tarch::la::Vector<TWO_POWER_D,int>& adjacentADERDGCellDescriptionsIndices,
    const int cellIndexL,
    const int cellIndexR,
    const int faceIndexL,
    const int faceIndexR,
    const int normalNonZero
)
{
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

//  for (
//      ADERDGCellDescriptionHeap::HeapEntries::iterator p = ADERDGCellDescriptionHeap::getInstance().getData( adjacentADERDGCellDescriptionsIndices[cellIndex] ).begin();
//      p != ADERDGCellDescriptionHeap::getInstance().getData( adjacentADERDGCellDescriptionsIndices[cellIndex] ).end();
//      p++
//  ) {
  const auto numberOfADERDGCellDescriptions = ADERDGCellDescriptionHeap::getInstance().getData( adjacentADERDGCellDescriptionsIndices[cellIndex] ).size();
  const peano::datatraversal::autotuning::MethodTrace methodTrace = peano::datatraversal::autotuning::UserDefined0; // Dominic, please use a different UserDefined per mapping/event. There should be enough by now.
  const int  grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfADERDGCellDescriptions,methodTrace);
  pfor(i,0,numberOfADERDGCellDescriptions,grainSize)
    // todo ugly
    // This is not beautiful and should be replaced by a reference next. I just
    // use it to mirror the aforementioned realisation. Dominic, please change
    // successively to a simpler scheme with just references. Pointers are
    // ugly.
    records::ADERDGCellDescription* p = &(ADERDGCellDescriptionHeap::getInstance().getData( adjacentADERDGCellDescriptionsIndices[cellIndex] )[i]);

    exahype::solvers::Solve& solve   = _localState.getSolveRegistry()       [ p->getSolveNumber() ];
    exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers  [ solve.getSolverNumber() ];

    // Lock the critical multithreading area.
    // Note that two boundary vertices can operate on the same face at the same time.
    bool riemannSolveNotPerformed = false;
    tarch::multicore::Lock lock(_semaphore);
    riemannSolveNotPerformed = !p->getRiemannSolvePerformed(faceIndex);
    if(riemannSolveNotPerformed) {
      p->setRiemannSolvePerformed(faceIndex,true);
    }
    lock.free();

    // Apply the boundary conditions.
    // todo Copy and periodic boundary conditions should be configured by additional mappings after
    // the initial grid refinement and the initialisation of the cell descriptions.
    // The configuration should only involve an edit of the index maps generated
    // by the multiscalelinkedcell toolbox. In case of consequent mesh-refinement,
    // these indices should be propagated down to the finer cells.
    // The resulting Riemann problems are then simply solved
    // by exahype::mappings::RiemannSolver.
    if (riemannSolveNotPerformed) {
      // @todo 03/02/16:Dominic Etienne Charrier
      // Change to solver->getUnknownsPerFace()
      const int numberOfFaceDof = solver->getNumberOfVariables() * tarch::la::aPowI(DIMENSIONS-1,solver->getNodesPerCoordinateAxis());

      double * Qhbnd = &(DataHeap::getInstance().getData(p->getExtrapolatedPredictor())[faceIndex * numberOfFaceDof]._persistentRecords._u);
      double * Fhbnd = &(DataHeap::getInstance().getData(p->getFluctuation())          [faceIndex * numberOfFaceDof]._persistentRecords._u);

      // todo Dominic Charrier
      // Do not solve a Riemann problem here:
      // Invoke user defined boundary condition function
      // At the moment, we simply copy the cell solution to the boundary.

      synchroniseTimeStepping(*p);

      logDebug("touchVertexLastTime(...)::debug::before::dt_min(previous ) of State*",_localState.getPreviousMinTimeStepSize());
      logDebug("touchVertexLastTime(...)::debug::before::dt_min(corrector) of Solve*",solve.getCorrectorTimeStepSize());
      logDebug("touchVertexLastTime(...)::debug::before::Qhbnd[0]*",Qhbnd[0]);
      logDebug("touchVertexLastTime(...)::debug::before::Fhbnd[0]",Fhbnd[0]);

      solver->riemannSolver(
          Fhbnd,
          Fhbnd,
          Qhbnd,
          Qhbnd,
          p->getCorrectorTimeStepSize(),//solve.getCorrectorTimeStepSize(),//_localState.getPreviousMinTimeStepSize(),
          normalNonZero);

      logDebug("touchVertexLastTime(...)::debug::after::Qhbnd[0]*",Qhbnd[0]);
      logDebug("touchVertexLastTime(...)::debug::after::Fhbnd[0]",Fhbnd[0]);
    }
//  }
  endpfor
  peano::datatraversal::autotuning::Oracle::getInstance().parallelSectionHasTerminated(methodTrace);
}

void exahype::mappings::BoundaryConditions::synchroniseTimeStepping(records::ADERDGCellDescription& p) {
  const exahype::solvers::Solve& solve = _localState.getSolveRegistry()[ p.getSolveNumber() ];
  if (solve.getTimeStepping()==exahype::solvers::Solve::GLOBAL) {
    p.setCorrectorTimeStamp   (solve.getCorrectorTimeStamp   ());
    p.setCorrectorTimeStepSize(solve.getCorrectorTimeStepSize());
    p.setPredictorTimeStamp   (solve.getPredictorTimeStamp   ());
    p.setPredictorTimeStepSize(solve.getPredictorTimeStepSize());

    assertionNumericalEquals1(p.getCorrectorTimeStamp()   ,solve.getCorrectorTimeStamp(),   1e-12); // todo precision
    assertionNumericalEquals1(p.getCorrectorTimeStepSize(),solve.getCorrectorTimeStepSize(),1e-12);
    assertionNumericalEquals1(p.getPredictorTimeStamp()   ,solve.getPredictorTimeStamp(),   1e-12);
    assertionNumericalEquals1(p.getPredictorTimeStepSize(),solve.getPredictorTimeStepSize(),1e-12);
  }
  if (!solve.isCorrectorTimeLagging()) {
    p.setCorrectorTimeStamp   (p.getPredictorTimeStamp   ());
    p.setCorrectorTimeStepSize(p.getPredictorTimeStepSize());
  }

#if defined(Debug) || defined(Asserts)
  if (solve.getTimeStepping()==exahype::solvers::Solve::GLOBAL && !solve.isCorrectorTimeLagging()) {
    // Note that the solve time stamps and time step sizes are not modified if corrector time lagging
    // is deactivated. Thus, solve.getPredictor... and solve.getCorrector... are not the same in general
    // for any value of solve.isCorrectorTimeLagging().
    assertionNumericalEquals1(p.getPredictorTimeStamp()   ,solve.getPredictorTimeStamp(),   1e-12); // todo precision
    assertionNumericalEquals1(p.getPredictorTimeStepSize(),solve.getPredictorTimeStepSize(),1e-12);
    assertionNumericalEquals1(p.getCorrectorTimeStamp()   ,solve.getPredictorTimeStamp(),   1e-12);
    assertionNumericalEquals1(p.getCorrectorTimeStepSize(),solve.getPredictorTimeStepSize(),1e-12);
  }
#endif
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
  // do nothing
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
