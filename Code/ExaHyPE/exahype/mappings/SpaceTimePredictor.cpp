#include "exahype/mappings/SpaceTimePredictor.h"

#include "peano/utils/Globals.h"

#include "tarch/multicore/Loop.h"
#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/Solver.h"


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
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::SpaceTimePredictor::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::SpaceTimePredictor::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::OnlyLeaves,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::SpaceTimePredictor::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::SpaceTimePredictor::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::SpaceTimePredictor::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


tarch::logging::Log                exahype::mappings::SpaceTimePredictor::_log( "exahype::mappings::SpaceTimePredictor" ); 


exahype::mappings::SpaceTimePredictor::SpaceTimePredictor() {
  // do nothing
}


exahype::mappings::SpaceTimePredictor::~SpaceTimePredictor() {
  // do nothing
}


#if defined(SharedMemoryParallelisation)
exahype::mappings::SpaceTimePredictor::SpaceTimePredictor(const SpaceTimePredictor&  masterThread):
      _localState(masterThread._localState) {
}


void exahype::mappings::SpaceTimePredictor::mergeWithWorkerThread(const SpaceTimePredictor& workerThread) {
  // do nothing
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
  // do nothing
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
  // do nothing
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
  // do nothing
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
  // do nothing
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
  // do nothing
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
  // do nothing
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
  // do nothing
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
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::prepareSendToNeighbour(
    exahype::Vertex&  vertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::prepareCopyToRemoteNode(
    exahype::Vertex&  localVertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::prepareCopyToRemoteNode(
    exahype::Cell&  localCell,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex&  localVertex,
    const exahype::Vertex&  masterOrWorkerVertex,
    int                                       fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&  x,
    const tarch::la::Vector<DIMENSIONS,double>&  h,
    int                                       level
) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell&  localCell,
    const exahype::Cell&  masterOrWorkerCell,
    int                                       fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                       level
) {
  // do nothing
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
  // do nothing
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
  // do nothing
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
  // do nothing
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
  // do nothing
}


void exahype::mappings::SpaceTimePredictor::mergeWithWorker(
    exahype::Cell&           localCell,
    const exahype::Cell&     receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                          level
) {
  // do nothing
}


void exahype::mappings::SpaceTimePredictor::mergeWithWorker(
    exahype::Vertex&        localVertex,
    const exahype::Vertex&  receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
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
  // do nothing
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
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::enterCell(
    exahype::Cell&                            fineGridCell,
    exahype::Vertex * const                   fineGridVertices,
    const peano::grid::VertexEnumerator&      fineGridVerticesEnumerator,
    exahype::Vertex * const                   coarseGridVertices,
    const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
    exahype::Cell&                            coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&  fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "enterCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );

  // todo 17/02/16:Dominic Etienne Charrier: Assertion not valid anymore if adaptive mesh refinement is switches on.
  assertion1(
      fineGridCell.isRefined() || !ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.getADERDGCellDescriptionsIndex()).empty(),
      fineGridCell.toString()
  );

  const int numberOfADERDGCellDescriptions = static_cast<int>( ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.getADERDGCellDescriptionsIndex()).size() );
  const peano::datatraversal::autotuning::MethodTrace methodTrace = peano::datatraversal::autotuning::UserDefined4; // Dominic, please use a different UserDefined per mapping/event. There should be enough by now.
  const int  grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfADERDGCellDescriptions,methodTrace);
  pfor(i,0,numberOfADERDGCellDescriptions,grainSize)
    // todo ugly
    // This is not beautiful and should be replaced by a reference next. I just
    // use it to mirror the aforementioned realisation. Dominic, please change
    // successively to a simpler scheme with just references. Pointers are
    // ugly.
    records::ADERDGCellDescription* p = &(ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.getADERDGCellDescriptionsIndex())[i]);

    exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers [ p->getSolverNumber() ];

    // space-time DoF (basisSize**(DIMENSIONS+1))
    double * lQi = &(DataHeap::getInstance().getData(p->getSpaceTimePredictor()) [0]._persistentRecords._u);
    double * lFi = &(DataHeap::getInstance().getData(p->getSpaceTimeVolumeFlux())[0]._persistentRecords._u);

    // volume DoF (basisSize**(DIMENSIONS))
    double * luh  = &(DataHeap::getInstance().getData(p->getSolution())  [0]._persistentRecords._u);
    double * lQhi = &(DataHeap::getInstance().getData(p->getPredictor()) [0]._persistentRecords._u);
    double * lFhi = &(DataHeap::getInstance().getData(p->getVolumeFlux())[0]._persistentRecords._u);

    // face DoF (basisSize**(DIMENSIONS-1))
    double * lQhbnd = &(DataHeap::getInstance().getData(p->getExtrapolatedPredictor())[0]._persistentRecords._u);
    double * lFhbnd = &(DataHeap::getInstance().getData(p->getFluctuation())          [0]._persistentRecords._u);

    solver->spaceTimePredictor(
        lQi,
        lFi,
        lQhi,
        lFhi,
        lQhbnd,  // da kommt was drauf
        lFhbnd,  // da kommt was drauf
        luh,
        fineGridVerticesEnumerator.getCellSize(),
        p->getPredictorTimeStepSize()
    );

    logDebug("enterCell(...)::debug::after::dt_max",_localState.getCurrentMinTimeStepSize());
    logDebug("enterCell(...)::debug::after::luh[0]",luh[0]);
    logDebug("enterCell(...)::debug::after::lQi[0]",lQi[0]);
    logDebug("enterCell(...)::debug::after::lFi[0]",lFi[0]);
    logDebug("enterCell(...)::debug::after::lQhi[0]",lQhi[0]);
    logDebug("enterCell(...)::debug::after::lFhi[0]",lFhi[0]);
    logDebug("enterCell(...)::debug::after::lQhbnd[0]",lQhbnd[0]);
    logDebug("enterCell(...)::debug::after::lFhbnd[0]",lFhbnd[0]);
  endpfor
  peano::datatraversal::autotuning::Oracle::getInstance().parallelSectionHasTerminated(methodTrace);

  logTraceOutWith1Argument( "enterCell(...)", fineGridCell );
}


void exahype::mappings::SpaceTimePredictor::leaveCell(
    exahype::Cell&           fineGridCell,
    exahype::Vertex * const  fineGridVertices,
    const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&     fineGridPositionOfCell
) {
  // do nothing
}


void exahype::mappings::SpaceTimePredictor::beginIteration(
    exahype::State&  solverState
) {
  logTraceInWith1Argument( "beginIteration(State)", solverState );

  _localState = solverState;

  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


void exahype::mappings::SpaceTimePredictor::endIteration(
    exahype::State&  solverState
) {
  // do nothing
}



void exahype::mappings::SpaceTimePredictor::descend(
    exahype::Cell * const          fineGridCells,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell
) {
  // do nothing
}


void exahype::mappings::SpaceTimePredictor::ascend(
    exahype::Cell * const    fineGridCells,
    exahype::Vertex * const  fineGridVertices,
    const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell
) {
  // do nothing
}

const exahype::State& exahype::mappings::SpaceTimePredictor::getState() const {
  return _localState;
}
