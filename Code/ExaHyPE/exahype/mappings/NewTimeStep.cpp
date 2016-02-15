#include "exahype/mappings/NewTimeStep.h"

#include "exahype/solvers/Solve.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification   exahype::mappings::NewTimeStep::communicationSpecification() {
  return peano::CommunicationSpecification(peano::CommunicationSpecification::SendDataAndStateBeforeFirstTouchVertexFirstTime,peano::CommunicationSpecification::SendDataAndStateAfterLastTouchVertexLastTime,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::NewTimeStep::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::NewTimeStep::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::NewTimeStep::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::OnlyLeaves,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::NewTimeStep::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::NewTimeStep::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::NewTimeStep::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


tarch::logging::Log                exahype::mappings::NewTimeStep::_log( "exahype::mappings::NewTimeStep" ); 


exahype::mappings::NewTimeStep::NewTimeStep() {
  // do nothing
}


exahype::mappings::NewTimeStep::~NewTimeStep() {
  // do nothing
}


#if defined(SharedMemoryParallelisation)
exahype::mappings::NewTimeStep::NewTimeStep(const NewTimeStep&  masterThread):
    _localState(masterThread._localState)
{
  // do nothing
}


void exahype::mappings::NewTimeStep::mergeWithWorkerThread(const NewTimeStep& workerThread) {
  // do nothing
}
#endif


void exahype::mappings::NewTimeStep::createHangingVertex(
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


void exahype::mappings::NewTimeStep::destroyHangingVertex(
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


void exahype::mappings::NewTimeStep::createInnerVertex(
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


void exahype::mappings::NewTimeStep::createBoundaryVertex(
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


void exahype::mappings::NewTimeStep::destroyVertex(
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


void exahype::mappings::NewTimeStep::createCell(
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


void exahype::mappings::NewTimeStep::destroyCell(
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
void exahype::mappings::NewTimeStep::mergeWithNeighbour(
  exahype::Vertex&  vertex,
  const exahype::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
  // do nothing
}

void exahype::mappings::NewTimeStep::prepareSendToNeighbour(
  exahype::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
  // do nothing
}

void exahype::mappings::NewTimeStep::prepareCopyToRemoteNode(
  exahype::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
  // do nothing
}

void exahype::mappings::NewTimeStep::prepareCopyToRemoteNode(
  exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
) {
  // do nothing
}

void exahype::mappings::NewTimeStep::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Vertex&  localVertex,
  const exahype::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
  // do nothing
}

void exahype::mappings::NewTimeStep::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Cell&  localCell,
  const exahype::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                       level
) {
  // do nothing
}

bool exahype::mappings::NewTimeStep::prepareSendToWorker(
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

void exahype::mappings::NewTimeStep::prepareSendToMaster(
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


void exahype::mappings::NewTimeStep::mergeWithMaster(
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


void exahype::mappings::NewTimeStep::receiveDataFromMaster(
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


void exahype::mappings::NewTimeStep::mergeWithWorker(
  exahype::Cell&           localCell, 
  const exahype::Cell&     receivedMasterCell,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level
) {
  // do nothing
}


void exahype::mappings::NewTimeStep::mergeWithWorker(
  exahype::Vertex&        localVertex,
  const exahype::Vertex&  receivedMasterVertex,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
  // do nothing
}
#endif

void exahype::mappings::NewTimeStep::touchVertexFirstTime(
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

void exahype::mappings::NewTimeStep::touchVertexLastTime(
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

void exahype::mappings::NewTimeStep::enterCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  for (
    ADERDGCellDescriptionHeap::HeapEntries::iterator p = ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.getADERDGCellDescriptionsIndex()).begin();
    p != ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.getADERDGCellDescriptionsIndex()).end();
    p++
  ) {
    const exahype::State::shared_ptr_Solve solve = _localState.getSolveRegistry()[ p->getSolveNumber() ];

    if (solve->getTimeStepping()==exahype::solvers::Solve::GLOBAL) {
      p->setCorrectorTimeStamp   (solve->getCorrectorTimeStamp   ());
      p->setCorrectorTimeStepSize(solve->getCorrectorTimeStepSize());
      p->setPredictorTimeStamp   (solve->getPredictorTimeStamp   ());
      p->setPredictorTimeStepSize(solve->getPredictorTimeStepSize());
    }
    if (!solve->isCorrectorTimeLagging()) {
      p->setCorrectorTimeStamp   (p->getPredictorTimeStamp   ());
      p->setCorrectorTimeStepSize(p->getPredictorTimeStepSize());
    }

    // todo remove
    assertion(solve->getTimeStepping()==exahype::solvers::Solve::GLOBAL);
    assertion(solve->isCorrectorTimeLagging()==true);
    assertion(p->getCorrectorTimeStamp()   ==solve->getCorrectorTimeStamp());
    assertion(p->getCorrectorTimeStepSize()==solve->getCorrectorTimeStepSize());
  }
}


void exahype::mappings::NewTimeStep::leaveCell(
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


void exahype::mappings::NewTimeStep::beginIteration(
  exahype::State&  solverState
) {
  logTraceInWith1Argument( "beginIteration(State)", solverState );

  _localState = solverState;

  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


void exahype::mappings::NewTimeStep::endIteration(
  exahype::State&  solverState
) {
  // do nothing
}



void exahype::mappings::NewTimeStep::descend(
  exahype::Cell * const          fineGridCells,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell
) {
  // do nothing
}


void exahype::mappings::NewTimeStep::ascend(
  exahype::Cell * const    fineGridCells,
  exahype::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell
) {
  // do nothing
}

const exahype::State& exahype::mappings::NewTimeStep::getState() const {
  return _localState;
}
