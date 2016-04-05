#include "exahype/mappings/GlobalTimeStepComputation.h"

#include "peano/utils/Globals.h"

#include "tarch/multicore/Loop.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/Solver.h"

#include "tarch/multicore/Lock.h"

#include <limits>

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::GlobalTimeStepComputation::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::
          SendDataAndStateAfterLastTouchVertexLastTime,
      false);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification exahype::mappings::GlobalTimeStepComputation::
    touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification exahype::mappings::GlobalTimeStepComputation::
    touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::GlobalTimeStepComputation::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::GlobalTimeStepComputation::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::GlobalTimeStepComputation::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::GlobalTimeStepComputation::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::GlobalTimeStepComputation::_log( "exahype::mappings::GlobalTimeStepComputation" );

exahype::mappings::GlobalTimeStepComputation::GlobalTimeStepComputation() {
  // do nothing
}

exahype::mappings::GlobalTimeStepComputation::~GlobalTimeStepComputation() {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::prepareEmptyLocalTimeStepData() {
  _minTimeStepSizes.resize( exahype::solvers::RegisteredSolvers.size() );

  for (int i=0; i<static_cast<int>( exahype::solvers::RegisteredSolvers.size() ); i++) {
    _minTimeStepSizes[i] = std::numeric_limits<double>::max();
  }
}


void exahype::mappings::GlobalTimeStepComputation::mergeLocalTimeStepDataIntoSolvers() {
  for (int i=0; i<static_cast<int>( exahype::solvers::RegisteredSolvers.size() ); i++) {
    exahype::solvers::Solver* solver =
      exahype::solvers::RegisteredSolvers[i];

    // might be too restrictive later
    // assertion( _minTimeStepSizes[i] < std::numeric_limits<double>::max() );
//    logInfo( "mergeLocalTimeStepDataIntoSolvers()", "solver " << i << " is updated with time step size " << _minTimeStepSizes[i] );
    solver->updateMinNextPredictorTimeStepSize( _minTimeStepSizes[i] );
  }
}


#if defined(SharedMemoryParallelisation)
exahype::mappings::GlobalTimeStepComputation::GlobalTimeStepComputation(
  const GlobalTimeStepComputation& masterThread
) {
  prepareEmptyLocalTimeStepData();
}


void exahype::mappings::GlobalTimeStepComputation::mergeWithWorkerThread(
  const GlobalTimeStepComputation& workerThread
) {
  for (int i=0; i<static_cast<int>( exahype::solvers::RegisteredSolvers.size() ); i++) {
    _minTimeStepSizes[i] = std::min( _minTimeStepSizes[i], workerThread._minTimeStepSizes[i] );
  }
}
#endif


void exahype::mappings::GlobalTimeStepComputation::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::GlobalTimeStepComputation::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::GlobalTimeStepComputation::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  // do nothing
  return true;
}

void exahype::mappings::GlobalTimeStepComputation::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell
) {
/*
  double minSolverTimeStampAndTimeStepSize[2];

  for (
    std::vector<exahype::solvers::Solver*>::const_iterator p = exahype::solvers::RegisteredSolvers.begin();
    p != exahype::solvers::RegisteredSolvers.end();
    p++
  ) {
    minSolverTimeStampAndTimeStepSize[0] = (*p)->getMinCorrectorTimeStamp();
    minSolverTimeStampAndTimeStepSize[1] = (*p)->getMinCorrectorTimeStepSize();

    MPI_Send( minSolverTimeStampAndTimeStepSize, 2, MPI_DOUBLE, tarch::parallel::NodePool::getInstance().getMasterRank(), x, tarch::parallel::Node::getInstance().getCommunicator() );
  }
*/
}


void exahype::mappings::GlobalTimeStepComputation::mergeWithMaster(
    const exahype::Cell& workerGridCell,
    exahype::Vertex* const workerGridVertices,
    const peano::grid::VertexEnumerator& workerEnumerator,
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker, const exahype::State& workerState,
    exahype::State& masterState) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::GlobalTimeStepComputation::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::enterCell(
  exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
  const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
  exahype::Vertex* const coarseGridVertices,
  const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
  exahype::Cell& coarseGridCell,
  const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell
) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(fineGridCell.getADERDGCellDescriptionsIndex())) {
      const int numberOfADERDGCellDescriptions = static_cast<int>(
          ADERDGCellDescriptionHeap::getInstance()
      .getData(fineGridCell.getADERDGCellDescriptionsIndex())
      .size());

      // please use a different UserDefined per mapping/event
      const peano::datatraversal::autotuning::MethodTrace methodTrace =
          peano::datatraversal::autotuning::UserDefined1;
      const int grainSize =
          peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
              numberOfADERDGCellDescriptions, methodTrace);
      pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
      records::ADERDGCellDescription& p =
          ADERDGCellDescriptionHeap::getInstance().getData(
              fineGridCell.getADERDGCellDescriptionsIndex())[i];

      if (p.getType()==exahype::Cell::RealCell) {
        exahype::solvers::Solver* solver =
            exahype::solvers::RegisteredSolvers[p.getSolverNumber()];

        double* luh = DataHeap::getInstance().getData(p.getSolution()).data();

        double admissibleTimeStepSize = solver->
            stableTimeStepSize(luh, fineGridVerticesEnumerator.getCellSize());

        logDebug("enterCell(...)::dt_adm", admissibleTimeStepSize);  // Delta tp+1

        // direct update of the cell description time steps
        p.setCorrectorTimeStamp(p.getPredictorTimeStamp());
        p.setCorrectorTimeStepSize(p.getPredictorTimeStepSize());
        p.setPredictorTimeStamp(p.getPredictorTimeStamp() + admissibleTimeStepSize);
        p.setPredictorTimeStepSize(admissibleTimeStepSize);
        //p.setNextPredictorTimeStepSize(admissibleTimeStepSize);

        // todo 16/02/27:Dominic Etienne Charrier
        // in case we use optimistic time stepping:
        // if last predictor time step size is larger
        // as admissibleTimeStepSize + tolerance:
        // make sure that corrector time step size
        // will equal predictor time step size in next
        // sweep.
        // Extra attention must be paid to time stamps.
        // All this should be done by the solver.

        // indirect update of the solver time step sizes
        //  tarch::multicore::Lock lock(_semaphore);

        _minTimeStepSizes[ p.getSolverNumber() ] = std::min( admissibleTimeStepSize, _minTimeStepSizes[ p.getSolverNumber() ] );
        //logInfo( "enterCell()", "update local entry " << p->getSolverNumber() << " with " << admissibleTimeStepSize << ", set " << _minTimeStepSizes[ p->getSolverNumber() ] );
        //  lock.free();
      }
      endpfor peano::datatraversal::autotuning::Oracle::getInstance()
      .parallelSectionHasTerminated(methodTrace);
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::GlobalTimeStepComputation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}


void exahype::mappings::GlobalTimeStepComputation::beginIteration(
  exahype::State& solverState
) {
  prepareEmptyLocalTimeStepData();
}


void exahype::mappings::GlobalTimeStepComputation::endIteration(
  exahype::State& solverState
) {
  mergeLocalTimeStepDataIntoSolvers();
}


void exahype::mappings::GlobalTimeStepComputation::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::GlobalTimeStepComputation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
