#include "mpibalancing/GreedyBalancing.h"

#include "tarch/Assertions.h"
#include "tarch/la/ScalarOperations.h"
#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "peano/parallel/loadbalancing/Oracle.h"


tarch::logging::Log mpibalancing::GreedyBalancing::_log( "mpibalancing::GreedyBalancing" );


bool mpibalancing::GreedyBalancing::_forkHasFailed = false;


mpibalancing::GreedyBalancing::GreedyBalancing(int coarsestLevelWithRealWork, int coarsestRegularInnerAndOuterGridLevel):
  _coarsestLevelWithRealWork(coarsestLevelWithRealWork),
  _coarsestRegularInnerAndOuterGridLevel(coarsestRegularInnerAndOuterGridLevel) {
  _finestLevelToForkAggressively = _coarsestLevelWithRealWork;
  int activeNodes                = 1;

  while ( activeNodes<tarch::parallel::Node::getInstance().getNumberOfNodes() ) {
    activeNodes *= THREE_POWER_D;
    _finestLevelToForkAggressively++;
  }
}


mpibalancing::GreedyBalancing::~GreedyBalancing() {
}


void mpibalancing::GreedyBalancing::receivedStartCommand( const int commandFromMaster ) {
}


int mpibalancing::GreedyBalancing::getCommandForWorker( int workerRank, bool forkIsAllowed, bool joinIsAllowed ) {
  logTraceInWith3Arguments( "getCommandForWorker(int,bool)", workerRank, forkIsAllowed, joinIsAllowed);

  int result = peano::parallel::loadbalancing::ForkGreedy;

  if (_forkHasFailed) {
    result = peano::parallel::loadbalancing::Continue;
  }
  else if (_workersLevel.count(workerRank)==1) {
    const int workersLevel = _workersLevel.count(workerRank);

    if (workersLevel>=_finestLevelToForkAggressively) {
      result = peano::parallel::loadbalancing::ForkOnce;
    }
  }

  logTraceOutWith1Argument( "getCommandForWorker(int,bool)", result );
  return result;
}


void mpibalancing::GreedyBalancing::receivedTerminateCommand(
  int     workerRank,
  double  workerNumberOfInnerVertices,
  double  workerNumberOfBoundaryVertices,
  double  workerNumberOfOuterVertices,
  double  workerNumberOfInnerCells,
  double  workerNumberOfOuterCells,
  int     workerMaxLevel,
  int     currentLevel,
  const tarch::la::Vector<DIMENSIONS,double>& boundingBoxOffset,
  const tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize,
  bool    workerCouldNotEraseDueToDecomposition
) {
  _workersLevel[workerRank] = currentLevel;
}


void mpibalancing::GreedyBalancing::plotStatistics() {
}


peano::parallel::loadbalancing::OracleForOnePhase* mpibalancing::GreedyBalancing::createNewOracle(int adapterNumber) const {
  return new GreedyBalancing(_coarsestLevelWithRealWork, _coarsestRegularInnerAndOuterGridLevel);
}


void mpibalancing::GreedyBalancing::forkFailed() {
  if (!_forkHasFailed) {
    logInfo(
      "forkFailed()",
      "oracle was informed that fork has failed. No further fork attempts"
    );
  }
  _forkHasFailed = true;
}


int mpibalancing::GreedyBalancing::getCoarsestRegularInnerAndOuterGridLevel() const {
  return _coarsestRegularInnerAndOuterGridLevel;
}
