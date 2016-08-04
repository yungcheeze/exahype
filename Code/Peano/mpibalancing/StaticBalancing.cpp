#include "mpibalancing/StaticBalancing.h"

#include "tarch/Assertions.h"
#include "tarch/la/ScalarOperations.h"
#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "peano/parallel/loadbalancing/Oracle.h"


tarch::logging::Log mpibalancing::StaticBalancing::_log( "mpibalancing::StaticBalancing" );


bool mpibalancing::StaticBalancing::_forkHasFailed = false;


/*
int mpibalancing::StaticBalancing::getDefaultFinestAdministrativeLevel() {
  int result = 0;
  int regularCellsPerLevel = 1;
  while ( regularCellsPerLevel < tarch::parallel::Node::getInstance().getNumberOfNodes() ) {
    result++;
    regularCellsPerLevel += tarch::la::aPowI(DIMENSIONS*result,3);
  }
  result--;
  return result;
}
*/


mpibalancing::StaticBalancing::StaticBalancing(bool joinsAllowed, int coarsestRegularInnerAndOuterGridLevel):
  _coarsestRegularInnerAndOuterGridLevel(coarsestRegularInnerAndOuterGridLevel),
  _joinsAllowed(joinsAllowed),
  _criticalWorker(),
  _maxForksOnCriticalWorker(THREE_POWER_D) {
  _weightMap.insert( std::pair<int,double>(tarch::parallel::Node::getInstance().getRank(), 1.0) );
  _workerCouldNotEraseDueToDecomposition.insert( std::pair<int,bool>(tarch::parallel::Node::getInstance().getRank(), false) );
}


mpibalancing::StaticBalancing::~StaticBalancing() {
}


double mpibalancing::StaticBalancing::getMaximumWeightOfWorkers() const {
  double maximumWeight = std::numeric_limits<double>::min();
  for ( std::map<int,double>::const_iterator p=_weightMap.begin(); p!=_weightMap.end(); p++ ) {
    if ( p->second > maximumWeight ) {
      maximumWeight = p->second;
    }
  }

  assertion1( maximumWeight>=1.0, maximumWeight );

  return maximumWeight;
}


double mpibalancing::StaticBalancing::getMinimumWeightOfWorkers() const {
  double minimumWeight  = std::numeric_limits<double>::max();
  for ( std::map<int,double>::const_iterator p=_weightMap.begin(); p!=_weightMap.end(); p++ ) {
    if ( p->second < minimumWeight && p->second > 1.0) {
      minimumWeight = p->second;
    }
  }
  assertion1( minimumWeight>=1.0, minimumWeight );
  return minimumWeight;
}


void mpibalancing::StaticBalancing::identifyCriticalPathes( int commandFromMaster ) {
  /**
   * We consider the CritialPathThreshold upper percent of the workers to be
   * critical.
   */
  const double CritialPathThreshold = 0.1;

  /**
   * If minimum and critical weight differ by less than LocalMinimumThreshold
   * percent, we assume that we've ran into a local minimum.
   */
  const double LocalMinimumThreshold = 0.1;

  _criticalWorker.clear();

  if (
    commandFromMaster>=peano::parallel::loadbalancing::ForkOnce ||
    commandFromMaster==peano::parallel::loadbalancing::Continue
  ) {
    double maximumWeight = getMaximumWeightOfWorkers();
    double minimumWeight = getMinimumWeightOfWorkers();

    assertion2( maximumWeight>=1.0, maximumWeight, minimumWeight );    assertion2( minimumWeight>=1.0, maximumWeight, minimumWeight );

    if (maximumWeight <= (1.0+LocalMinimumThreshold)*minimumWeight) {
      maximumWeight = 0.0;
      logDebug(
        "receivedStartCommand(int)",
        "identified local minimal, fork all workers"
      );
    }

    for ( std::map<int,double>::const_iterator p=_weightMap.begin(); p!=_weightMap.end(); p++ ) {
      if ( p->second >= (1.0-CritialPathThreshold) * maximumWeight ) {
        _criticalWorker.insert(p->first);
      }
    }
  }
  else if(
    commandFromMaster==peano::parallel::loadbalancing::ForkAllChildrenAndBecomeAdministrativeRank  ||
    commandFromMaster==peano::parallel::loadbalancing::ForkGreedy
  ) {
    for ( std::map<int,double>::const_iterator p=_weightMap.begin(); p!=_weightMap.end(); p++ ) {
      _criticalWorker.insert( p->first );
    }
  }
}


void mpibalancing::StaticBalancing::computeMaxForksOnCriticalWorker( const int commandFromMaster ) {
  if ( commandFromMaster==peano::parallel::loadbalancing::ForkGreedy ) {
    _maxForksOnCriticalWorker = peano::parallel::loadbalancing::ForkGreedy;
  }
  else if ( commandFromMaster>=peano::parallel::loadbalancing::ForkOnce ) {
    _maxForksOnCriticalWorker = getMinimumWeightOfWorkers() < std::numeric_limits<double>::max() ?
      static_cast<int>(std::ceil(
        commandFromMaster*( 1.0-getMinimumWeightOfWorkers()/getMaximumWeightOfWorkers() )
      )) : 0;
    assertion4(
      _maxForksOnCriticalWorker>=0,
      commandFromMaster, _maxForksOnCriticalWorker,
      getMinimumWeightOfWorkers(), getMaximumWeightOfWorkers()
    );

    if ( _maxForksOnCriticalWorker>static_cast<int>(commandFromMaster) ) {
      _maxForksOnCriticalWorker = static_cast<int>(commandFromMaster);
      logDebug( "receivedStartCommand(LoadBalancingFlag)", "manually reduced forks to " << _maxForksOnCriticalWorker << " due to restriction from master" );
    }

    logInfo(
      "receivedStartCommand(int)",
      _maxForksOnCriticalWorker << " forks should be done next on " << _criticalWorker.size() <<
      " worker(s) given the master fork command " << peano::parallel::loadbalancing::convertLoadBalancingFlagToString(commandFromMaster) <<
      ". Number of weight entries (workers+1)=" << _weightMap.size() <<
      ". Fork has failed before and vetos forks=" << _forkHasFailed <<
      ". Load balancing is activated=" << peano::parallel::loadbalancing::Oracle::getInstance().isLoadBalancingActivated()
    );
  }
}


void mpibalancing::StaticBalancing::receivedStartCommand( const int commandFromMaster ) {
  logTraceInWith1Argument("receivedStartCommand(LoadBalancingFlag)", peano::parallel::loadbalancing::convertLoadBalancingFlagToString(commandFromMaster));

  identifyCriticalPathes( commandFromMaster );
  computeMaxForksOnCriticalWorker( commandFromMaster );

  _weightMap[tarch::parallel::Node::getInstance().getRank()] = 1.0;

  logTraceOut("receivedStartCommand(LoadBalancingFlag)" );
}


int mpibalancing::StaticBalancing::getCommandForWorker( int workerRank, bool forkIsAllowed, bool joinIsAllowed ) {
  logTraceInWith4Arguments( "getCommandForWorker(int,bool)", workerRank, forkIsAllowed, joinIsAllowed, _joinsAllowed );
  
  int result = peano::parallel::loadbalancing::Continue;
  if (
    tarch::parallel::Node::getInstance().isGlobalMaster()
  ) {
    result = peano::parallel::loadbalancing::ForkAllChildrenAndBecomeAdministrativeRank;
  }
  else {
    if (_joinsAllowed && _workerCouldNotEraseDueToDecomposition[workerRank] && joinIsAllowed) {
      _forkHasFailed = false;
      result         = peano::parallel::loadbalancing::Join;
      logInfo(
        "forkFailed()",
        "reset fork-has-failed flag as child command is " << peano::parallel::loadbalancing::convertLoadBalancingFlagToString(result)
      );
    }
    else if (_joinsAllowed && _workerCouldNotEraseDueToDecomposition[workerRank] && !joinIsAllowed) {
      _forkHasFailed = false;
      result         = peano::parallel::loadbalancing::ContinueButTryToJoinWorkers;
      logInfo(
        "forkFailed()",
        "reset fork-has-failed flag as child command is " << peano::parallel::loadbalancing::convertLoadBalancingFlagToString(result)
      );
    }
    else if (
      _criticalWorker.count(workerRank)>0 &&
      forkIsAllowed &&
      !_forkHasFailed &&
      _maxForksOnCriticalWorker==0
    ) {
      result = peano::parallel::loadbalancing::Continue;
      assertion( result!=peano::parallel::loadbalancing::UndefinedLoadBalancingFlag );
    }
    else if (
      _criticalWorker.count(workerRank)>0 &&
      forkIsAllowed &&
      !_forkHasFailed
    ) {
      result = _maxForksOnCriticalWorker;
      assertion( result!=peano::parallel::loadbalancing::UndefinedLoadBalancingFlag );
    }
  }

  assertion( result!=peano::parallel::loadbalancing::UndefinedLoadBalancingFlag );
  logTraceOutWith1Argument( "getCommandForWorker(int,bool)", result );
  return result;
}


void mpibalancing::StaticBalancing::receivedTerminateCommand(
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
  _workerCouldNotEraseDueToDecomposition[workerRank] = workerCouldNotEraseDueToDecomposition;
  _weightMap[workerRank]                             = workerNumberOfInnerCells > 0.0 ? workerNumberOfInnerCells : 1.0;

/*
  if (parentCellLocalWorkload>_weightMap[tarch::parallel::Node::getInstance().getRank()]) {
    _weightMap[tarch::parallel::Node::getInstance().getRank()] = parentCellLocalWorkload;
  }
*/
}


void mpibalancing::StaticBalancing::plotStatistics() {
}


peano::parallel::loadbalancing::OracleForOnePhase* mpibalancing::StaticBalancing::createNewOracle(int adapterNumber) const {
  return new StaticBalancing(_joinsAllowed, _coarsestRegularInnerAndOuterGridLevel);
}


void mpibalancing::StaticBalancing::forkFailed() {
  if (!_forkHasFailed) {
    logInfo(
      "forkFailed()",
      "oracle was informed that fork has failed. No further fork attempts in this iteration"
    );
  }
  _forkHasFailed = true;
}


int mpibalancing::StaticBalancing::getCoarsestRegularInnerAndOuterGridLevel() const {
  return _coarsestRegularInnerAndOuterGridLevel;
}
