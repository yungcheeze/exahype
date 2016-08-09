#include "HotspotBalancing.h"
#include "tarch/Assertions.h"
#include "tarch/la/ScalarOperations.h"
#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "peano/parallel/loadbalancing/Oracle.h"


tarch::logging::Log mpibalancing::HotspotBalancing::_log( "mpibalancing::HotspotBalancing" );


bool mpibalancing::HotspotBalancing::_forkHasFailed = false;


mpibalancing::HotspotBalancing::HotspotBalancing(bool joinsAllowed, int coarsestRegularInnerAndOuterGridLevel):
  _coarsestRegularInnerAndOuterGridLevel(coarsestRegularInnerAndOuterGridLevel),
  _joinsAllowed(joinsAllowed),
  _criticalWorker(),
  _maxForksOnCriticalWorker(THREE_POWER_D) {
  _weightMap.insert( std::pair<int,double>(tarch::parallel::Node::getInstance().getRank(), 1.0) );
  _workerCouldNotEraseDueToDecomposition.insert( std::pair<int,bool>(tarch::parallel::Node::getInstance().getRank(), false) );
}


mpibalancing::HotspotBalancing::~HotspotBalancing() {
}


double mpibalancing::HotspotBalancing::getMaximumWeightOfWorkers() const {
  double maximumWeight = std::numeric_limits<double>::min();
  for ( std::map<int,double>::const_iterator p=_weightMap.begin(); p!=_weightMap.end(); p++ ) {
    if ( p->second > maximumWeight ) {
      maximumWeight = p->second;
    }
  }

  assertion1( maximumWeight>=1.0, maximumWeight );

  return maximumWeight;
}


double mpibalancing::HotspotBalancing::getMinimumWeightOfWorkers() const {
  double minimumWeight  = std::numeric_limits<double>::max();
  for ( std::map<int,double>::const_iterator p=_weightMap.begin(); p!=_weightMap.end(); p++ ) {
    if ( p->second < minimumWeight ) {
      minimumWeight = p->second;
    }
  }
  assertion1( minimumWeight>=1.0, minimumWeight );
  return minimumWeight;
}


void mpibalancing::HotspotBalancing::identifyCriticalPathes( peano::parallel::loadbalancing::LoadBalancingFlag commandFromMaster ) {
  /**
   * We consider the CritialPathThreshold upper percent of the workers to be
   * critical.
   */
  const double CritialPathThreshold = 0.1;

  _criticalWorker.clear();

  if (
    commandFromMaster>=peano::parallel::loadbalancing::LoadBalancingFlag::ForkOnce ||
    commandFromMaster==peano::parallel::loadbalancing::LoadBalancingFlag::Continue
  ) {
    double maximumWeight = getMaximumWeightOfWorkers();
    double minimumWeight = getMinimumWeightOfWorkers();
    assertion(minimumWeight<=maximumWeight);

    assertion2( maximumWeight>=1.0, maximumWeight, minimumWeight );    assertion2( minimumWeight>=1.0, maximumWeight, minimumWeight );

    for ( std::map<int,double>::const_iterator p=_weightMap.begin(); p!=_weightMap.end(); p++ ) {
      if ( p->second >= (1.0-CritialPathThreshold) * maximumWeight ) {
        _criticalWorker.insert(p->first);
      }
    }

    if (!_weightMap.empty()) {
      std::ostringstream msg;
      msg << "min weight=" << minimumWeight << ", max weight=" << maximumWeight << ". Critical workers are/is";
      for ( auto p: _weightMap ) {
        if (p!=*_weightMap.begin()) {
          msg << ",";
        }
        msg << " (" << p.first << "," << p.second << ")";
      }
      logInfo( "receivedStartCommand(int)", msg.str() );
    }
  }
}


void mpibalancing::HotspotBalancing::computeMaxForksOnCriticalWorker( peano::parallel::loadbalancing::LoadBalancingFlag commandFromMaster ) {
  if ( commandFromMaster>=peano::parallel::loadbalancing::LoadBalancingFlag::ForkOnce ) {
    _maxForksOnCriticalWorker = getMinimumWeightOfWorkers() < std::numeric_limits<double>::max() ?
      static_cast<int>(std::ceil(
        static_cast<int>(commandFromMaster)*( 1.0-getMinimumWeightOfWorkers()/getMaximumWeightOfWorkers() )
      )) : 0;

    if ( _maxForksOnCriticalWorker>static_cast<int>(commandFromMaster) ) {
      _maxForksOnCriticalWorker = static_cast<int>(commandFromMaster);
      logDebug( "receivedStartCommand(LoadBalancingFlag)", "manually reduced forks to " << _maxForksOnCriticalWorker << " due to restriction from master" );
    }
    else if (_maxForksOnCriticalWorker>static_cast<int>(peano::parallel::loadbalancing::LoadBalancingFlag::ForkGreedy)) {
      _maxForksOnCriticalWorker = static_cast<int>(peano::parallel::loadbalancing::LoadBalancingFlag::ForkGreedy)/2;
    }

    logInfo(
      "receivedStartCommand(int)",
      _maxForksOnCriticalWorker << " forks should be done next on " << _criticalWorker.size()-1 <<
      " worker(s) plus the actual rank given the master fork command " << peano::parallel::loadbalancing::convertLoadBalancingFlagToString(commandFromMaster) <<
      ". Number of weight entries (workers+1)=" << _weightMap.size() <<
      ". Fork has failed before and vetos forks=" << _forkHasFailed <<
      ". Load balancing is activated=" << peano::parallel::loadbalancing::Oracle::getInstance().isLoadBalancingActivated()
    );
  }
}


void mpibalancing::HotspotBalancing::receivedStartCommand( peano::parallel::loadbalancing::LoadBalancingFlag commandFromMaster ) {
  logTraceInWith1Argument("receivedStartCommand(LoadBalancingFlag)", peano::parallel::loadbalancing::convertLoadBalancingFlagToString(commandFromMaster));

  identifyCriticalPathes( commandFromMaster );
  computeMaxForksOnCriticalWorker( commandFromMaster );

  _weightMap[tarch::parallel::Node::getInstance().getRank()] = 1.0;

  logTraceOut("receivedStartCommand(LoadBalancingFlag)" );
}


peano::parallel::loadbalancing::LoadBalancingFlag  mpibalancing::HotspotBalancing::getCommandForWorker( int workerRank, bool forkIsAllowed, bool joinIsAllowed ) {
  logTraceInWith4Arguments( "getCommandForWorker(int,bool)", workerRank, forkIsAllowed, joinIsAllowed, _joinsAllowed );
  
  peano::parallel::loadbalancing::LoadBalancingFlag  result = peano::parallel::loadbalancing::LoadBalancingFlag::Continue;

  if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    result = peano::parallel::loadbalancing::LoadBalancingFlag::ForkAllChildrenAndBecomeAdministrativeRank;
  }
  else if (_joinsAllowed && _workerCouldNotEraseDueToDecomposition[workerRank] && joinIsAllowed) {
    _forkHasFailed = false;
    result         = peano::parallel::loadbalancing::LoadBalancingFlag::Join;
  }
  else if (_joinsAllowed && _workerCouldNotEraseDueToDecomposition[workerRank] && !joinIsAllowed) {
    _forkHasFailed = false;
    result         = peano::parallel::loadbalancing::LoadBalancingFlag::ContinueButTryToJoinWorkers;
  }
  else if (
    _criticalWorker.count(workerRank)>0 &&
    forkIsAllowed &&
    !_forkHasFailed &&
    _maxForksOnCriticalWorker==0
  ) {
    result = peano::parallel::loadbalancing::LoadBalancingFlag::Continue;
    assertion( result!=peano::parallel::loadbalancing::LoadBalancingFlag::UndefinedLoadBalancingFlag );
  }
  else if (
    _criticalWorker.count(workerRank)>0 &&
    forkIsAllowed &&
    !_forkHasFailed
  ) {
    result = static_cast<peano::parallel::loadbalancing::LoadBalancingFlag>(_maxForksOnCriticalWorker);
    assertion( result!=peano::parallel::loadbalancing::LoadBalancingFlag::UndefinedLoadBalancingFlag );
  }

  // @todo Raus
  logInfo(
    "getCommandForWorker()",
    "send rank " << workerRank << " flag " << peano::parallel::loadbalancing::convertLoadBalancingFlagToString(result)
  );


  assertion( result!=peano::parallel::loadbalancing::LoadBalancingFlag::UndefinedLoadBalancingFlag );
  logTraceOutWith1Argument( "getCommandForWorker(int,bool)", result );
  return result;
}


/*
void mpibalancing::HotspotBalancing::receivedTerminateCommand(
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
}
*/


void mpibalancing::HotspotBalancing::plotStatistics() {
}


peano::parallel::loadbalancing::OracleForOnePhase* mpibalancing::HotspotBalancing::createNewOracle(int adapterNumber) const {
  return new HotspotBalancing(_joinsAllowed, _coarsestRegularInnerAndOuterGridLevel);
}


void mpibalancing::HotspotBalancing::forkFailed() {
  if (!_forkHasFailed) {
    logInfo(
      "forkFailed()",
      "oracle was informed that fork has failed. No further fork attempts in this iteration"
    );
  }
  _forkHasFailed = true;
}


int mpibalancing::HotspotBalancing::getCoarsestRegularInnerAndOuterGridLevel() const {
  return _coarsestRegularInnerAndOuterGridLevel;
}
