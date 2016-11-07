#include "HotspotBalancing.h"
#include "tarch/Assertions.h"
#include "tarch/la/ScalarOperations.h"
#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "peano/parallel/loadbalancing/Oracle.h"


tarch::logging::Log  mpibalancing::HotspotBalancing::_log( "mpibalancing::HotspotBalancing" );
int                  mpibalancing::HotspotBalancing::_loadBalancingTag = -1;


bool                        mpibalancing::HotspotBalancing::_forkHasFailed = false;
std::map<int,double>        mpibalancing::HotspotBalancing::_weightMap;
std::map<int,bool>          mpibalancing::HotspotBalancing::_workerCouldNotEraseDueToDecomposition;
int                         mpibalancing::HotspotBalancing::_regularLevelAlongBoundary = 0;


mpibalancing::HotspotBalancing::HotspotBalancing(bool joinsAllowed, int coarsestRegularInnerAndOuterGridLevel):
  _joinsAllowed(joinsAllowed),
  _criticalWorker(),
  _maxForksOnCriticalWorker(THREE_POWER_D) {
  _workerCouldNotEraseDueToDecomposition.insert( std::pair<int,bool>(tarch::parallel::Node::getInstance().getRank(), false) );
  _regularLevelAlongBoundary = coarsestRegularInnerAndOuterGridLevel;

  if (_loadBalancingTag<0) {
    _loadBalancingTag = tarch::parallel::Node::reserveFreeTag("mpibalancing::HotspotBalancing");
    assertion(_loadBalancingTag>=0);
  }
}


mpibalancing::HotspotBalancing::~HotspotBalancing() {
}


double mpibalancing::HotspotBalancing::getMaximumWeightOfWorkers() {
  double maximumWeight = std::numeric_limits<double>::min();
  for ( std::map<int,double>::const_iterator p=_weightMap.begin(); p!=_weightMap.end(); p++ ) {
    if ( p->second > maximumWeight ) {
      maximumWeight = p->second;
    }
  }
  return maximumWeight;
}


double mpibalancing::HotspotBalancing::getMinimumWeightOfWorkers() {
  double minimumWeight  = std::numeric_limits<double>::max();
  for ( std::map<int,double>::const_iterator p=_weightMap.begin(); p!=_weightMap.end(); p++ ) {
    if ( p->second < minimumWeight ) {
      minimumWeight = p->second;
    }
  }
  return minimumWeight;
}


void mpibalancing::HotspotBalancing::identifyCriticalPathes( peano::parallel::loadbalancing::LoadBalancingFlag commandFromMaster ) {
  /**
   * We consider the CritialPathThreshold upper percent of the workers to be
   * critical.
   */
  const double CritialPathThreshold = 0.1;

  _criticalWorker.clear();
  assertion(_criticalWorker.size()==0);

  if (
    (
    commandFromMaster>=peano::parallel::loadbalancing::LoadBalancingFlag::ForkOnce ||
    commandFromMaster==peano::parallel::loadbalancing::LoadBalancingFlag::Continue
    )
    &&
    !_forkHasFailed
  ) {
    double maximumWeight = getMaximumWeightOfWorkers();

    for ( std::map<int,double>::const_iterator p=_weightMap.begin(); p!=_weightMap.end(); p++ ) {
      if ( p->second >= (1.0-CritialPathThreshold) * maximumWeight ) {
        _criticalWorker.insert(p->first);
        assertion(_criticalWorker.size()>0);
      }
    }
  }
}


void mpibalancing::HotspotBalancing::computeMaxForksOnCriticalWorker( peano::parallel::loadbalancing::LoadBalancingFlag commandFromMaster ) {
  if ( !_criticalWorker.empty() ) {
    _maxForksOnCriticalWorker =
      static_cast<int>(std::ceil(
        static_cast<int>(commandFromMaster)*( 1.0-getMinimumWeightOfWorkers()/getMaximumWeightOfWorkers() )
      ));

    if ( _maxForksOnCriticalWorker>static_cast<int>(commandFromMaster) ) {
      _maxForksOnCriticalWorker = static_cast<int>(commandFromMaster);
      logDebug( "receivedStartCommand(LoadBalancingFlag)", "manually reduced forks to " << _maxForksOnCriticalWorker << " due to restriction from master" );
    }
    else if (_maxForksOnCriticalWorker>static_cast<int>(peano::parallel::loadbalancing::LoadBalancingFlag::ForkGreedy)) {
      _maxForksOnCriticalWorker = static_cast<int>(peano::parallel::loadbalancing::LoadBalancingFlag::ForkGreedy)/2;
    }
    else if (_maxForksOnCriticalWorker<=0) {
      _maxForksOnCriticalWorker = 1;
    }

    if (_maxForksOnCriticalWorker>0 && !_forkHasFailed) {
      logInfo(
        "computeMaxForksOnCriticalWorker(LoadBalancingFlag)",
        _maxForksOnCriticalWorker << " forks should be done (-2 = continue) next on " << _criticalWorker.size() <<
        " worker(s) given the master fork command " << peano::parallel::loadbalancing::convertLoadBalancingFlagToString(commandFromMaster) <<
        ". Number of weight entries (workers+1)=" << _weightMap.size() <<
        ". Fork has failed before and vetos forks=" << _forkHasFailed <<
        ". max weight=" << getMaximumWeightOfWorkers() <<
        ". min weight=" << getMinimumWeightOfWorkers() <<
        ". Load balancing is activated=" << peano::parallel::loadbalancing::Oracle::getInstance().isLoadBalancingActivated()
      );
      std::ostringstream msg;
      msg << "critical workers are/is";
      for ( auto p: _criticalWorker ) {
        if (p!=*_criticalWorker.begin()) {
          msg << ",";
        }
        msg << " (" << p << "," << _weightMap[p] << ")";
      }
      msg << ". local weight= " << _weightMap[tarch::parallel::Node::getInstance().getRank()];
      logInfo( "identifyCriticalPathes(LoadBalancingFlag)", msg.str() );
    }
  }
}


void mpibalancing::HotspotBalancing::receivedStartCommand( peano::parallel::loadbalancing::LoadBalancingFlag commandFromMaster ) {
  logTraceInWith1Argument("receivedStartCommand(LoadBalancingFlag)", peano::parallel::loadbalancing::convertLoadBalancingFlagToString(commandFromMaster));

  identifyCriticalPathes( commandFromMaster );
  computeMaxForksOnCriticalWorker( commandFromMaster );

  logTraceOut("receivedStartCommand(LoadBalancingFlag)" );
}


peano::parallel::loadbalancing::LoadBalancingFlag  mpibalancing::HotspotBalancing::getCommandForWorker( int workerRank, bool forkIsAllowed, bool joinIsAllowed ) {
  logTraceInWith4Arguments( "getCommandForWorker(int,bool)", workerRank, forkIsAllowed, joinIsAllowed, _joinsAllowed );
  
  peano::parallel::loadbalancing::LoadBalancingFlag  result = peano::parallel::loadbalancing::LoadBalancingFlag::Continue;

  if (
    tarch::parallel::Node::getInstance().isGlobalMaster()
    &&
    forkIsAllowed
  ) {
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

  assertion( result!=peano::parallel::loadbalancing::LoadBalancingFlag::UndefinedLoadBalancingFlag );
  logTraceOutWith1Argument( "getCommandForWorker(int,bool)", static_cast<int>(result) );
  return result;
}


void mpibalancing::HotspotBalancing::mergeWithMaster(
  int     workerRank,
  bool    workerCouldNotEraseDueToDecomposition
) {
  double workerWeight;
  MPI_Recv( &workerWeight, 1, MPI_DOUBLE, workerRank, _loadBalancingTag, tarch::parallel::Node::getInstance().getCommunicator(), MPI_STATUS_IGNORE );
  _workerCouldNotEraseDueToDecomposition[workerRank] = workerCouldNotEraseDueToDecomposition;
  _weightMap[workerRank]                             = workerWeight > 1.0 ? workerWeight : 1.0;
}


void mpibalancing::HotspotBalancing::setLocalWeightAndNotifyMaster(
  double localWeight
) {
  _weightMap[tarch::parallel::Node::getInstance().getRank()] = localWeight;

  double ranksWeight = getMaximumWeightOfWorkers();

  MPI_Send( &ranksWeight, 1, MPI_DOUBLE, tarch::parallel::NodePool::getInstance().getMasterRank(), _loadBalancingTag, tarch::parallel::Node::getInstance().getCommunicator() );
}


void mpibalancing::HotspotBalancing::plotStatistics() {
}


peano::parallel::loadbalancing::OracleForOnePhase* mpibalancing::HotspotBalancing::createNewOracle(int adapterNumber) const {
  return new HotspotBalancing(_joinsAllowed, _regularLevelAlongBoundary);
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


int mpibalancing::HotspotBalancing::getRegularLevelAlongBoundary() const {
  return _regularLevelAlongBoundary;
}
