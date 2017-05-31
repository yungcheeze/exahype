#include "mpibalancing/SFCDiffusionNodePoolStrategy.h"
#include "tarch/parallel/Node.h"
#include "tarch/compiler/CompilerSpecificSettings.h"

#include <sstream>
#include <limits>


tarch::logging::Log mpibalancing::SFCDiffusionNodePoolStrategy::_log( "mpibalancing::SFCDiffusionNodePoolStrategy" );


mpibalancing::SFCDiffusionNodePoolStrategy::SFCDiffusionNodePoolStrategy(int mpiRanksPerNode, double waitTimeOutSec):
  NodePoolStrategy(),
  _tag(-1),
  _nodes(),
  _waitTimeOut(waitTimeOutSec),
  _ranksPerNode(mpiRanksPerNode) {

  assertion(_ranksPerNode>0);
}


mpibalancing::SFCDiffusionNodePoolStrategy::~SFCDiffusionNodePoolStrategy() {
}


void mpibalancing::SFCDiffusionNodePoolStrategy::fillWorkerRequestQueue(RequestQueue& queue) {
  #ifdef Parallel
  assertion( _tag >= 0 );

  const std::clock_t waitTimeoutTimeStamp = clock() + static_cast<std::clock_t>(std::floor(_waitTimeOut * CLOCKS_PER_SEC));

  bool continueToWait = true;
  while ( continueToWait ) {
    while (tarch::parallel::messages::WorkerRequestMessage::isMessageInQueue(_tag, true)) {
      tarch::parallel::messages::WorkerRequestMessage message;
      message.receive(MPI_ANY_SOURCE,_tag, true, SendAndReceiveLoadBalancingMessagesBlocking);
      queue.push_back( message );
    }

    continueToWait =
        hasIdleNode(-1) &&
        (static_cast<int>(queue.size())+1 < getNumberOfRegisteredNodes()-getNumberOfIdleNodes()) &&
        (clock() < waitTimeoutTimeStamp);
  }
  #endif
}



tarch::parallel::messages::WorkerRequestMessage mpibalancing::SFCDiffusionNodePoolStrategy::extractElementFromRequestQueue(RequestQueue& queue) {
  assertion( !queue.empty() );

  RequestQueue::iterator pResultInQueue = queue.begin();
  tarch::parallel::messages::WorkerRequestMessage result = *pResultInQueue;
  queue.erase(pResultInQueue);

  return result;
}



void mpibalancing::SFCDiffusionNodePoolStrategy::addNode(const tarch::parallel::messages::RegisterAtNodePoolMessage& node) {
  #ifdef Parallel
  logTraceInWith1Argument( "addNode(...)", node.getSenderRank() );

  while (static_cast<int>(_nodes.size())<node.getSenderRank()) {
    _nodes.push_back( NodePoolListEntry() );
  }

  _nodes[node.getSenderRank()] = NodePoolListEntry(
    tarch::parallel::StringTools::convert(node.getNodeName())
  );

  logTraceOutWith1Argument( "addNode(...)", newEntry.toString() );
  #endif
}

void mpibalancing::SFCDiffusionNodePoolStrategy::removeNode( int rank ) {
  assertion( isRegisteredNode(rank) );
  assertion(static_cast<int>(_nodes.size())<rank);
  // @todo
//  _nodes[rank].kill();

}


bool mpibalancing::SFCDiffusionNodePoolStrategy::hasIdleNode(int forMaster) const {
/*
  if (_nodes.empty()) return false;

  if (forMaster==AnyMaster) {
    return _nodes.front().isIdle() || _nodes.back().isIdle();
  }
  else {
    return _nodes.front().isIdle();
  }
*/
  // @todo
  return false;
}


int mpibalancing::SFCDiffusionNodePoolStrategy::removeNextIdleNode() {
  // @todo
  return 0;
}


int mpibalancing::SFCDiffusionNodePoolStrategy::getNumberOfIdleNodes() const {
/*
  int result = 0;
  NodeContainer::const_iterator p = _nodes.begin();
  while (p != _nodes.end()&& p->isIdle() ) {
    p++;
    result++;
  }
  return result;
*/
  // @todo
  return 0;
}


void mpibalancing::SFCDiffusionNodePoolStrategy::setNodeIdle( int rank ) {
/*
  for (
    NodeContainer::iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      p->deActivate();
    }
  }

  _nodes.sort();
*/
  // @todo
}


bool mpibalancing::SFCDiffusionNodePoolStrategy::isRegisteredNode(int rank) const {
/*
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      return true;
    }
  }
*/
  // @todo
  return false;
}


bool mpibalancing::SFCDiffusionNodePoolStrategy::isIdleNode(int rank) const {
/*
  assertion1( isRegisteredNode(rank), rank );
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank && p->isIdle() ) {
      return true;
    }
  }
*/
  // @todo
  return false;
}


int mpibalancing::SFCDiffusionNodePoolStrategy::getNumberOfRegisteredNodes() const {
  // @todo
  return 0;
//  return static_cast<int>( _nodes.size() );
}


std::string mpibalancing::SFCDiffusionNodePoolStrategy::toString() const {
  std::ostringstream result;
/*
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    p->toString(result);
  }
*/
  // @todo
  return result.str();
}


void mpibalancing::SFCDiffusionNodePoolStrategy::reserveParticularNode(int rank) {
/*
  assertion( rank>=0 );
  assertion1( isIdleNode(rank), rank );

  for (
    NodeContainer::iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      p->activate();
    }
  }

  _nodes.sort();
*/
  // @todo
}


int mpibalancing::SFCDiffusionNodePoolStrategy::reserveNode(int forMaster) {
/*
  assertion1(hasIdleNode(forMaster),forMaster);

  NodePoolListEntry result;

  if (_nodes.front().isIdle()) {
    result = _nodes.front();
    _nodes.pop_front();
  }
  else {
    result = _nodes.back();
    _nodes.pop_back();
  }

  result.activate();
  _nodes.push_back(result);
  _nodes.sort();

  for (NodeContainer::iterator p=_nodes.begin(); p!=_nodes.end(); p++ ) {
    if (p->getRank()==forMaster) {
      p->addNewWorker();
    }
  }

  updateNodeWeights();

  return result.getRank();
*/
  // @todo
  return 0;
}


void mpibalancing::SFCDiffusionNodePoolStrategy::setNodePoolTag(int tag) {
  _tag = tag;
}



mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::NodePoolListEntry():
  _state(State::Undef),
  _name("undef") {
}


mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::NodePoolListEntry(
  const std::string& name
  ):
  _state(State::Undef),
  _name(name) {
}


mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::~NodePoolListEntry() {
}


std::string mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::getNodeName() const {
  return _name;
}


std::string mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::toString() const {
  std::ostringstream out;
  toString(out);
  return out.str();
}


void mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::toString(std::ostream& out) const {
  out << "(state:";
  switch (_state) {
    case State::Undef:                 out << "undef";                   break;
    case State::WorkingPrimaryNode:    out << "working-primary-node";    break;
    case State::WorkingSecondaryNode:  out << "working-secondary-node";  break;
    case State::IdlePrimaryNode:       out << "idle-primary-node";       break;
    case State::IdleSecondaryNode:     out << "idle-secondary-node";     break;
  }
  out << ",name:" << _name << ")";
}


void mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::activate() {
/*
  assertionEquals1( _state, IDLE, toString() );
  _state = State.Working;
*/
}


void mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::deActivate() {
//  _state         = IDLE;
}
/*


bool mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::isIdle() const {
//  return _state != State.Working;
  return false;
}

*/
