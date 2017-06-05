#include "mpibalancing/SFCDiffusionNodePoolStrategy.h"
#include "tarch/parallel/Node.h"
#include "tarch/compiler/CompilerSpecificSettings.h"

#include <sstream>
#include <limits>


tarch::logging::Log mpibalancing::SFCDiffusionNodePoolStrategy::_log( "mpibalancing::SFCDiffusionNodePoolStrategy" );


mpibalancing::SFCDiffusionNodePoolStrategy::SFCDiffusionNodePoolStrategy(int mpiRanksPerNode, int primaryMPIRanksPerNode, double waitTimeOutSec):
  NodePoolStrategy(),
  _tag(-1),
  _nodes(),
  _waitTimeOut(waitTimeOutSec),
  _mpiRanksPerNode(mpiRanksPerNode),
  _primaryMPIRanksPerNode(primaryMPIRanksPerNode),
  _nodePoolState(NodePoolState::DeployingIdlePrimaryNodes) {

  assertion( mpiRanksPerNode>0 );
  assertion( primaryMPIRanksPerNode<=mpiRanksPerNode );
}


mpibalancing::SFCDiffusionNodePoolStrategy::~SFCDiffusionNodePoolStrategy() {
}


void mpibalancing::SFCDiffusionNodePoolStrategy::fillWorkerRequestQueue(RequestQueue& queue) {
  if (queue.empty() && _nodePoolState==NodePoolState::DeployingAlsoSecondaryNodes) {
    logInfo(
      "fillWorkerRequestQueue(RequestQueue)",
      "running out of ranks. Answered all pending MPI questions before but new requests keep on dropping in. Stop to deliver MPI ranks"
    );
    _nodePoolState = NodePoolState::NoNodesLeft;
  }

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

  int totalNumberOfRequestedWorkers = 0;
  for (auto m: queue) {
    totalNumberOfRequestedWorkers += m.getNumberOfRequestedWorkers();
  }
  if (
    totalNumberOfRequestedWorkers>getNumberOfIdlePrimaryNodes()
    &&
    _nodePoolState!=NodePoolState::NoNodesLeft
  ) {
    _nodePoolState = NodePoolState::DeployingAlsoSecondaryNodes;
    logInfo(
      "fillWorkerRequestQueue(RequestQueue)",
      "have " << totalNumberOfRequestedWorkers <<
      " worker requests but only " << getNumberOfIdlePrimaryNodes() <<
      " primary node(s), i.e. code is running out of idle nodes. Start to deploy secondary nodes"
    );
    queue = sortRequestQueue( queue );
  }
}


mpibalancing::SFCDiffusionNodePoolStrategy::RequestQueue mpibalancing::SFCDiffusionNodePoolStrategy::sortRequestQueue( const RequestQueue&  queue ) {
  assertion( _nodePoolState == NodePoolState::DeployingAlsoSecondaryNodes );
  RequestQueue result;

  // Compute priorities of the nodes which equals the total number of
  // requests issued through all ranks of a node.
  std::map<std::string, int>  nodeToPriorityMap;
  for (auto p: queue) {
    assertion( p.getSenderRank()>=0 );
    assertion( p.getSenderRank()<static_cast<int>(_nodes.size()) );
    const std::string nodeName = _nodes[p.getSenderRank()].getNodeName();

    if (nodeToPriorityMap.count(nodeName)) {
      nodeToPriorityMap.insert( std::pair<std::string,int>(nodeName,0) );
    }

    nodeToPriorityMap[nodeName] += p.getNumberOfRequestedWorkers();
    assertion( nodeToPriorityMap[nodeName]>0 );
  }
  logDebug(
    "sortRequestQueue(Queue)",
    "node priority map has " << nodeToPriorityMap.size() << " entries. First entry=" <<
    nodeToPriorityMap.begin()->first << " with priority " << nodeToPriorityMap.begin()->second
  );


  // Sort the queue's requests according to their priorities.
  while ( !nodeToPriorityMap.empty() ) {
    // Find the one node with the maximum weight
    std::string maxNode     = "undef";
    int         maxPriority = -1;
    for (auto p: nodeToPriorityMap) {
      if (p.second>maxPriority) {
        maxPriority = p.second;
        maxNode     = p.first;
      }
    }
    logDebug( "sortRequestQueue(Queue)", "extract data of node " << maxNode );
    assertionEquals(nodeToPriorityMap.count(maxNode),1);
    nodeToPriorityMap.erase(maxNode);

    // First sort in all those nodes that are inside
    for(auto p: queue) {
      if (
        _nodes[p.getSenderRank()].getNodeName().compare(maxNode)==0
        &&
        !isFirstOrLastPrimaryMPIRankOnANode(p.getSenderRank())
      ) {
        logDebug( "sortRequestQueue(Queue)", "prioritise rank " << p.getSenderRank() );
        result.push_back(p);
      }
    }

    // Now insert all those guys that are outside along the SFC
    for(auto p: queue) {
      if (
        _nodes[p.getSenderRank()].getNodeName().compare(maxNode)==0
        &&
        isFirstOrLastPrimaryMPIRankOnANode(p.getSenderRank())
      ) {
        result.push_back(p);
      }
    }
  }

  assertionEquals( result.size(), queue.size() );
  return result;
}


tarch::parallel::messages::WorkerRequestMessage mpibalancing::SFCDiffusionNodePoolStrategy::extractElementFromRequestQueue(RequestQueue& queue) {
  assertion( !queue.empty() );

  RequestQueue::iterator pResultInQueue = queue.begin();
  tarch::parallel::messages::WorkerRequestMessage result = *pResultInQueue;
  queue.erase(pResultInQueue);

  return result;
}


bool mpibalancing::SFCDiffusionNodePoolStrategy::isFirstOrLastPrimaryMPIRankOnANode(int rank) const {
  return (rank % _mpiRanksPerNode == 0)
      || (rank % _mpiRanksPerNode == _primaryMPIRanksPerNode-1);
}


bool mpibalancing::SFCDiffusionNodePoolStrategy::isPrimaryMPIRank(int rank) const {
  return rank % _mpiRanksPerNode < _primaryMPIRanksPerNode;
}


void mpibalancing::SFCDiffusionNodePoolStrategy::addNode(const tarch::parallel::messages::RegisterAtNodePoolMessage& node) {
  #ifdef Parallel
  logTraceInWith1Argument( "addNode(...)", node.getSenderRank() );

  while (static_cast<int>(_nodes.size())<=node.getSenderRank()) {
    _nodes.push_back( NodePoolListEntry() );
  }

  const std::string name = tarch::parallel::StringTools::convert(node.getNodeName());
  assertion(_primaryMPIRanksPerNode!=0);
  _nodes[node.getSenderRank()] = NodePoolListEntry(
    name,
    isPrimaryMPIRank(node.getSenderRank())
  );

  logTraceOutWith1Argument( "addNode(...)", newEntry.toString() );
  #endif
}


void mpibalancing::SFCDiffusionNodePoolStrategy::removeNode( int rank ) {
  assertion( isRegisteredNode(rank) );
  _nodes[rank] = NodePoolListEntry(); // overwrite with invalidated entry
}


int mpibalancing::SFCDiffusionNodePoolStrategy::getNumberOfIdlePrimaryNodes() const {
  int result = 0;
  for (auto node: _nodes) {
    if (node.isIdlePrimaryNode()) {
      result++;
    }
  }
  return result;
}


bool mpibalancing::SFCDiffusionNodePoolStrategy::hasIdleNode(int forMaster) const {
  switch (_nodePoolState) {
    case NodePoolState::DeployingIdlePrimaryNodes:
      {
        for (auto node: _nodes) {
          if (node.isIdlePrimaryNode()) {
            return true;
          }
        }
      }
      break;
    case NodePoolState::DeployingAlsoSecondaryNodes:
    {
      for (auto node: _nodes) {
        if (node.isIdlePrimaryNode() || node.isIdleSecondaryNode()) {
          return true;
        }
      }
    }
    break;
    case NodePoolState::NoNodesLeft:
      break;
  }
  return false;
}


int mpibalancing::SFCDiffusionNodePoolStrategy::getNumberOfIdleNodes() const {
  int result = 0;
  for (auto node: _nodes) {
    if (node.isIdlePrimaryNode() || node.isIdleSecondaryNode()) {
      result++;
    }
  }
  return result;
}


void mpibalancing::SFCDiffusionNodePoolStrategy::setNodeIdle( int rank ) {
  assertion( isRegisteredNode(rank) );
  _nodes[rank].deActivate();
  if (isPrimaryMPIRank(rank)) {
    _nodePoolState = NodePoolState::DeployingIdlePrimaryNodes;
  }
  else if (_nodePoolState == NodePoolState::NoNodesLeft) {
    _nodePoolState = NodePoolState::DeployingAlsoSecondaryNodes;
  }
}


bool mpibalancing::SFCDiffusionNodePoolStrategy::isRegisteredNode(int rank) const {
  return static_cast<int>(_nodes.size())>rank
      && _nodes[rank].isRegistered();
}


bool mpibalancing::SFCDiffusionNodePoolStrategy::isIdleNode(int rank) const {
  return (static_cast<int>(_nodes.size())>rank)
      && (_nodes[rank].isIdlePrimaryNode() || _nodes[rank].isIdleSecondaryNode());
}


int mpibalancing::SFCDiffusionNodePoolStrategy::getNumberOfRegisteredNodes() const {
  int result = 0;
  for (auto node: _nodes) {
    if (node.isRegistered()) {
      result++;
    }
  }
  return result;
}


std::string mpibalancing::SFCDiffusionNodePoolStrategy::toString() const {
  std::ostringstream result;
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    p->toString(result);
  }
  return result.str();
}


void mpibalancing::SFCDiffusionNodePoolStrategy::reserveParticularNode(int rank) {
  assertion( isRegisteredNode(rank) );
  assertion( _nodes[rank].isIdlePrimaryNode() || _nodes[rank].isIdleSecondaryNode() );

  _nodes[rank].activate();
}


int mpibalancing::SFCDiffusionNodePoolStrategy::reserveNode(int forMaster) {
  switch (_nodePoolState) {
    case NodePoolState::DeployingIdlePrimaryNodes:
      {
        for (int i=0; i<static_cast<int>(_nodes.size());i++) {
          if (_nodes[i].isIdlePrimaryNode()) {
            _nodes[i].activate();
            return i;
          }
        }
      }
      break;
    case NodePoolState::DeployingAlsoSecondaryNodes:
      {
        assertion( hasIdleNode(forMaster) );
        assertion( !_nodes[forMaster].isIdlePrimaryNode() );
        assertion( !_nodes[forMaster].isIdleSecondaryNode() );

        for (int i=1; i>0; i++) { // not an endless loop
          if (
            forMaster+i < tarch::parallel::Node::getInstance().getNumberOfNodes()
            &&
            (_nodes[forMaster+i].isIdlePrimaryNode() || _nodes[forMaster+i].isIdleSecondaryNode())
          ) {
            _nodes[forMaster+i].activate();
            return forMaster+i;
          }
          if (
            forMaster-i >= 0
            &&
            ( _nodes[forMaster-i].isIdlePrimaryNode() || _nodes[forMaster-i].isIdleSecondaryNode())
          ) {
            _nodes[forMaster-i].activate();
            return forMaster-i;
          }
        }
      }
      break;
    case NodePoolState::NoNodesLeft:
      break;
  }
  assertion(false);
  return -1;
}


void mpibalancing::SFCDiffusionNodePoolStrategy::setNodePoolTag(int tag) {
  _tag = tag;
}


mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::NodePoolListEntry():
  _state(State::Undef),
  _name("undef") {
}


mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::NodePoolListEntry(
  const std::string& name,
  bool isPrimaryNode
  ):
  _state( isPrimaryNode ? State::WorkingPrimaryNode : State::WorkingSecondaryNode ),
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
  assertion( _state==State::IdlePrimaryNode || _state==State::IdleSecondaryNode );
  _state = _state==State::IdlePrimaryNode ? State::WorkingPrimaryNode : State::WorkingSecondaryNode;
}


void mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::deActivate() {
  assertion( _state==State::WorkingPrimaryNode || _state==State::WorkingSecondaryNode );
  _state = _state==State::WorkingPrimaryNode ? State::IdlePrimaryNode : State::IdleSecondaryNode;
}


bool mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::isRegistered() const {
  return _state!=State::Undef;
}


bool mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::isIdlePrimaryNode() const {
  return _state==State::IdlePrimaryNode;
}


bool mpibalancing::SFCDiffusionNodePoolStrategy::NodePoolListEntry::isIdleSecondaryNode() const {
  return _state==State::IdleSecondaryNode;
}
