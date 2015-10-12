#include "ExaHyPE/repositories/RepositorySTDStack.h"

#include "tarch/Assertions.h"
#include "tarch/timing/Watch.h"

#include "tarch/compiler/CompilerSpecificSettings.h"

#ifdef Parallel
#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "peano/parallel/SendReceiveBufferPool.h"
#include "peano/parallel/loadbalancing/Oracle.h"
#endif

#include "peano/datatraversal/autotuning/Oracle.h"

#include "tarch/compiler/CompilerSpecificSettings.h"

#if !defined(CompilerICC)
#include "peano/grid/Grid.cpph"
#endif


tarch::logging::Log ExaHyPE::repositories::RepositorySTDStack::_log( "ExaHyPE::repositories::RepositorySTDStack" );


ExaHyPE::repositories::RepositorySTDStack::RepositorySTDStack(
  peano::geometry::Geometry&                   geometry,
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset
):
  _geometry(geometry),
  _cellStack(),
  _vertexStack(),
  _solverState(),
  _gridWithCreateGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositorySTDStack(...)" );
  _repositoryState.setAction( ExaHyPE::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(ExaHyPE::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(ExaHyPE::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositorySTDStack(...)" );
}



ExaHyPE::repositories::RepositorySTDStack::RepositorySTDStack(
  peano::geometry::Geometry&                   geometry
):
  _geometry(geometry),
  _cellStack(),
  _vertexStack(),
  _solverState(),
  _gridWithCreateGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositorySTDStack(Geometry&)" );

  _repositoryState.setAction( ExaHyPE::records::RepositoryState::Terminate );
  
  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(ExaHyPE::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(ExaHyPE::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositorySTDStack(Geometry&)" );
}
    
   
ExaHyPE::repositories::RepositorySTDStack::~RepositorySTDStack() {
  assertionMsg( _repositoryState.getAction() == ExaHyPE::records::RepositoryState::Terminate, "terminate() must be called before destroying repository." );
}


void ExaHyPE::repositories::RepositorySTDStack::restart(
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
  int                                          domainLevel,
  const tarch::la::Vector<DIMENSIONS,int>&     positionOfCentralElementWithRespectToCoarserRemoteLevel
) {
  logTraceInWith4Arguments( "restart(...)", domainSize, domainOffset, domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel );
  #ifdef Parallel
  assertion( !tarch::parallel::Node::getInstance().isGlobalMaster());
  #endif
  
  logInfo( "restart(...)", "start node for subdomain " << domainOffset << "x" << domainSize << " on level " << domainLevel );
  
  assertion( _repositoryState.getAction() == ExaHyPE::records::RepositoryState::Terminate );

  _vertexStack.clear();
  _cellStack.clear();

  _gridWithCreateGrid.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);


  _solverState.restart();

  logTraceOut( "restart(...)" );
}


void ExaHyPE::repositories::RepositorySTDStack::terminate() {
  logTraceIn( "terminate()" );
  
  _repositoryState.setAction( ExaHyPE::records::RepositoryState::Terminate );
  
  #ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
      _repositoryState,
      peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
    );
  }
  peano::parallel::SendReceiveBufferPool::getInstance().terminate();
  #endif

  _gridWithCreateGrid.terminate();


  logTraceOut( "terminate()" );
}


ExaHyPE::State& ExaHyPE::repositories::RepositorySTDStack::getState() {
  return _solverState;
}


const ExaHyPE::State& ExaHyPE::repositories::RepositorySTDStack::getState() const {
  return _solverState;
}


void ExaHyPE::repositories::RepositorySTDStack::iterate(int numberOfIterations, bool exchangeBoundaryVertices) {
  tarch::timing::Watch watch( "ExaHyPE::repositories::RepositorySTDStack", "iterate(bool)", false);
  
  #ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    _repositoryState.setNumberOfIterations(numberOfIterations);
    _repositoryState.setExchangeBoundaryVertices(exchangeBoundaryVertices);
    tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
      _repositoryState,
      peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
    );
  }
  else {
    assertionEquals( numberOfIterations, 1 );
    numberOfIterations = _repositoryState.getNumberOfIterations();
  }
  
  peano::parallel::SendReceiveBufferPool::getInstance().exchangeBoundaryVertices(_repositoryState.getExchangeBoundaryVertices());

  if ( numberOfIterations > 1 && ( peano::parallel::loadbalancing::Oracle::getInstance().isLoadBalancingActivated() || _solverState.isInvolvedInJoinOrFork() )) {
    logWarning( "iterate()", "iterate invoked for multiple traversals though load balancing is switched on or grid is not balanced globally. Use activateLoadBalancing(false) to deactivate the load balancing before" );
  }

  peano::datatraversal::autotuning::Oracle::getInstance().switchToOracle(_repositoryState.getAction());

  peano::parallel::loadbalancing::Oracle::getInstance().switchToOracle(_repositoryState.getAction());
  
  _solverState.currentlyRunsMultipleIterations(_repositoryState.getNumberOfIterations()>1);
  #else
  peano::datatraversal::autotuning::Oracle::getInstance().switchToOracle(_repositoryState.getAction());
  #endif
  
  for (int i=0; i<numberOfIterations; i++) {
    switch ( _repositoryState.getAction()) {
      case ExaHyPE::records::RepositoryState::UseAdapterCreateGrid: watch.startTimer(); _gridWithCreateGrid.iterate(); watch.stopTimer(); _measureCreateGridCPUTime.setValue( watch.getCPUTime() ); _measureCreateGridCalendarTime.setValue( watch.getCalendarTime() ); break;

      case ExaHyPE::records::RepositoryState::Terminate:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case ExaHyPE::records::RepositoryState::NumberOfAdapters:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case ExaHyPE::records::RepositoryState::RunOnAllNodes:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case ExaHyPE::records::RepositoryState::ReadCheckpoint:
        assertionMsg( false, "not implemented yet" );
        break;
      case ExaHyPE::records::RepositoryState::WriteCheckpoint:
        assertionMsg( false, "not implemented yet" );
        break;
    }
  }
  
  #ifdef Parallel
  if (_solverState.isJoiningWithMaster()) {
    _repositoryState.setAction( ExaHyPE::records::RepositoryState::Terminate );
  }
  #endif
}

 void ExaHyPE::repositories::RepositorySTDStack::switchToCreateGrid() { _repositoryState.setAction(ExaHyPE::records::RepositoryState::UseAdapterCreateGrid); }



 bool ExaHyPE::repositories::RepositorySTDStack::isActiveAdapterCreateGrid() const { return _repositoryState.getAction() == ExaHyPE::records::RepositoryState::UseAdapterCreateGrid; }



peano::grid::Checkpoint<ExaHyPE::Vertex, ExaHyPE::Cell>* ExaHyPE::repositories::RepositorySTDStack::createEmptyCheckpoint() {
  return new peano::grid::Checkpoint<ExaHyPE::Vertex, ExaHyPE::Cell>();
} 


void ExaHyPE::repositories::RepositorySTDStack::writeCheckpoint(peano::grid::Checkpoint<ExaHyPE::Vertex, ExaHyPE::Cell> * const checkpoint) {
  _solverState.writeToCheckpoint( *checkpoint );
  _vertexStack.writeToCheckpoint( *checkpoint );
  _cellStack.writeToCheckpoint( *checkpoint );
} 


void ExaHyPE::repositories::RepositorySTDStack::readCheckpoint( peano::grid::Checkpoint<ExaHyPE::Vertex, ExaHyPE::Cell> const * const checkpoint ) {
  assertionMsg( checkpoint->isValid(), "checkpoint has to be valid if you call this operation" );

  _solverState.readFromCheckpoint( *checkpoint );
  _vertexStack.readFromCheckpoint( *checkpoint );
  _cellStack.readFromCheckpoint( *checkpoint );
}


void ExaHyPE::repositories::RepositorySTDStack::setMaximumMemoryFootprintForTemporaryRegularGrids(double value) {
  _regularGridContainer.setMaximumMemoryFootprintForTemporaryRegularGrids(value);
}


#ifdef Parallel
void ExaHyPE::repositories::RepositorySTDStack::runGlobalStep() {
  assertion(tarch::parallel::Node::getInstance().isGlobalMaster());

  ExaHyPE::records::RepositoryState intermediateStateForWorkingNodes;
  intermediateStateForWorkingNodes.setAction( ExaHyPE::records::RepositoryState::RunOnAllNodes );
  
  tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
    intermediateStateForWorkingNodes,
    peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
  );
  tarch::parallel::NodePool::getInstance().activateIdleNodes();
}


ExaHyPE::repositories::RepositorySTDStack::ContinueCommand ExaHyPE::repositories::RepositorySTDStack::continueToIterate() {
  logTraceIn( "continueToIterate()" );

  assertion( !tarch::parallel::Node::getInstance().isGlobalMaster());

  ContinueCommand result;
  if ( _solverState.hasJoinedWithMaster() ) {
    result = Terminate;
  }
  else {
    int masterNode = tarch::parallel::Node::getInstance().getGlobalMasterRank();
    assertion( masterNode != -1 );

    _repositoryState.receive( masterNode, peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag(), true, ReceiveIterationControlMessagesBlocking );

    result = Continue;
    if (_repositoryState.getAction()==ExaHyPE::records::RepositoryState::Terminate) {
      result = Terminate;
    } 
    if (_repositoryState.getAction()==ExaHyPE::records::RepositoryState::RunOnAllNodes) {
      result = RunGlobalStep;
    } 
  }
   
  logTraceOutWith1Argument( "continueToIterate()", result );
  return result;
}
#endif


void ExaHyPE::repositories::RepositorySTDStack::logIterationStatistics() const {
  logInfo( "logIterationStatistics()", "|| adapter name \t || iterations \t || total CPU time [t]=s \t || average CPU time [t]=s \t || total user time [t]=s \t || average user time [t]=s  || CPU time properties  || user time properties " );  
   logInfo( "logIterationStatistics()", "| CreateGrid \t |  " << _measureCreateGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCreateGridCPUTime.getAccumulatedValue() << " \t |  " << _measureCreateGridCPUTime.getValue()  << " \t |  " << _measureCreateGridCalendarTime.getAccumulatedValue() << " \t |  " << _measureCreateGridCalendarTime.getValue() << " \t |  " << _measureCreateGridCPUTime.toString() << " \t |  " << _measureCreateGridCalendarTime.toString() );

}


void ExaHyPE::repositories::RepositorySTDStack::clearIterationStatistics() {
   _measureCreateGridCPUTime.erase();

   _measureCreateGridCalendarTime.erase();

}
