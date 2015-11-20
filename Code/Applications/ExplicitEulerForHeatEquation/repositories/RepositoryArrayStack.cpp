#include "ExplicitEulerForHeatEquation/repositories/RepositoryArrayStack.h"

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


tarch::logging::Log myproject::repositories::RepositoryArrayStack::_log( "myproject::repositories::RepositoryArrayStack" );


myproject::repositories::RepositoryArrayStack::RepositoryArrayStack(
  peano::geometry::Geometry&                   geometry,
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
  int                                          maximumSizeOfCellInOutStack,
  int                                          maximumSizeOfVertexInOutStack,
  int                                          maximumSizeOfVertexTemporaryStack
):
  _geometry(geometry),
  _cellStack(maximumSizeOfCellInOutStack),
  _vertexStack(maximumSizeOfVertexInOutStack, maximumSizeOfVertexTemporaryStack),
  _solverState(),
  _gridWithCreateGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithTimeStep(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCreateGridAndPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithTimeStepAndPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositoryArrayStack(...)" );
  
  _repositoryState.setAction( myproject::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(myproject::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(myproject::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositoryArrayStack(...)" );
}



myproject::repositories::RepositoryArrayStack::RepositoryArrayStack(
  peano::geometry::Geometry&                   geometry,
  int                                          maximumSizeOfCellInOutStack,
  int                                          maximumSizeOfVertexInOutStack,
  int                                          maximumSizeOfVertexTemporaryStack
):
  _geometry(geometry),
  _cellStack(maximumSizeOfCellInOutStack),
  _vertexStack(maximumSizeOfVertexInOutStack,maximumSizeOfVertexTemporaryStack),
  _solverState(),
  _gridWithCreateGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithTimeStep(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCreateGridAndPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithTimeStepAndPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositoryArrayStack(Geometry&)" );
  
  _repositoryState.setAction( myproject::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(myproject::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(myproject::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositoryArrayStack(Geometry&)" );
}
    
   
myproject::repositories::RepositoryArrayStack::~RepositoryArrayStack() {
  assertion( _repositoryState.getAction() == myproject::records::RepositoryState::Terminate );
}


void myproject::repositories::RepositoryArrayStack::restart(
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
  
  assertion( _repositoryState.getAction() == myproject::records::RepositoryState::Terminate );

  _vertexStack.clear();
  _cellStack.clear();

  _gridWithCreateGrid.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithTimeStep.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCreateGridAndPlot.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithTimeStepAndPlot.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);

 
   _solverState.restart();
 
  logTraceOut( "restart(...)" );
}


void myproject::repositories::RepositoryArrayStack::terminate() {
  logTraceIn( "terminate()" );
  
  _repositoryState.setAction( myproject::records::RepositoryState::Terminate );
  
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
  _gridWithTimeStep.terminate();
  _gridWithCreateGridAndPlot.terminate();
  _gridWithTimeStepAndPlot.terminate();

 
  logTraceOut( "terminate()" );
}


myproject::State& myproject::repositories::RepositoryArrayStack::getState() {
  return _solverState;
}


const myproject::State& myproject::repositories::RepositoryArrayStack::getState() const {
  return _solverState;
}

   
void myproject::repositories::RepositoryArrayStack::iterate(int numberOfIterations, bool exchangeBoundaryVertices) {
  tarch::timing::Watch watch( "myproject::repositories::RepositoryArrayStack", "iterate(bool)", false);
  
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
      case myproject::records::RepositoryState::UseAdapterCreateGrid: watch.startTimer(); _gridWithCreateGrid.iterate(); watch.stopTimer(); _measureCreateGridCPUTime.setValue( watch.getCPUTime() ); _measureCreateGridCalendarTime.setValue( watch.getCalendarTime() ); break;
      case myproject::records::RepositoryState::UseAdapterTimeStep: watch.startTimer(); _gridWithTimeStep.iterate(); watch.stopTimer(); _measureTimeStepCPUTime.setValue( watch.getCPUTime() ); _measureTimeStepCalendarTime.setValue( watch.getCalendarTime() ); break;
      case myproject::records::RepositoryState::UseAdapterCreateGridAndPlot: watch.startTimer(); _gridWithCreateGridAndPlot.iterate(); watch.stopTimer(); _measureCreateGridAndPlotCPUTime.setValue( watch.getCPUTime() ); _measureCreateGridAndPlotCalendarTime.setValue( watch.getCalendarTime() ); break;
      case myproject::records::RepositoryState::UseAdapterTimeStepAndPlot: watch.startTimer(); _gridWithTimeStepAndPlot.iterate(); watch.stopTimer(); _measureTimeStepAndPlotCPUTime.setValue( watch.getCPUTime() ); _measureTimeStepAndPlotCalendarTime.setValue( watch.getCalendarTime() ); break;

      case myproject::records::RepositoryState::Terminate:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case myproject::records::RepositoryState::NumberOfAdapters:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case myproject::records::RepositoryState::RunOnAllNodes:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case myproject::records::RepositoryState::ReadCheckpoint:
        assertionMsg( false, "not implemented yet" );
        break;
      case myproject::records::RepositoryState::WriteCheckpoint:
        assertionMsg( false, "not implemented yet" );
        break;
    }
  }
    
  #ifdef Parallel
  if (_solverState.isJoiningWithMaster()) {
    _repositoryState.setAction( myproject::records::RepositoryState::Terminate );
  }
  #endif
}

 void myproject::repositories::RepositoryArrayStack::switchToCreateGrid() { _repositoryState.setAction(myproject::records::RepositoryState::UseAdapterCreateGrid); }
 void myproject::repositories::RepositoryArrayStack::switchToTimeStep() { _repositoryState.setAction(myproject::records::RepositoryState::UseAdapterTimeStep); }
 void myproject::repositories::RepositoryArrayStack::switchToCreateGridAndPlot() { _repositoryState.setAction(myproject::records::RepositoryState::UseAdapterCreateGridAndPlot); }
 void myproject::repositories::RepositoryArrayStack::switchToTimeStepAndPlot() { _repositoryState.setAction(myproject::records::RepositoryState::UseAdapterTimeStepAndPlot); }



 bool myproject::repositories::RepositoryArrayStack::isActiveAdapterCreateGrid() const { return _repositoryState.getAction() == myproject::records::RepositoryState::UseAdapterCreateGrid; }
 bool myproject::repositories::RepositoryArrayStack::isActiveAdapterTimeStep() const { return _repositoryState.getAction() == myproject::records::RepositoryState::UseAdapterTimeStep; }
 bool myproject::repositories::RepositoryArrayStack::isActiveAdapterCreateGridAndPlot() const { return _repositoryState.getAction() == myproject::records::RepositoryState::UseAdapterCreateGridAndPlot; }
 bool myproject::repositories::RepositoryArrayStack::isActiveAdapterTimeStepAndPlot() const { return _repositoryState.getAction() == myproject::records::RepositoryState::UseAdapterTimeStepAndPlot; }



peano::grid::Checkpoint<myproject::Vertex, myproject::Cell>* myproject::repositories::RepositoryArrayStack::createEmptyCheckpoint() {
  return new peano::grid::Checkpoint<myproject::Vertex, myproject::Cell>();
} 


void myproject::repositories::RepositoryArrayStack::writeCheckpoint(peano::grid::Checkpoint<myproject::Vertex, myproject::Cell> * const checkpoint) {
  _solverState.writeToCheckpoint( *checkpoint );
  _vertexStack.writeToCheckpoint( *checkpoint );
  _cellStack.writeToCheckpoint( *checkpoint );
} 


void myproject::repositories::RepositoryArrayStack::setMaximumMemoryFootprintForTemporaryRegularGrids(double value) {
  _regularGridContainer.setMaximumMemoryFootprintForTemporaryRegularGrids(value);
}


void myproject::repositories::RepositoryArrayStack::readCheckpoint( peano::grid::Checkpoint<myproject::Vertex, myproject::Cell> const * const checkpoint ) {
  assertionMsg( checkpoint->isValid(), "checkpoint has to be valid if you call this operation" );

  _solverState.readFromCheckpoint( *checkpoint );
  _vertexStack.readFromCheckpoint( *checkpoint );
  _cellStack.readFromCheckpoint( *checkpoint );
}


#ifdef Parallel
void myproject::repositories::RepositoryArrayStack::runGlobalStep() {
  assertion(tarch::parallel::Node::getInstance().isGlobalMaster());

  myproject::records::RepositoryState intermediateStateForWorkingNodes;
  intermediateStateForWorkingNodes.setAction( myproject::records::RepositoryState::RunOnAllNodes );
  
  tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
    intermediateStateForWorkingNodes,
    peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
  );
  tarch::parallel::NodePool::getInstance().activateIdleNodes();
}


myproject::repositories::RepositoryArrayStack::ContinueCommand myproject::repositories::RepositoryArrayStack::continueToIterate() {
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
    if (_repositoryState.getAction()==myproject::records::RepositoryState::Terminate) {
      result = Terminate;
    } 
    if (_repositoryState.getAction()==myproject::records::RepositoryState::RunOnAllNodes) {
      result = RunGlobalStep;
    } 
  }
   
  logTraceOutWith1Argument( "continueToIterate()", result );
  return result;
}
#endif


void myproject::repositories::RepositoryArrayStack::logIterationStatistics() const {
  logInfo( "logIterationStatistics()", "|| adapter name \t || iterations \t || total CPU time [t]=s \t || average CPU time [t]=s \t || total user time [t]=s \t || average user time [t]=s  || CPU time properties  || user time properties " );  
   logInfo( "logIterationStatistics()", "| CreateGrid \t |  " << _measureCreateGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCreateGridCPUTime.getAccumulatedValue() << " \t |  " << _measureCreateGridCPUTime.getValue()  << " \t |  " << _measureCreateGridCalendarTime.getAccumulatedValue() << " \t |  " << _measureCreateGridCalendarTime.getValue() << " \t |  " << _measureCreateGridCPUTime.toString() << " \t |  " << _measureCreateGridCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| TimeStep \t |  " << _measureTimeStepCPUTime.getNumberOfMeasurements() << " \t |  " << _measureTimeStepCPUTime.getAccumulatedValue() << " \t |  " << _measureTimeStepCPUTime.getValue()  << " \t |  " << _measureTimeStepCalendarTime.getAccumulatedValue() << " \t |  " << _measureTimeStepCalendarTime.getValue() << " \t |  " << _measureTimeStepCPUTime.toString() << " \t |  " << _measureTimeStepCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| CreateGridAndPlot \t |  " << _measureCreateGridAndPlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCreateGridAndPlotCPUTime.getAccumulatedValue() << " \t |  " << _measureCreateGridAndPlotCPUTime.getValue()  << " \t |  " << _measureCreateGridAndPlotCalendarTime.getAccumulatedValue() << " \t |  " << _measureCreateGridAndPlotCalendarTime.getValue() << " \t |  " << _measureCreateGridAndPlotCPUTime.toString() << " \t |  " << _measureCreateGridAndPlotCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| TimeStepAndPlot \t |  " << _measureTimeStepAndPlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measureTimeStepAndPlotCPUTime.getAccumulatedValue() << " \t |  " << _measureTimeStepAndPlotCPUTime.getValue()  << " \t |  " << _measureTimeStepAndPlotCalendarTime.getAccumulatedValue() << " \t |  " << _measureTimeStepAndPlotCalendarTime.getValue() << " \t |  " << _measureTimeStepAndPlotCPUTime.toString() << " \t |  " << _measureTimeStepAndPlotCalendarTime.toString() );

}


void myproject::repositories::RepositoryArrayStack::clearIterationStatistics() {
   _measureCreateGridCPUTime.erase();
   _measureTimeStepCPUTime.erase();
   _measureCreateGridAndPlotCPUTime.erase();
   _measureTimeStepAndPlotCPUTime.erase();

   _measureCreateGridCalendarTime.erase();
   _measureTimeStepCalendarTime.erase();
   _measureCreateGridAndPlotCalendarTime.erase();
   _measureTimeStepAndPlotCalendarTime.erase();

}
