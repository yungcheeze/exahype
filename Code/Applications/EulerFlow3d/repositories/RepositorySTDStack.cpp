#include "EulerFlow3d/repositories/RepositorySTDStack.h"

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


tarch::logging::Log exahype::repositories::RepositorySTDStack::_log( "exahype::repositories::RepositorySTDStack" );


exahype::repositories::RepositorySTDStack::RepositorySTDStack(
  peano::geometry::Geometry&                   geometry,
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset
):
  _geometry(geometry),
  _cellStack(),
  _vertexStack(),
  _solverState(),
  _gridWithInitialGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGridExport(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPatchInit(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialCondition(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictor(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrector(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionExport(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositorySTDStack(...)" );
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositorySTDStack(...)" );
}



exahype::repositories::RepositorySTDStack::RepositorySTDStack(
  peano::geometry::Geometry&                   geometry
):
  _geometry(geometry),
  _cellStack(),
  _vertexStack(),
  _solverState(),
  _gridWithInitialGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGridExport(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPatchInit(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialCondition(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictor(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrector(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionExport(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositorySTDStack(Geometry&)" );

  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );
  
  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositorySTDStack(Geometry&)" );
}
    
   
exahype::repositories::RepositorySTDStack::~RepositorySTDStack() {
  assertionMsg( _repositoryState.getAction() == exahype::records::RepositoryState::Terminate, "terminate() must be called before destroying repository." );
}


void exahype::repositories::RepositorySTDStack::restart(
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
  
  assertion( _repositoryState.getAction() == exahype::records::RepositoryState::Terminate );

  _vertexStack.clear();
  _cellStack.clear();

  _gridWithInitialGrid.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithGridExport.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPatchInit.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithInitialCondition.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictor.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrector.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithSolutionExport.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);


  _solverState.restart();

  logTraceOut( "restart(...)" );
}


void exahype::repositories::RepositorySTDStack::terminate() {
  logTraceIn( "terminate()" );
  
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );
  
  #ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
      _repositoryState,
      peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
    );
  }
  peano::parallel::SendReceiveBufferPool::getInstance().terminate();
  #endif

  _gridWithInitialGrid.terminate();
  _gridWithGridExport.terminate();
  _gridWithPatchInit.terminate();
  _gridWithInitialCondition.terminate();
  _gridWithPredictor.terminate();
  _gridWithCorrector.terminate();
  _gridWithSolutionExport.terminate();


  logTraceOut( "terminate()" );
}


exahype::State& exahype::repositories::RepositorySTDStack::getState() {
  return _solverState;
}


const exahype::State& exahype::repositories::RepositorySTDStack::getState() const {
  return _solverState;
}


void exahype::repositories::RepositorySTDStack::iterate(int numberOfIterations, bool exchangeBoundaryVertices) {
  tarch::timing::Watch watch( "exahype::repositories::RepositorySTDStack", "iterate(bool)", false);
  
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
      case exahype::records::RepositoryState::UseAdapterInitialGrid: watch.startTimer(); _gridWithInitialGrid.iterate(); watch.stopTimer(); _measureInitialGridCPUTime.setValue( watch.getCPUTime() ); _measureInitialGridCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterGridExport: watch.startTimer(); _gridWithGridExport.iterate(); watch.stopTimer(); _measureGridExportCPUTime.setValue( watch.getCPUTime() ); _measureGridExportCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPatchInit: watch.startTimer(); _gridWithPatchInit.iterate(); watch.stopTimer(); _measurePatchInitCPUTime.setValue( watch.getCPUTime() ); _measurePatchInitCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterInitialCondition: watch.startTimer(); _gridWithInitialCondition.iterate(); watch.stopTimer(); _measureInitialConditionCPUTime.setValue( watch.getCPUTime() ); _measureInitialConditionCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictor: watch.startTimer(); _gridWithPredictor.iterate(); watch.stopTimer(); _measurePredictorCPUTime.setValue( watch.getCPUTime() ); _measurePredictorCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrector: watch.startTimer(); _gridWithCorrector.iterate(); watch.stopTimer(); _measureCorrectorCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterSolutionExport: watch.startTimer(); _gridWithSolutionExport.iterate(); watch.stopTimer(); _measureSolutionExportCPUTime.setValue( watch.getCPUTime() ); _measureSolutionExportCalendarTime.setValue( watch.getCalendarTime() ); break;

      case exahype::records::RepositoryState::Terminate:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case exahype::records::RepositoryState::NumberOfAdapters:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case exahype::records::RepositoryState::RunOnAllNodes:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case exahype::records::RepositoryState::ReadCheckpoint:
        assertionMsg( false, "not implemented yet" );
        break;
      case exahype::records::RepositoryState::WriteCheckpoint:
        assertionMsg( false, "not implemented yet" );
        break;
    }
  }
  
  #ifdef Parallel
  if (_solverState.isJoiningWithMaster()) {
    _repositoryState.setAction( exahype::records::RepositoryState::Terminate );
  }
  #endif
}

 void exahype::repositories::RepositorySTDStack::switchToInitialGrid() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterInitialGrid); }
 void exahype::repositories::RepositorySTDStack::switchToGridExport() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterGridExport); }
 void exahype::repositories::RepositorySTDStack::switchToPatchInit() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPatchInit); }
 void exahype::repositories::RepositorySTDStack::switchToInitialCondition() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterInitialCondition); }
 void exahype::repositories::RepositorySTDStack::switchToPredictor() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictor); }
 void exahype::repositories::RepositorySTDStack::switchToCorrector() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrector); }
 void exahype::repositories::RepositorySTDStack::switchToSolutionExport() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterSolutionExport); }



 bool exahype::repositories::RepositorySTDStack::isActiveAdapterInitialGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterInitialGrid; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterGridExport() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterGridExport; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPatchInit() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPatchInit; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterInitialCondition() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterInitialCondition; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictor() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictor; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrector() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrector; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterSolutionExport() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterSolutionExport; }



peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>* exahype::repositories::RepositorySTDStack::createEmptyCheckpoint() {
  return new peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>();
} 


void exahype::repositories::RepositorySTDStack::writeCheckpoint(peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> * const checkpoint) {
  _solverState.writeToCheckpoint( *checkpoint );
  _vertexStack.writeToCheckpoint( *checkpoint );
  _cellStack.writeToCheckpoint( *checkpoint );
} 


void exahype::repositories::RepositorySTDStack::readCheckpoint( peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> const * const checkpoint ) {
  assertionMsg( checkpoint->isValid(), "checkpoint has to be valid if you call this operation" );

  _solverState.readFromCheckpoint( *checkpoint );
  _vertexStack.readFromCheckpoint( *checkpoint );
  _cellStack.readFromCheckpoint( *checkpoint );
}


void exahype::repositories::RepositorySTDStack::setMaximumMemoryFootprintForTemporaryRegularGrids(double value) {
  _regularGridContainer.setMaximumMemoryFootprintForTemporaryRegularGrids(value);
}


#ifdef Parallel
void exahype::repositories::RepositorySTDStack::runGlobalStep() {
  assertion(tarch::parallel::Node::getInstance().isGlobalMaster());

  exahype::records::RepositoryState intermediateStateForWorkingNodes;
  intermediateStateForWorkingNodes.setAction( exahype::records::RepositoryState::RunOnAllNodes );
  
  tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
    intermediateStateForWorkingNodes,
    peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
  );
  tarch::parallel::NodePool::getInstance().activateIdleNodes();
}


exahype::repositories::RepositorySTDStack::ContinueCommand exahype::repositories::RepositorySTDStack::continueToIterate() {
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
    if (_repositoryState.getAction()==exahype::records::RepositoryState::Terminate) {
      result = Terminate;
    } 
    if (_repositoryState.getAction()==exahype::records::RepositoryState::RunOnAllNodes) {
      result = RunGlobalStep;
    } 
  }
   
  logTraceOutWith1Argument( "continueToIterate()", result );
  return result;
}
#endif


void exahype::repositories::RepositorySTDStack::logIterationStatistics() const {
  logInfo( "logIterationStatistics()", "|| adapter name \t || iterations \t || total CPU time [t]=s \t || average CPU time [t]=s \t || total user time [t]=s \t || average user time [t]=s  || CPU time properties  || user time properties " );  
   logInfo( "logIterationStatistics()", "| InitialGrid \t |  " << _measureInitialGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measureInitialGridCPUTime.getAccumulatedValue() << " \t |  " << _measureInitialGridCPUTime.getValue()  << " \t |  " << _measureInitialGridCalendarTime.getAccumulatedValue() << " \t |  " << _measureInitialGridCalendarTime.getValue() << " \t |  " << _measureInitialGridCPUTime.toString() << " \t |  " << _measureInitialGridCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| GridExport \t |  " << _measureGridExportCPUTime.getNumberOfMeasurements() << " \t |  " << _measureGridExportCPUTime.getAccumulatedValue() << " \t |  " << _measureGridExportCPUTime.getValue()  << " \t |  " << _measureGridExportCalendarTime.getAccumulatedValue() << " \t |  " << _measureGridExportCalendarTime.getValue() << " \t |  " << _measureGridExportCPUTime.toString() << " \t |  " << _measureGridExportCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| PatchInit \t |  " << _measurePatchInitCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePatchInitCPUTime.getAccumulatedValue() << " \t |  " << _measurePatchInitCPUTime.getValue()  << " \t |  " << _measurePatchInitCalendarTime.getAccumulatedValue() << " \t |  " << _measurePatchInitCalendarTime.getValue() << " \t |  " << _measurePatchInitCPUTime.toString() << " \t |  " << _measurePatchInitCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| InitialCondition \t |  " << _measureInitialConditionCPUTime.getNumberOfMeasurements() << " \t |  " << _measureInitialConditionCPUTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionCPUTime.getValue()  << " \t |  " << _measureInitialConditionCalendarTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionCalendarTime.getValue() << " \t |  " << _measureInitialConditionCPUTime.toString() << " \t |  " << _measureInitialConditionCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| Predictor \t |  " << _measurePredictorCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorCPUTime.getValue()  << " \t |  " << _measurePredictorCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorCalendarTime.getValue() << " \t |  " << _measurePredictorCPUTime.toString() << " \t |  " << _measurePredictorCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| Corrector \t |  " << _measureCorrectorCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorCPUTime.getValue()  << " \t |  " << _measureCorrectorCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorCalendarTime.getValue() << " \t |  " << _measureCorrectorCPUTime.toString() << " \t |  " << _measureCorrectorCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| SolutionExport \t |  " << _measureSolutionExportCPUTime.getNumberOfMeasurements() << " \t |  " << _measureSolutionExportCPUTime.getAccumulatedValue() << " \t |  " << _measureSolutionExportCPUTime.getValue()  << " \t |  " << _measureSolutionExportCalendarTime.getAccumulatedValue() << " \t |  " << _measureSolutionExportCalendarTime.getValue() << " \t |  " << _measureSolutionExportCPUTime.toString() << " \t |  " << _measureSolutionExportCalendarTime.toString() );

}


void exahype::repositories::RepositorySTDStack::clearIterationStatistics() {
   _measureInitialGridCPUTime.erase();
   _measureGridExportCPUTime.erase();
   _measurePatchInitCPUTime.erase();
   _measureInitialConditionCPUTime.erase();
   _measurePredictorCPUTime.erase();
   _measureCorrectorCPUTime.erase();
   _measureSolutionExportCPUTime.erase();

   _measureInitialGridCalendarTime.erase();
   _measureGridExportCalendarTime.erase();
   _measurePatchInitCalendarTime.erase();
   _measureInitialConditionCalendarTime.erase();
   _measurePredictorCalendarTime.erase();
   _measureCorrectorCalendarTime.erase();
   _measureSolutionExportCalendarTime.erase();

}
