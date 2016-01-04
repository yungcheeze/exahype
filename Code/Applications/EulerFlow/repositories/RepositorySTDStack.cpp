#include "EulerFlow/repositories/RepositorySTDStack.h"

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
  _gridWithPatchInitialisation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPatchInitialisationAndExport(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithFaceDataExchange(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndExportAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputationAndExport(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),

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
  _gridWithPatchInitialisation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPatchInitialisationAndExport(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithFaceDataExchange(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndExportAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputationAndExport(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),

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
  _gridWithPatchInitialisation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPatchInitialisationAndExport.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithFaceDataExchange.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithInitialConditionAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithInitialConditionAndExportAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictorAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputationAndExport.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);


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
  _gridWithPatchInitialisation.terminate();
  _gridWithPatchInitialisationAndExport.terminate();
  _gridWithFaceDataExchange.terminate();
  _gridWithInitialConditionAndGlobalTimeStepComputation.terminate();
  _gridWithInitialConditionAndExportAndGlobalTimeStepComputation.terminate();
  _gridWithPredictorAndGlobalTimeStepComputation.terminate();
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputation.terminate();
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputationAndExport.terminate();


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
      case exahype::records::RepositoryState::UseAdapterPatchInitialisation: watch.startTimer(); _gridWithPatchInitialisation.iterate(); watch.stopTimer(); _measurePatchInitialisationCPUTime.setValue( watch.getCPUTime() ); _measurePatchInitialisationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPatchInitialisationAndExport: watch.startTimer(); _gridWithPatchInitialisationAndExport.iterate(); watch.stopTimer(); _measurePatchInitialisationAndExportCPUTime.setValue( watch.getCPUTime() ); _measurePatchInitialisationAndExportCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterFaceDataExchange: watch.startTimer(); _gridWithFaceDataExchange.iterate(); watch.stopTimer(); _measureFaceDataExchangeCPUTime.setValue( watch.getCPUTime() ); _measureFaceDataExchangeCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterInitialConditionAndGlobalTimeStepComputation: watch.startTimer(); _gridWithInitialConditionAndGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measureInitialConditionAndGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measureInitialConditionAndGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterInitialConditionAndExportAndGlobalTimeStepComputation: watch.startTimer(); _gridWithInitialConditionAndExportAndGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measureInitialConditionAndExportAndGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measureInitialConditionAndExportAndGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictorAndGlobalTimeStepComputation: watch.startTimer(); _gridWithPredictorAndGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measurePredictorAndGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measurePredictorAndGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndGlobalTimeStepComputation: watch.startTimer(); _gridWithCorrectorAndPredictorAndGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorAndPredictorAndGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndGlobalTimeStepComputationAndExport: watch.startTimer(); _gridWithCorrectorAndPredictorAndGlobalTimeStepComputationAndExport.iterate(); watch.stopTimer(); _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCalendarTime.setValue( watch.getCalendarTime() ); break;

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
 void exahype::repositories::RepositorySTDStack::switchToPatchInitialisation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPatchInitialisation); }
 void exahype::repositories::RepositorySTDStack::switchToPatchInitialisationAndExport() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPatchInitialisationAndExport); }
 void exahype::repositories::RepositorySTDStack::switchToFaceDataExchange() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterFaceDataExchange); }
 void exahype::repositories::RepositorySTDStack::switchToInitialConditionAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterInitialConditionAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToInitialConditionAndExportAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterInitialConditionAndExportAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToPredictorAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictorAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToCorrectorAndPredictorAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToCorrectorAndPredictorAndGlobalTimeStepComputationAndExport() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndGlobalTimeStepComputationAndExport); }



 bool exahype::repositories::RepositorySTDStack::isActiveAdapterInitialGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterInitialGrid; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterGridExport() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterGridExport; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPatchInitialisation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPatchInitialisation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPatchInitialisationAndExport() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPatchInitialisationAndExport; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterFaceDataExchange() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterFaceDataExchange; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterInitialConditionAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterInitialConditionAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterInitialConditionAndExportAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterInitialConditionAndExportAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictorAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictorAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrectorAndPredictorAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrectorAndPredictorAndGlobalTimeStepComputationAndExport() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndGlobalTimeStepComputationAndExport; }



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
   logInfo( "logIterationStatistics()", "| PatchInitialisation \t |  " << _measurePatchInitialisationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePatchInitialisationCPUTime.getAccumulatedValue() << " \t |  " << _measurePatchInitialisationCPUTime.getValue()  << " \t |  " << _measurePatchInitialisationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePatchInitialisationCalendarTime.getValue() << " \t |  " << _measurePatchInitialisationCPUTime.toString() << " \t |  " << _measurePatchInitialisationCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| PatchInitialisationAndExport \t |  " << _measurePatchInitialisationAndExportCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePatchInitialisationAndExportCPUTime.getAccumulatedValue() << " \t |  " << _measurePatchInitialisationAndExportCPUTime.getValue()  << " \t |  " << _measurePatchInitialisationAndExportCalendarTime.getAccumulatedValue() << " \t |  " << _measurePatchInitialisationAndExportCalendarTime.getValue() << " \t |  " << _measurePatchInitialisationAndExportCPUTime.toString() << " \t |  " << _measurePatchInitialisationAndExportCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| FaceDataExchange \t |  " << _measureFaceDataExchangeCPUTime.getNumberOfMeasurements() << " \t |  " << _measureFaceDataExchangeCPUTime.getAccumulatedValue() << " \t |  " << _measureFaceDataExchangeCPUTime.getValue()  << " \t |  " << _measureFaceDataExchangeCalendarTime.getAccumulatedValue() << " \t |  " << _measureFaceDataExchangeCalendarTime.getValue() << " \t |  " << _measureFaceDataExchangeCPUTime.toString() << " \t |  " << _measureFaceDataExchangeCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| InitialConditionAndGlobalTimeStepComputation \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| InitialConditionAndExportAndGlobalTimeStepComputation \t |  " << _measureInitialConditionAndExportAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureInitialConditionAndExportAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndExportAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measureInitialConditionAndExportAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndExportAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measureInitialConditionAndExportAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measureInitialConditionAndExportAndGlobalTimeStepComputationCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| PredictorAndGlobalTimeStepComputation \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| CorrectorAndPredictorAndGlobalTimeStepComputation \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| CorrectorAndPredictorAndGlobalTimeStepComputationAndExport \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime.getValue()  << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCalendarTime.getValue() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime.toString() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCalendarTime.toString() );

}


void exahype::repositories::RepositorySTDStack::clearIterationStatistics() {
   _measureInitialGridCPUTime.erase();
   _measureGridExportCPUTime.erase();
   _measurePatchInitialisationCPUTime.erase();
   _measurePatchInitialisationAndExportCPUTime.erase();
   _measureFaceDataExchangeCPUTime.erase();
   _measureInitialConditionAndGlobalTimeStepComputationCPUTime.erase();
   _measureInitialConditionAndExportAndGlobalTimeStepComputationCPUTime.erase();
   _measurePredictorAndGlobalTimeStepComputationCPUTime.erase();
   _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime.erase();
   _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime.erase();

   _measureInitialGridCalendarTime.erase();
   _measureGridExportCalendarTime.erase();
   _measurePatchInitialisationCalendarTime.erase();
   _measurePatchInitialisationAndExportCalendarTime.erase();
   _measureFaceDataExchangeCalendarTime.erase();
   _measureInitialConditionAndGlobalTimeStepComputationCalendarTime.erase();
   _measureInitialConditionAndExportAndGlobalTimeStepComputationCalendarTime.erase();
   _measurePredictorAndGlobalTimeStepComputationCalendarTime.erase();
   _measureCorrectorAndPredictorAndGlobalTimeStepComputationCalendarTime.erase();
   _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCalendarTime.erase();

}
