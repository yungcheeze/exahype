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
  _gridWithPatchInitialisation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPatchInitialisationAndExport(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialCondition(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndExport(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictor(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithFaceDataExchange(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrector(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndExport(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPredictor(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPredictorAndExport(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputationAndExport(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
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
  _gridWithPatchInitialisation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPatchInitialisationAndExport(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialCondition(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndExport(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictor(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithFaceDataExchange(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrector(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndExport(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPredictor(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPredictorAndExport(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputationAndExport(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
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
  _gridWithPatchInitialisation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPatchInitialisationAndExport.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithInitialCondition.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithInitialConditionAndExport.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictor.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithFaceDataExchange.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrector.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrectorAndExport.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrectorAndPredictor.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrectorAndPredictorAndExport.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputationAndExport.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
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
  _gridWithPatchInitialisation.terminate();
  _gridWithPatchInitialisationAndExport.terminate();
  _gridWithInitialCondition.terminate();
  _gridWithInitialConditionAndExport.terminate();
  _gridWithGlobalTimeStepComputation.terminate();
  _gridWithPredictor.terminate();
  _gridWithFaceDataExchange.terminate();
  _gridWithCorrector.terminate();
  _gridWithCorrectorAndExport.terminate();
  _gridWithCorrectorAndPredictor.terminate();
  _gridWithCorrectorAndPredictorAndExport.terminate();
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputation.terminate();
  _gridWithCorrectorAndPredictorAndGlobalTimeStepComputationAndExport.terminate();
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
      case exahype::records::RepositoryState::UseAdapterPatchInitialisation: watch.startTimer(); _gridWithPatchInitialisation.iterate(); watch.stopTimer(); _measurePatchInitialisationCPUTime.setValue( watch.getCPUTime() ); _measurePatchInitialisationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPatchInitialisationAndExport: watch.startTimer(); _gridWithPatchInitialisationAndExport.iterate(); watch.stopTimer(); _measurePatchInitialisationAndExportCPUTime.setValue( watch.getCPUTime() ); _measurePatchInitialisationAndExportCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterInitialCondition: watch.startTimer(); _gridWithInitialCondition.iterate(); watch.stopTimer(); _measureInitialConditionCPUTime.setValue( watch.getCPUTime() ); _measureInitialConditionCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterInitialConditionAndExport: watch.startTimer(); _gridWithInitialConditionAndExport.iterate(); watch.stopTimer(); _measureInitialConditionAndExportCPUTime.setValue( watch.getCPUTime() ); _measureInitialConditionAndExportCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterGlobalTimeStepComputation: watch.startTimer(); _gridWithGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measureGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measureGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictor: watch.startTimer(); _gridWithPredictor.iterate(); watch.stopTimer(); _measurePredictorCPUTime.setValue( watch.getCPUTime() ); _measurePredictorCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterFaceDataExchange: watch.startTimer(); _gridWithFaceDataExchange.iterate(); watch.stopTimer(); _measureFaceDataExchangeCPUTime.setValue( watch.getCPUTime() ); _measureFaceDataExchangeCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrector: watch.startTimer(); _gridWithCorrector.iterate(); watch.stopTimer(); _measureCorrectorCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrectorAndExport: watch.startTimer(); _gridWithCorrectorAndExport.iterate(); watch.stopTimer(); _measureCorrectorAndExportCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorAndExportCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrectorAndPredictor: watch.startTimer(); _gridWithCorrectorAndPredictor.iterate(); watch.stopTimer(); _measureCorrectorAndPredictorCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorAndPredictorCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndExport: watch.startTimer(); _gridWithCorrectorAndPredictorAndExport.iterate(); watch.stopTimer(); _measureCorrectorAndPredictorAndExportCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorAndPredictorAndExportCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndGlobalTimeStepComputation: watch.startTimer(); _gridWithCorrectorAndPredictorAndGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorAndPredictorAndGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndGlobalTimeStepComputationAndExport: watch.startTimer(); _gridWithCorrectorAndPredictorAndGlobalTimeStepComputationAndExport.iterate(); watch.stopTimer(); _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCalendarTime.setValue( watch.getCalendarTime() ); break;
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
 void exahype::repositories::RepositorySTDStack::switchToPatchInitialisation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPatchInitialisation); }
 void exahype::repositories::RepositorySTDStack::switchToPatchInitialisationAndExport() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPatchInitialisationAndExport); }
 void exahype::repositories::RepositorySTDStack::switchToInitialCondition() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterInitialCondition); }
 void exahype::repositories::RepositorySTDStack::switchToInitialConditionAndExport() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterInitialConditionAndExport); }
 void exahype::repositories::RepositorySTDStack::switchToGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToPredictor() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictor); }
 void exahype::repositories::RepositorySTDStack::switchToFaceDataExchange() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterFaceDataExchange); }
 void exahype::repositories::RepositorySTDStack::switchToCorrector() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrector); }
 void exahype::repositories::RepositorySTDStack::switchToCorrectorAndExport() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrectorAndExport); }
 void exahype::repositories::RepositorySTDStack::switchToCorrectorAndPredictor() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrectorAndPredictor); }
 void exahype::repositories::RepositorySTDStack::switchToCorrectorAndPredictorAndExport() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndExport); }
 void exahype::repositories::RepositorySTDStack::switchToCorrectorAndPredictorAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToCorrectorAndPredictorAndGlobalTimeStepComputationAndExport() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndGlobalTimeStepComputationAndExport); }
 void exahype::repositories::RepositorySTDStack::switchToSolutionExport() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterSolutionExport); }



 bool exahype::repositories::RepositorySTDStack::isActiveAdapterInitialGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterInitialGrid; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterGridExport() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterGridExport; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPatchInitialisation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPatchInitialisation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPatchInitialisationAndExport() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPatchInitialisationAndExport; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterInitialCondition() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterInitialCondition; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterInitialConditionAndExport() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterInitialConditionAndExport; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictor() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictor; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterFaceDataExchange() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterFaceDataExchange; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrector() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrector; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrectorAndExport() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrectorAndExport; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrectorAndPredictor() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrectorAndPredictor; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrectorAndPredictorAndExport() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndExport; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrectorAndPredictorAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrectorAndPredictorAndGlobalTimeStepComputationAndExport() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrectorAndPredictorAndGlobalTimeStepComputationAndExport; }
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
   logInfo( "logIterationStatistics()", "| PatchInitialisation \t |  " << _measurePatchInitialisationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePatchInitialisationCPUTime.getAccumulatedValue() << " \t |  " << _measurePatchInitialisationCPUTime.getValue()  << " \t |  " << _measurePatchInitialisationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePatchInitialisationCalendarTime.getValue() << " \t |  " << _measurePatchInitialisationCPUTime.toString() << " \t |  " << _measurePatchInitialisationCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| PatchInitialisationAndExport \t |  " << _measurePatchInitialisationAndExportCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePatchInitialisationAndExportCPUTime.getAccumulatedValue() << " \t |  " << _measurePatchInitialisationAndExportCPUTime.getValue()  << " \t |  " << _measurePatchInitialisationAndExportCalendarTime.getAccumulatedValue() << " \t |  " << _measurePatchInitialisationAndExportCalendarTime.getValue() << " \t |  " << _measurePatchInitialisationAndExportCPUTime.toString() << " \t |  " << _measurePatchInitialisationAndExportCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| InitialCondition \t |  " << _measureInitialConditionCPUTime.getNumberOfMeasurements() << " \t |  " << _measureInitialConditionCPUTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionCPUTime.getValue()  << " \t |  " << _measureInitialConditionCalendarTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionCalendarTime.getValue() << " \t |  " << _measureInitialConditionCPUTime.toString() << " \t |  " << _measureInitialConditionCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| InitialConditionAndExport \t |  " << _measureInitialConditionAndExportCPUTime.getNumberOfMeasurements() << " \t |  " << _measureInitialConditionAndExportCPUTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndExportCPUTime.getValue()  << " \t |  " << _measureInitialConditionAndExportCalendarTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndExportCalendarTime.getValue() << " \t |  " << _measureInitialConditionAndExportCPUTime.toString() << " \t |  " << _measureInitialConditionAndExportCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| GlobalTimeStepComputation \t |  " << _measureGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measureGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measureGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measureGlobalTimeStepComputationCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| Predictor \t |  " << _measurePredictorCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorCPUTime.getValue()  << " \t |  " << _measurePredictorCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorCalendarTime.getValue() << " \t |  " << _measurePredictorCPUTime.toString() << " \t |  " << _measurePredictorCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| FaceDataExchange \t |  " << _measureFaceDataExchangeCPUTime.getNumberOfMeasurements() << " \t |  " << _measureFaceDataExchangeCPUTime.getAccumulatedValue() << " \t |  " << _measureFaceDataExchangeCPUTime.getValue()  << " \t |  " << _measureFaceDataExchangeCalendarTime.getAccumulatedValue() << " \t |  " << _measureFaceDataExchangeCalendarTime.getValue() << " \t |  " << _measureFaceDataExchangeCPUTime.toString() << " \t |  " << _measureFaceDataExchangeCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| Corrector \t |  " << _measureCorrectorCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorCPUTime.getValue()  << " \t |  " << _measureCorrectorCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorCalendarTime.getValue() << " \t |  " << _measureCorrectorCPUTime.toString() << " \t |  " << _measureCorrectorCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| CorrectorAndExport \t |  " << _measureCorrectorAndExportCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorAndExportCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndExportCPUTime.getValue()  << " \t |  " << _measureCorrectorAndExportCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndExportCalendarTime.getValue() << " \t |  " << _measureCorrectorAndExportCPUTime.toString() << " \t |  " << _measureCorrectorAndExportCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| CorrectorAndPredictor \t |  " << _measureCorrectorAndPredictorCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorAndPredictorCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPredictorCPUTime.getValue()  << " \t |  " << _measureCorrectorAndPredictorCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPredictorCalendarTime.getValue() << " \t |  " << _measureCorrectorAndPredictorCPUTime.toString() << " \t |  " << _measureCorrectorAndPredictorCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| CorrectorAndPredictorAndExport \t |  " << _measureCorrectorAndPredictorAndExportCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorAndPredictorAndExportCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPredictorAndExportCPUTime.getValue()  << " \t |  " << _measureCorrectorAndPredictorAndExportCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPredictorAndExportCalendarTime.getValue() << " \t |  " << _measureCorrectorAndPredictorAndExportCPUTime.toString() << " \t |  " << _measureCorrectorAndPredictorAndExportCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| CorrectorAndPredictorAndGlobalTimeStepComputation \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| CorrectorAndPredictorAndGlobalTimeStepComputationAndExport \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime.getValue()  << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCalendarTime.getValue() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime.toString() << " \t |  " << _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCalendarTime.toString() );
   logInfo( "logIterationStatistics()", "| SolutionExport \t |  " << _measureSolutionExportCPUTime.getNumberOfMeasurements() << " \t |  " << _measureSolutionExportCPUTime.getAccumulatedValue() << " \t |  " << _measureSolutionExportCPUTime.getValue()  << " \t |  " << _measureSolutionExportCalendarTime.getAccumulatedValue() << " \t |  " << _measureSolutionExportCalendarTime.getValue() << " \t |  " << _measureSolutionExportCPUTime.toString() << " \t |  " << _measureSolutionExportCalendarTime.toString() );

}


void exahype::repositories::RepositorySTDStack::clearIterationStatistics() {
   _measureInitialGridCPUTime.erase();
   _measureGridExportCPUTime.erase();
   _measurePatchInitialisationCPUTime.erase();
   _measurePatchInitialisationAndExportCPUTime.erase();
   _measureInitialConditionCPUTime.erase();
   _measureInitialConditionAndExportCPUTime.erase();
   _measureGlobalTimeStepComputationCPUTime.erase();
   _measurePredictorCPUTime.erase();
   _measureFaceDataExchangeCPUTime.erase();
   _measureCorrectorCPUTime.erase();
   _measureCorrectorAndExportCPUTime.erase();
   _measureCorrectorAndPredictorCPUTime.erase();
   _measureCorrectorAndPredictorAndExportCPUTime.erase();
   _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime.erase();
   _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime.erase();
   _measureSolutionExportCPUTime.erase();

   _measureInitialGridCalendarTime.erase();
   _measureGridExportCalendarTime.erase();
   _measurePatchInitialisationCalendarTime.erase();
   _measurePatchInitialisationAndExportCalendarTime.erase();
   _measureInitialConditionCalendarTime.erase();
   _measureInitialConditionAndExportCalendarTime.erase();
   _measureGlobalTimeStepComputationCalendarTime.erase();
   _measurePredictorCalendarTime.erase();
   _measureFaceDataExchangeCalendarTime.erase();
   _measureCorrectorCalendarTime.erase();
   _measureCorrectorAndExportCalendarTime.erase();
   _measureCorrectorAndPredictorCalendarTime.erase();
   _measureCorrectorAndPredictorAndExportCalendarTime.erase();
   _measureCorrectorAndPredictorAndGlobalTimeStepComputationCalendarTime.erase();
   _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCalendarTime.erase();
   _measureSolutionExportCalendarTime.erase();

}
