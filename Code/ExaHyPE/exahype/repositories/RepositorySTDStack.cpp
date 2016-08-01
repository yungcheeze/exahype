#include "exahype/repositories/RepositorySTDStack.h"

#include "tarch/Assertions.h"
#include "tarch/timing/Watch.h"

#include "tarch/compiler/CompilerSpecificSettings.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"

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
  _gridWithAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionAdjustmentAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndPlotAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStepAndPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithFaceDataExchange(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictor(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorRerun(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrector(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),

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
  _gridWithAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionAdjustmentAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndPlotAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStepAndPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithFaceDataExchange(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictor(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorRerun(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrector(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),

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
  
  logInfo( "restart(...)", "start node for subdomain " << domainOffset << "x" << domainSize << " on level " << domainLevel << " with master " << tarch::parallel::NodePool::getInstance().getMasterRank() );
  
  assertion( _repositoryState.getAction() == exahype::records::RepositoryState::Terminate );

  _vertexStack.clear();
  _cellStack.clear();

  _gridWithAugmentedAMRGrid.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlotAugmentedAMRGrid.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithSolutionAdjustmentAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictorAndPlotAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictorAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithADERDGTimeStep.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithADERDGTimeStepAndPlot.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlotAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithFaceDataExchange.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictor.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictorRerun.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrector.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlot.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);


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

  _gridWithAugmentedAMRGrid.terminate();
  _gridWithPlotAugmentedAMRGrid.terminate();
  _gridWithSolutionAdjustmentAndGlobalTimeStepComputation.terminate();
  _gridWithPredictorAndPlotAndGlobalTimeStepComputation.terminate();
  _gridWithPredictorAndGlobalTimeStepComputation.terminate();
  _gridWithADERDGTimeStep.terminate();
  _gridWithADERDGTimeStepAndPlot.terminate();
  _gridWithGlobalTimeStepComputation.terminate();
  _gridWithPlotAndGlobalTimeStepComputation.terminate();
  _gridWithFaceDataExchange.terminate();
  _gridWithPredictor.terminate();
  _gridWithPredictorRerun.terminate();
  _gridWithCorrector.terminate();
  _gridWithPlot.terminate();


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
      case exahype::records::RepositoryState::UseAdapterAugmentedAMRGrid: watch.startTimer(); _gridWithAugmentedAMRGrid.iterate(); watch.stopTimer(); _measureAugmentedAMRGridCPUTime.setValue( watch.getCPUTime() ); _measureAugmentedAMRGridCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid: watch.startTimer(); _gridWithPlotAugmentedAMRGrid.iterate(); watch.stopTimer(); _measurePlotAugmentedAMRGridCPUTime.setValue( watch.getCPUTime() ); _measurePlotAugmentedAMRGridCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterSolutionAdjustmentAndGlobalTimeStepComputation: watch.startTimer(); _gridWithSolutionAdjustmentAndGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measureSolutionAdjustmentAndGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictorAndPlotAndGlobalTimeStepComputation: watch.startTimer(); _gridWithPredictorAndPlotAndGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictorAndGlobalTimeStepComputation: watch.startTimer(); _gridWithPredictorAndGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measurePredictorAndGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measurePredictorAndGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterADERDGTimeStep: watch.startTimer(); _gridWithADERDGTimeStep.iterate(); watch.stopTimer(); _measureADERDGTimeStepCPUTime.setValue( watch.getCPUTime() ); _measureADERDGTimeStepCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterADERDGTimeStepAndPlot: watch.startTimer(); _gridWithADERDGTimeStepAndPlot.iterate(); watch.stopTimer(); _measureADERDGTimeStepAndPlotCPUTime.setValue( watch.getCPUTime() ); _measureADERDGTimeStepAndPlotCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterGlobalTimeStepComputation: watch.startTimer(); _gridWithGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measureGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measureGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPlotAndGlobalTimeStepComputation: watch.startTimer(); _gridWithPlotAndGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measurePlotAndGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measurePlotAndGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterFaceDataExchange: watch.startTimer(); _gridWithFaceDataExchange.iterate(); watch.stopTimer(); _measureFaceDataExchangeCPUTime.setValue( watch.getCPUTime() ); _measureFaceDataExchangeCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictor: watch.startTimer(); _gridWithPredictor.iterate(); watch.stopTimer(); _measurePredictorCPUTime.setValue( watch.getCPUTime() ); _measurePredictorCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictorRerun: watch.startTimer(); _gridWithPredictorRerun.iterate(); watch.stopTimer(); _measurePredictorRerunCPUTime.setValue( watch.getCPUTime() ); _measurePredictorRerunCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrector: watch.startTimer(); _gridWithCorrector.iterate(); watch.stopTimer(); _measureCorrectorCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPlot: watch.startTimer(); _gridWithPlot.iterate(); watch.stopTimer(); _measurePlotCPUTime.setValue( watch.getCPUTime() ); _measurePlotCalendarTime.setValue( watch.getCalendarTime() ); break;

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

 void exahype::repositories::RepositorySTDStack::switchToAugmentedAMRGrid() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterAugmentedAMRGrid); }
 void exahype::repositories::RepositorySTDStack::switchToPlotAugmentedAMRGrid() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid); }
 void exahype::repositories::RepositorySTDStack::switchToSolutionAdjustmentAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterSolutionAdjustmentAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToPredictorAndPlotAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictorAndPlotAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToPredictorAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictorAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToADERDGTimeStep() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterADERDGTimeStep); }
 void exahype::repositories::RepositorySTDStack::switchToADERDGTimeStepAndPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterADERDGTimeStepAndPlot); }
 void exahype::repositories::RepositorySTDStack::switchToGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToPlotAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlotAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToFaceDataExchange() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterFaceDataExchange); }
 void exahype::repositories::RepositorySTDStack::switchToPredictor() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictor); }
 void exahype::repositories::RepositorySTDStack::switchToPredictorRerun() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictorRerun); }
 void exahype::repositories::RepositorySTDStack::switchToCorrector() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrector); }
 void exahype::repositories::RepositorySTDStack::switchToPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlot); }



 bool exahype::repositories::RepositorySTDStack::isActiveAdapterAugmentedAMRGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterAugmentedAMRGrid; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPlotAugmentedAMRGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterSolutionAdjustmentAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterSolutionAdjustmentAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictorAndPlotAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictorAndPlotAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictorAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictorAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterADERDGTimeStep() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterADERDGTimeStep; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterADERDGTimeStepAndPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterADERDGTimeStepAndPlot; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPlotAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlotAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterFaceDataExchange() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterFaceDataExchange; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictor() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictor; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictorRerun() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictorRerun; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrector() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrector; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlot; }



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


void exahype::repositories::RepositorySTDStack::logIterationStatistics(bool logAllAdapters) const {
  logInfo( "logIterationStatistics()", "|| adapter name \t || iterations \t || total CPU time [t]=s \t || average CPU time [t]=s \t || total user time [t]=s \t || average user time [t]=s  || CPU time properties  || user time properties " );  
   if (logAllAdapters || _measureAugmentedAMRGridCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| AugmentedAMRGrid \t |  " << _measureAugmentedAMRGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measureAugmentedAMRGridCPUTime.getAccumulatedValue() << " \t |  " << _measureAugmentedAMRGridCPUTime.getValue()  << " \t |  " << _measureAugmentedAMRGridCalendarTime.getAccumulatedValue() << " \t |  " << _measureAugmentedAMRGridCalendarTime.getValue() << " \t |  " << _measureAugmentedAMRGridCPUTime.toString() << " \t |  " << _measureAugmentedAMRGridCalendarTime.toString() );
   if (logAllAdapters || _measurePlotAugmentedAMRGridCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PlotAugmentedAMRGrid \t |  " << _measurePlotAugmentedAMRGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.getValue()  << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.getValue() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.toString() << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.toString() );
   if (logAllAdapters || _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| SolutionAdjustmentAndGlobalTimeStepComputation \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictorAndPlotAndGlobalTimeStepComputation \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictorAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictorAndGlobalTimeStepComputation \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCalendarTime.toString() );
   if (logAllAdapters || _measureADERDGTimeStepCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| ADERDGTimeStep \t |  " << _measureADERDGTimeStepCPUTime.getNumberOfMeasurements() << " \t |  " << _measureADERDGTimeStepCPUTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepCPUTime.getValue()  << " \t |  " << _measureADERDGTimeStepCalendarTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepCalendarTime.getValue() << " \t |  " << _measureADERDGTimeStepCPUTime.toString() << " \t |  " << _measureADERDGTimeStepCalendarTime.toString() );
   if (logAllAdapters || _measureADERDGTimeStepAndPlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| ADERDGTimeStepAndPlot \t |  " << _measureADERDGTimeStepAndPlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measureADERDGTimeStepAndPlotCPUTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepAndPlotCPUTime.getValue()  << " \t |  " << _measureADERDGTimeStepAndPlotCalendarTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepAndPlotCalendarTime.getValue() << " \t |  " << _measureADERDGTimeStepAndPlotCPUTime.toString() << " \t |  " << _measureADERDGTimeStepAndPlotCalendarTime.toString() );
   if (logAllAdapters || _measureGlobalTimeStepComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| GlobalTimeStepComputation \t |  " << _measureGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measureGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measureGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measureGlobalTimeStepComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePlotAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PlotAndGlobalTimeStepComputation \t |  " << _measurePlotAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measurePlotAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measurePlotAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measurePlotAndGlobalTimeStepComputationCalendarTime.toString() );
   if (logAllAdapters || _measureFaceDataExchangeCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| FaceDataExchange \t |  " << _measureFaceDataExchangeCPUTime.getNumberOfMeasurements() << " \t |  " << _measureFaceDataExchangeCPUTime.getAccumulatedValue() << " \t |  " << _measureFaceDataExchangeCPUTime.getValue()  << " \t |  " << _measureFaceDataExchangeCalendarTime.getAccumulatedValue() << " \t |  " << _measureFaceDataExchangeCalendarTime.getValue() << " \t |  " << _measureFaceDataExchangeCPUTime.toString() << " \t |  " << _measureFaceDataExchangeCalendarTime.toString() );
   if (logAllAdapters || _measurePredictorCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Predictor \t |  " << _measurePredictorCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorCPUTime.getValue()  << " \t |  " << _measurePredictorCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorCalendarTime.getValue() << " \t |  " << _measurePredictorCPUTime.toString() << " \t |  " << _measurePredictorCalendarTime.toString() );
   if (logAllAdapters || _measurePredictorRerunCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictorRerun \t |  " << _measurePredictorRerunCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorRerunCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorRerunCPUTime.getValue()  << " \t |  " << _measurePredictorRerunCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorRerunCalendarTime.getValue() << " \t |  " << _measurePredictorRerunCPUTime.toString() << " \t |  " << _measurePredictorRerunCalendarTime.toString() );
   if (logAllAdapters || _measureCorrectorCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Corrector \t |  " << _measureCorrectorCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorCPUTime.getValue()  << " \t |  " << _measureCorrectorCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorCalendarTime.getValue() << " \t |  " << _measureCorrectorCPUTime.toString() << " \t |  " << _measureCorrectorCalendarTime.toString() );
   if (logAllAdapters || _measurePlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Plot \t |  " << _measurePlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotCPUTime.getValue()  << " \t |  " << _measurePlotCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotCalendarTime.getValue() << " \t |  " << _measurePlotCPUTime.toString() << " \t |  " << _measurePlotCalendarTime.toString() );

}


void exahype::repositories::RepositorySTDStack::clearIterationStatistics() {
   _measureAugmentedAMRGridCPUTime.erase();
   _measurePlotAugmentedAMRGridCPUTime.erase();
   _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.erase();
   _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.erase();
   _measurePredictorAndGlobalTimeStepComputationCPUTime.erase();
   _measureADERDGTimeStepCPUTime.erase();
   _measureADERDGTimeStepAndPlotCPUTime.erase();
   _measureGlobalTimeStepComputationCPUTime.erase();
   _measurePlotAndGlobalTimeStepComputationCPUTime.erase();
   _measureFaceDataExchangeCPUTime.erase();
   _measurePredictorCPUTime.erase();
   _measurePredictorRerunCPUTime.erase();
   _measureCorrectorCPUTime.erase();
   _measurePlotCPUTime.erase();

   _measureAugmentedAMRGridCalendarTime.erase();
   _measurePlotAugmentedAMRGridCalendarTime.erase();
   _measureSolutionAdjustmentAndGlobalTimeStepComputationCalendarTime.erase();
   _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.erase();
   _measurePredictorAndGlobalTimeStepComputationCalendarTime.erase();
   _measureADERDGTimeStepCalendarTime.erase();
   _measureADERDGTimeStepAndPlotCalendarTime.erase();
   _measureGlobalTimeStepComputationCalendarTime.erase();
   _measurePlotAndGlobalTimeStepComputationCalendarTime.erase();
   _measureFaceDataExchangeCalendarTime.erase();
   _measurePredictorCalendarTime.erase();
   _measurePredictorRerunCalendarTime.erase();
   _measureCorrectorCalendarTime.erase();
   _measurePlotCalendarTime.erase();

}
