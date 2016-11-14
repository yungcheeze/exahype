#include "exahype/repositories/RepositorySTDStack.h"

#include "tarch/Assertions.h"
#include "tarch/timing/Watch.h"

#include "tarch/compiler/CompilerSpecificSettings.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"

#ifdef Parallel
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
  _gridWithMeshRefinement(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndPlotAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndPlotAndTimeStepSizeComputation2d(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGridErasing(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAndADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionRerun(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithLimiterStatusSpreading(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithReinitialisation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionRecomputationAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithNeighbourDataMerging(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPrediction(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionUpdate(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAndSolutionUpdate(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
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
  _gridWithMeshRefinement(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndPlotAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndPlotAndTimeStepSizeComputation2d(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGridErasing(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAndADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionRerun(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithLimiterStatusSpreading(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithReinitialisation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionRecomputationAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithNeighbourDataMerging(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPrediction(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionUpdate(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAndSolutionUpdate(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
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

  _gridWithMeshRefinement.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlotAugmentedAMRGrid.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithInitialConditionAndTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionAndPlotAndTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionAndPlotAndTimeStepSizeComputation2d.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionAndTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithGridErasing.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithADERDGTimeStep.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlotAndADERDGTimeStep.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionRerun.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithLimiterStatusSpreading.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithReinitialisation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithSolutionRecomputationAndTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithNeighbourDataMerging.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPrediction.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithSolutionUpdate.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlotAndSolutionUpdate.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
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

  _gridWithMeshRefinement.terminate();
  _gridWithPlotAugmentedAMRGrid.terminate();
  _gridWithInitialConditionAndTimeStepSizeComputation.terminate();
  _gridWithPredictionAndPlotAndTimeStepSizeComputation.terminate();
  _gridWithPredictionAndPlotAndTimeStepSizeComputation2d.terminate();
  _gridWithPredictionAndTimeStepSizeComputation.terminate();
  _gridWithGridErasing.terminate();
  _gridWithADERDGTimeStep.terminate();
  _gridWithPlotAndADERDGTimeStep.terminate();
  _gridWithPredictionRerun.terminate();
  _gridWithLimiterStatusSpreading.terminate();
  _gridWithReinitialisation.terminate();
  _gridWithSolutionRecomputationAndTimeStepSizeComputation.terminate();
  _gridWithNeighbourDataMerging.terminate();
  _gridWithPrediction.terminate();
  _gridWithSolutionUpdate.terminate();
  _gridWithPlotAndSolutionUpdate.terminate();
  _gridWithPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation.terminate();
  _gridWithTimeStepSizeComputation.terminate();
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

  if ( numberOfIterations > 1 && _solverState.isInvolvedInJoinOrFork() ) {
    logWarning( "iterate()", "iterate invoked for multiple traversals though load balancing still does redistribute data" );
  }
  bool switchedLoadBalancingTemporarilyOff = false;
  if ( numberOfIterations > 1 && peano::parallel::loadbalancing::Oracle::getInstance().isLoadBalancingActivated() ) {
    switchedLoadBalancingTemporarilyOff = true;
    peano::parallel::loadbalancing::Oracle::getInstance().activateLoadBalancing(false);
  }

  peano::datatraversal::autotuning::Oracle::getInstance().switchToOracle(_repositoryState.getAction());

  peano::parallel::loadbalancing::Oracle::getInstance().switchToOracle(_repositoryState.getAction());
  
  _solverState.currentlyRunsMultipleIterations(_repositoryState.getNumberOfIterations()>1);
  #else
  peano::datatraversal::autotuning::Oracle::getInstance().switchToOracle(_repositoryState.getAction());
  #endif
  
  for (int i=0; i<numberOfIterations; i++) {
    switch ( _repositoryState.getAction()) {
      case exahype::records::RepositoryState::UseAdapterMeshRefinement: watch.startTimer(); _gridWithMeshRefinement.iterate(); watch.stopTimer(); _measureMeshRefinementCPUTime.setValue( watch.getCPUTime() ); _measureMeshRefinementCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid: watch.startTimer(); _gridWithPlotAugmentedAMRGrid.iterate(); watch.stopTimer(); _measurePlotAugmentedAMRGridCPUTime.setValue( watch.getCPUTime() ); _measurePlotAugmentedAMRGridCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterInitialConditionAndTimeStepSizeComputation: watch.startTimer(); _gridWithInitialConditionAndTimeStepSizeComputation.iterate(); watch.stopTimer(); _measureInitialConditionAndTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measureInitialConditionAndTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionAndPlotAndTimeStepSizeComputation: watch.startTimer(); _gridWithPredictionAndPlotAndTimeStepSizeComputation.iterate(); watch.stopTimer(); _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measurePredictionAndPlotAndTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionAndPlotAndTimeStepSizeComputation2d: watch.startTimer(); _gridWithPredictionAndPlotAndTimeStepSizeComputation2d.iterate(); watch.stopTimer(); _measurePredictionAndPlotAndTimeStepSizeComputation2dCPUTime.setValue( watch.getCPUTime() ); _measurePredictionAndPlotAndTimeStepSizeComputation2dCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionAndTimeStepSizeComputation: watch.startTimer(); _gridWithPredictionAndTimeStepSizeComputation.iterate(); watch.stopTimer(); _measurePredictionAndTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measurePredictionAndTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterGridErasing: watch.startTimer(); _gridWithGridErasing.iterate(); watch.stopTimer(); _measureGridErasingCPUTime.setValue( watch.getCPUTime() ); _measureGridErasingCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterADERDGTimeStep: watch.startTimer(); _gridWithADERDGTimeStep.iterate(); watch.stopTimer(); _measureADERDGTimeStepCPUTime.setValue( watch.getCPUTime() ); _measureADERDGTimeStepCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPlotAndADERDGTimeStep: watch.startTimer(); _gridWithPlotAndADERDGTimeStep.iterate(); watch.stopTimer(); _measurePlotAndADERDGTimeStepCPUTime.setValue( watch.getCPUTime() ); _measurePlotAndADERDGTimeStepCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionRerun: watch.startTimer(); _gridWithPredictionRerun.iterate(); watch.stopTimer(); _measurePredictionRerunCPUTime.setValue( watch.getCPUTime() ); _measurePredictionRerunCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterLimiterStatusSpreading: watch.startTimer(); _gridWithLimiterStatusSpreading.iterate(); watch.stopTimer(); _measureLimiterStatusSpreadingCPUTime.setValue( watch.getCPUTime() ); _measureLimiterStatusSpreadingCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterReinitialisation: watch.startTimer(); _gridWithReinitialisation.iterate(); watch.stopTimer(); _measureReinitialisationCPUTime.setValue( watch.getCPUTime() ); _measureReinitialisationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterSolutionRecomputationAndTimeStepSizeComputation: watch.startTimer(); _gridWithSolutionRecomputationAndTimeStepSizeComputation.iterate(); watch.stopTimer(); _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measureSolutionRecomputationAndTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterNeighbourDataMerging: watch.startTimer(); _gridWithNeighbourDataMerging.iterate(); watch.stopTimer(); _measureNeighbourDataMergingCPUTime.setValue( watch.getCPUTime() ); _measureNeighbourDataMergingCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPrediction: watch.startTimer(); _gridWithPrediction.iterate(); watch.stopTimer(); _measurePredictionCPUTime.setValue( watch.getCPUTime() ); _measurePredictionCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterSolutionUpdate: watch.startTimer(); _gridWithSolutionUpdate.iterate(); watch.stopTimer(); _measureSolutionUpdateCPUTime.setValue( watch.getCPUTime() ); _measureSolutionUpdateCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPlotAndSolutionUpdate: watch.startTimer(); _gridWithPlotAndSolutionUpdate.iterate(); watch.stopTimer(); _measurePlotAndSolutionUpdateCPUTime.setValue( watch.getCPUTime() ); _measurePlotAndSolutionUpdateCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation: watch.startTimer(); _gridWithPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation.iterate(); watch.stopTimer(); _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterTimeStepSizeComputation: watch.startTimer(); _gridWithTimeStepSizeComputation.iterate(); watch.stopTimer(); _measureTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measureTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
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
    #ifdef Parallel
    if ( switchedLoadBalancingTemporarilyOff && i==numberOfIterations-1) {
      peano::parallel::loadbalancing::Oracle::getInstance().activateLoadBalancing(true);
    }
    #endif
  }
  
  #ifdef Parallel
  if (_solverState.isJoiningWithMaster()) {
    _repositoryState.setAction( exahype::records::RepositoryState::Terminate );
  }
  #endif
}

 void exahype::repositories::RepositorySTDStack::switchToMeshRefinement() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterMeshRefinement); }
 void exahype::repositories::RepositorySTDStack::switchToPlotAugmentedAMRGrid() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid); }
 void exahype::repositories::RepositorySTDStack::switchToInitialConditionAndTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterInitialConditionAndTimeStepSizeComputation); }
 void exahype::repositories::RepositorySTDStack::switchToPredictionAndPlotAndTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionAndPlotAndTimeStepSizeComputation); }
 void exahype::repositories::RepositorySTDStack::switchToPredictionAndPlotAndTimeStepSizeComputation2d() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionAndPlotAndTimeStepSizeComputation2d); }
 void exahype::repositories::RepositorySTDStack::switchToPredictionAndTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionAndTimeStepSizeComputation); }
 void exahype::repositories::RepositorySTDStack::switchToGridErasing() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterGridErasing); }
 void exahype::repositories::RepositorySTDStack::switchToADERDGTimeStep() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterADERDGTimeStep); }
 void exahype::repositories::RepositorySTDStack::switchToPlotAndADERDGTimeStep() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlotAndADERDGTimeStep); }
 void exahype::repositories::RepositorySTDStack::switchToPredictionRerun() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionRerun); }
 void exahype::repositories::RepositorySTDStack::switchToLimiterStatusSpreading() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterLimiterStatusSpreading); }
 void exahype::repositories::RepositorySTDStack::switchToReinitialisation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterReinitialisation); }
 void exahype::repositories::RepositorySTDStack::switchToSolutionRecomputationAndTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterSolutionRecomputationAndTimeStepSizeComputation); }
 void exahype::repositories::RepositorySTDStack::switchToNeighbourDataMerging() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterNeighbourDataMerging); }
 void exahype::repositories::RepositorySTDStack::switchToPrediction() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPrediction); }
 void exahype::repositories::RepositorySTDStack::switchToSolutionUpdate() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterSolutionUpdate); }
 void exahype::repositories::RepositorySTDStack::switchToPlotAndSolutionUpdate() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlotAndSolutionUpdate); }
 void exahype::repositories::RepositorySTDStack::switchToPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation); }
 void exahype::repositories::RepositorySTDStack::switchToTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterTimeStepSizeComputation); }
 void exahype::repositories::RepositorySTDStack::switchToPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlot); }



 bool exahype::repositories::RepositorySTDStack::isActiveAdapterMeshRefinement() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterMeshRefinement; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPlotAugmentedAMRGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterInitialConditionAndTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterInitialConditionAndTimeStepSizeComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictionAndPlotAndTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionAndPlotAndTimeStepSizeComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictionAndPlotAndTimeStepSizeComputation2d() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionAndPlotAndTimeStepSizeComputation2d; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictionAndTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionAndTimeStepSizeComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterGridErasing() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterGridErasing; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterADERDGTimeStep() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterADERDGTimeStep; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPlotAndADERDGTimeStep() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlotAndADERDGTimeStep; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictionRerun() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionRerun; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterLimiterStatusSpreading() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterLimiterStatusSpreading; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterReinitialisation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterReinitialisation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterSolutionRecomputationAndTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterSolutionRecomputationAndTimeStepSizeComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterNeighbourDataMerging() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterNeighbourDataMerging; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPrediction() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPrediction; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterSolutionUpdate() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterSolutionUpdate; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPlotAndSolutionUpdate() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlotAndSolutionUpdate; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterTimeStepSizeComputation; }
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
   if (logAllAdapters || _measureMeshRefinementCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| MeshRefinement \t |  " << _measureMeshRefinementCPUTime.getNumberOfMeasurements() << " \t |  " << _measureMeshRefinementCPUTime.getAccumulatedValue() << " \t |  " << _measureMeshRefinementCPUTime.getValue()  << " \t |  " << _measureMeshRefinementCalendarTime.getAccumulatedValue() << " \t |  " << _measureMeshRefinementCalendarTime.getValue() << " \t |  " << _measureMeshRefinementCPUTime.toString() << " \t |  " << _measureMeshRefinementCalendarTime.toString() );
   if (logAllAdapters || _measurePlotAugmentedAMRGridCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PlotAugmentedAMRGrid \t |  " << _measurePlotAugmentedAMRGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.getValue()  << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.getValue() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.toString() << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.toString() );
   if (logAllAdapters || _measureInitialConditionAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| InitialConditionAndTimeStepSizeComputation \t |  " << _measureInitialConditionAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionAndPlotAndTimeStepSizeComputation \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionAndPlotAndTimeStepSizeComputation2dCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionAndPlotAndTimeStepSizeComputation2d \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputation2dCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputation2dCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputation2dCPUTime.getValue()  << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputation2dCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputation2dCalendarTime.getValue() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputation2dCPUTime.toString() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputation2dCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionAndTimeStepSizeComputation \t |  " << _measurePredictionAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionAndTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measurePredictionAndTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measurePredictionAndTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measurePredictionAndTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measureGridErasingCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| GridErasing \t |  " << _measureGridErasingCPUTime.getNumberOfMeasurements() << " \t |  " << _measureGridErasingCPUTime.getAccumulatedValue() << " \t |  " << _measureGridErasingCPUTime.getValue()  << " \t |  " << _measureGridErasingCalendarTime.getAccumulatedValue() << " \t |  " << _measureGridErasingCalendarTime.getValue() << " \t |  " << _measureGridErasingCPUTime.toString() << " \t |  " << _measureGridErasingCalendarTime.toString() );
   if (logAllAdapters || _measureADERDGTimeStepCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| ADERDGTimeStep \t |  " << _measureADERDGTimeStepCPUTime.getNumberOfMeasurements() << " \t |  " << _measureADERDGTimeStepCPUTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepCPUTime.getValue()  << " \t |  " << _measureADERDGTimeStepCalendarTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepCalendarTime.getValue() << " \t |  " << _measureADERDGTimeStepCPUTime.toString() << " \t |  " << _measureADERDGTimeStepCalendarTime.toString() );
   if (logAllAdapters || _measurePlotAndADERDGTimeStepCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PlotAndADERDGTimeStep \t |  " << _measurePlotAndADERDGTimeStepCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotAndADERDGTimeStepCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotAndADERDGTimeStepCPUTime.getValue()  << " \t |  " << _measurePlotAndADERDGTimeStepCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotAndADERDGTimeStepCalendarTime.getValue() << " \t |  " << _measurePlotAndADERDGTimeStepCPUTime.toString() << " \t |  " << _measurePlotAndADERDGTimeStepCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionRerunCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionRerun \t |  " << _measurePredictionRerunCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionRerunCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionRerunCPUTime.getValue()  << " \t |  " << _measurePredictionRerunCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionRerunCalendarTime.getValue() << " \t |  " << _measurePredictionRerunCPUTime.toString() << " \t |  " << _measurePredictionRerunCalendarTime.toString() );
   if (logAllAdapters || _measureLimiterStatusSpreadingCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| LimiterStatusSpreading \t |  " << _measureLimiterStatusSpreadingCPUTime.getNumberOfMeasurements() << " \t |  " << _measureLimiterStatusSpreadingCPUTime.getAccumulatedValue() << " \t |  " << _measureLimiterStatusSpreadingCPUTime.getValue()  << " \t |  " << _measureLimiterStatusSpreadingCalendarTime.getAccumulatedValue() << " \t |  " << _measureLimiterStatusSpreadingCalendarTime.getValue() << " \t |  " << _measureLimiterStatusSpreadingCPUTime.toString() << " \t |  " << _measureLimiterStatusSpreadingCalendarTime.toString() );
   if (logAllAdapters || _measureReinitialisationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Reinitialisation \t |  " << _measureReinitialisationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureReinitialisationCPUTime.getAccumulatedValue() << " \t |  " << _measureReinitialisationCPUTime.getValue()  << " \t |  " << _measureReinitialisationCalendarTime.getAccumulatedValue() << " \t |  " << _measureReinitialisationCalendarTime.getValue() << " \t |  " << _measureReinitialisationCPUTime.toString() << " \t |  " << _measureReinitialisationCalendarTime.toString() );
   if (logAllAdapters || _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| SolutionRecomputationAndTimeStepSizeComputation \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measureNeighbourDataMergingCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| NeighbourDataMerging \t |  " << _measureNeighbourDataMergingCPUTime.getNumberOfMeasurements() << " \t |  " << _measureNeighbourDataMergingCPUTime.getAccumulatedValue() << " \t |  " << _measureNeighbourDataMergingCPUTime.getValue()  << " \t |  " << _measureNeighbourDataMergingCalendarTime.getAccumulatedValue() << " \t |  " << _measureNeighbourDataMergingCalendarTime.getValue() << " \t |  " << _measureNeighbourDataMergingCPUTime.toString() << " \t |  " << _measureNeighbourDataMergingCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Prediction \t |  " << _measurePredictionCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionCPUTime.getValue()  << " \t |  " << _measurePredictionCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionCalendarTime.getValue() << " \t |  " << _measurePredictionCPUTime.toString() << " \t |  " << _measurePredictionCalendarTime.toString() );
   if (logAllAdapters || _measureSolutionUpdateCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| SolutionUpdate \t |  " << _measureSolutionUpdateCPUTime.getNumberOfMeasurements() << " \t |  " << _measureSolutionUpdateCPUTime.getAccumulatedValue() << " \t |  " << _measureSolutionUpdateCPUTime.getValue()  << " \t |  " << _measureSolutionUpdateCalendarTime.getAccumulatedValue() << " \t |  " << _measureSolutionUpdateCalendarTime.getValue() << " \t |  " << _measureSolutionUpdateCPUTime.toString() << " \t |  " << _measureSolutionUpdateCalendarTime.toString() );
   if (logAllAdapters || _measurePlotAndSolutionUpdateCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PlotAndSolutionUpdate \t |  " << _measurePlotAndSolutionUpdateCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotAndSolutionUpdateCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotAndSolutionUpdateCPUTime.getValue()  << " \t |  " << _measurePlotAndSolutionUpdateCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotAndSolutionUpdateCalendarTime.getValue() << " \t |  " << _measurePlotAndSolutionUpdateCPUTime.toString() << " \t |  " << _measurePlotAndSolutionUpdateCalendarTime.toString() );
   if (logAllAdapters || _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measureTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| TimeStepSizeComputation \t |  " << _measureTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measureTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measureTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measureTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Plot \t |  " << _measurePlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotCPUTime.getValue()  << " \t |  " << _measurePlotCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotCalendarTime.getValue() << " \t |  " << _measurePlotCPUTime.toString() << " \t |  " << _measurePlotCalendarTime.toString() );

}


void exahype::repositories::RepositorySTDStack::clearIterationStatistics() {
   _measureMeshRefinementCPUTime.erase();
   _measurePlotAugmentedAMRGridCPUTime.erase();
   _measureInitialConditionAndTimeStepSizeComputationCPUTime.erase();
   _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.erase();
   _measurePredictionAndPlotAndTimeStepSizeComputation2dCPUTime.erase();
   _measurePredictionAndTimeStepSizeComputationCPUTime.erase();
   _measureGridErasingCPUTime.erase();
   _measureADERDGTimeStepCPUTime.erase();
   _measurePlotAndADERDGTimeStepCPUTime.erase();
   _measurePredictionRerunCPUTime.erase();
   _measureLimiterStatusSpreadingCPUTime.erase();
   _measureReinitialisationCPUTime.erase();
   _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.erase();
   _measureNeighbourDataMergingCPUTime.erase();
   _measurePredictionCPUTime.erase();
   _measureSolutionUpdateCPUTime.erase();
   _measurePlotAndSolutionUpdateCPUTime.erase();
   _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.erase();
   _measureTimeStepSizeComputationCPUTime.erase();
   _measurePlotCPUTime.erase();

   _measureMeshRefinementCalendarTime.erase();
   _measurePlotAugmentedAMRGridCalendarTime.erase();
   _measureInitialConditionAndTimeStepSizeComputationCalendarTime.erase();
   _measurePredictionAndPlotAndTimeStepSizeComputationCalendarTime.erase();
   _measurePredictionAndPlotAndTimeStepSizeComputation2dCalendarTime.erase();
   _measurePredictionAndTimeStepSizeComputationCalendarTime.erase();
   _measureGridErasingCalendarTime.erase();
   _measureADERDGTimeStepCalendarTime.erase();
   _measurePlotAndADERDGTimeStepCalendarTime.erase();
   _measurePredictionRerunCalendarTime.erase();
   _measureLimiterStatusSpreadingCalendarTime.erase();
   _measureReinitialisationCalendarTime.erase();
   _measureSolutionRecomputationAndTimeStepSizeComputationCalendarTime.erase();
   _measureNeighbourDataMergingCalendarTime.erase();
   _measurePredictionCalendarTime.erase();
   _measureSolutionUpdateCalendarTime.erase();
   _measurePlotAndSolutionUpdateCalendarTime.erase();
   _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCalendarTime.erase();
   _measureTimeStepSizeComputationCalendarTime.erase();
   _measurePlotCalendarTime.erase();

}
