#include "exahype/repositories/RepositoryArrayStack.h"

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


tarch::logging::Log exahype::repositories::RepositoryArrayStack::_log( "exahype::repositories::RepositoryArrayStack" );


exahype::repositories::RepositoryArrayStack::RepositoryArrayStack(
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
  _gridWithAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionAdjustmentAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndPlotAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStepAndPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorRerun(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithRiemannSolver(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictor(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrector(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositoryArrayStack(...)" );
  
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositoryArrayStack(...)" );
}



exahype::repositories::RepositoryArrayStack::RepositoryArrayStack(
  peano::geometry::Geometry&                   geometry,
  int                                          maximumSizeOfCellInOutStack,
  int                                          maximumSizeOfVertexInOutStack,
  int                                          maximumSizeOfVertexTemporaryStack
):
  _geometry(geometry),
  _cellStack(maximumSizeOfCellInOutStack),
  _vertexStack(maximumSizeOfVertexInOutStack,maximumSizeOfVertexTemporaryStack),
  _solverState(),
  _gridWithAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionAdjustmentAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndPlotAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStepAndPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorRerun(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithRiemannSolver(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictor(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrector(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositoryArrayStack(Geometry&)" );
  
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositoryArrayStack(Geometry&)" );
}
    
   
exahype::repositories::RepositoryArrayStack::~RepositoryArrayStack() {
  assertion( _repositoryState.getAction() == exahype::records::RepositoryState::Terminate );
}


void exahype::repositories::RepositoryArrayStack::restart(
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

  _gridWithAugmentedAMRGrid.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlotAugmentedAMRGrid.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithSolutionAdjustmentAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictorAndPlotAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictorAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithADERDGTimeStep.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithADERDGTimeStepAndPlot.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictorRerun.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithRiemannSolver.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictor.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrector.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrectorAndPlot.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlot.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);

 
   _solverState.restart();
 
  logTraceOut( "restart(...)" );
}


void exahype::repositories::RepositoryArrayStack::terminate() {
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
  _gridWithPredictorRerun.terminate();
  _gridWithRiemannSolver.terminate();
  _gridWithPredictor.terminate();
  _gridWithCorrector.terminate();
  _gridWithCorrectorAndPlot.terminate();
  _gridWithPlot.terminate();

 
  logTraceOut( "terminate()" );
}


exahype::State& exahype::repositories::RepositoryArrayStack::getState() {
  return _solverState;
}


const exahype::State& exahype::repositories::RepositoryArrayStack::getState() const {
  return _solverState;
}

   
void exahype::repositories::RepositoryArrayStack::iterate(int numberOfIterations, bool exchangeBoundaryVertices) {
  tarch::timing::Watch watch( "exahype::repositories::RepositoryArrayStack", "iterate(bool)", false);
  
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
      case exahype::records::RepositoryState::UseAdapterPredictorRerun: watch.startTimer(); _gridWithPredictorRerun.iterate(); watch.stopTimer(); _measurePredictorRerunCPUTime.setValue( watch.getCPUTime() ); _measurePredictorRerunCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterRiemannSolver: watch.startTimer(); _gridWithRiemannSolver.iterate(); watch.stopTimer(); _measureRiemannSolverCPUTime.setValue( watch.getCPUTime() ); _measureRiemannSolverCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictor: watch.startTimer(); _gridWithPredictor.iterate(); watch.stopTimer(); _measurePredictorCPUTime.setValue( watch.getCPUTime() ); _measurePredictorCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrector: watch.startTimer(); _gridWithCorrector.iterate(); watch.stopTimer(); _measureCorrectorCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrectorAndPlot: watch.startTimer(); _gridWithCorrectorAndPlot.iterate(); watch.stopTimer(); _measureCorrectorAndPlotCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorAndPlotCalendarTime.setValue( watch.getCalendarTime() ); break;
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

 void exahype::repositories::RepositoryArrayStack::switchToAugmentedAMRGrid() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterAugmentedAMRGrid); }
 void exahype::repositories::RepositoryArrayStack::switchToPlotAugmentedAMRGrid() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid); }
 void exahype::repositories::RepositoryArrayStack::switchToSolutionAdjustmentAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterSolutionAdjustmentAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictorAndPlotAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictorAndPlotAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictorAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictorAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositoryArrayStack::switchToADERDGTimeStep() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterADERDGTimeStep); }
 void exahype::repositories::RepositoryArrayStack::switchToADERDGTimeStepAndPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterADERDGTimeStepAndPlot); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictorRerun() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictorRerun); }
 void exahype::repositories::RepositoryArrayStack::switchToRiemannSolver() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterRiemannSolver); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictor() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictor); }
 void exahype::repositories::RepositoryArrayStack::switchToCorrector() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrector); }
 void exahype::repositories::RepositoryArrayStack::switchToCorrectorAndPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrectorAndPlot); }
 void exahype::repositories::RepositoryArrayStack::switchToPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlot); }



 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterAugmentedAMRGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterAugmentedAMRGrid; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPlotAugmentedAMRGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterSolutionAdjustmentAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterSolutionAdjustmentAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictorAndPlotAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictorAndPlotAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictorAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictorAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterADERDGTimeStep() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterADERDGTimeStep; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterADERDGTimeStepAndPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterADERDGTimeStepAndPlot; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictorRerun() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictorRerun; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterRiemannSolver() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterRiemannSolver; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictor() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictor; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterCorrector() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrector; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterCorrectorAndPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrectorAndPlot; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlot; }



peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>* exahype::repositories::RepositoryArrayStack::createEmptyCheckpoint() {
  return new peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>();
} 


void exahype::repositories::RepositoryArrayStack::writeCheckpoint(peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> * const checkpoint) {
  _solverState.writeToCheckpoint( *checkpoint );
  _vertexStack.writeToCheckpoint( *checkpoint );
  _cellStack.writeToCheckpoint( *checkpoint );
} 


void exahype::repositories::RepositoryArrayStack::setMaximumMemoryFootprintForTemporaryRegularGrids(double value) {
  _regularGridContainer.setMaximumMemoryFootprintForTemporaryRegularGrids(value);
}


void exahype::repositories::RepositoryArrayStack::readCheckpoint( peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> const * const checkpoint ) {
  assertionMsg( checkpoint->isValid(), "checkpoint has to be valid if you call this operation" );

  _solverState.readFromCheckpoint( *checkpoint );
  _vertexStack.readFromCheckpoint( *checkpoint );
  _cellStack.readFromCheckpoint( *checkpoint );
}


#ifdef Parallel
void exahype::repositories::RepositoryArrayStack::runGlobalStep() {
  assertion(tarch::parallel::Node::getInstance().isGlobalMaster());

  exahype::records::RepositoryState intermediateStateForWorkingNodes;
  intermediateStateForWorkingNodes.setAction( exahype::records::RepositoryState::RunOnAllNodes );
  
  tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
    intermediateStateForWorkingNodes,
    peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
  );
  tarch::parallel::NodePool::getInstance().activateIdleNodes();
}


exahype::repositories::RepositoryArrayStack::ContinueCommand exahype::repositories::RepositoryArrayStack::continueToIterate() {
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


void exahype::repositories::RepositoryArrayStack::logIterationStatistics(bool logAllAdapters) const {
  logInfo( "logIterationStatistics()", "|| adapter name \t || iterations \t || total CPU time [t]=s \t || average CPU time [t]=s \t || total user time [t]=s \t || average user time [t]=s  || CPU time properties  || user time properties " );  
   if (logAllAdapters || _measureAugmentedAMRGridCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| AugmentedAMRGrid \t |  " << _measureAugmentedAMRGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measureAugmentedAMRGridCPUTime.getAccumulatedValue() << " \t |  " << _measureAugmentedAMRGridCPUTime.getValue()  << " \t |  " << _measureAugmentedAMRGridCalendarTime.getAccumulatedValue() << " \t |  " << _measureAugmentedAMRGridCalendarTime.getValue() << " \t |  " << _measureAugmentedAMRGridCPUTime.toString() << " \t |  " << _measureAugmentedAMRGridCalendarTime.toString() );
   if (logAllAdapters || _measurePlotAugmentedAMRGridCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PlotAugmentedAMRGrid \t |  " << _measurePlotAugmentedAMRGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.getValue()  << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.getValue() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.toString() << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.toString() );
   if (logAllAdapters || _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| SolutionAdjustmentAndGlobalTimeStepComputation \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measureSolutionAdjustmentAndGlobalTimeStepComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictorAndPlotAndGlobalTimeStepComputation \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictorAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictorAndGlobalTimeStepComputation \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCalendarTime.toString() );
   if (logAllAdapters || _measureADERDGTimeStepCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| ADERDGTimeStep \t |  " << _measureADERDGTimeStepCPUTime.getNumberOfMeasurements() << " \t |  " << _measureADERDGTimeStepCPUTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepCPUTime.getValue()  << " \t |  " << _measureADERDGTimeStepCalendarTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepCalendarTime.getValue() << " \t |  " << _measureADERDGTimeStepCPUTime.toString() << " \t |  " << _measureADERDGTimeStepCalendarTime.toString() );
   if (logAllAdapters || _measureADERDGTimeStepAndPlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| ADERDGTimeStepAndPlot \t |  " << _measureADERDGTimeStepAndPlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measureADERDGTimeStepAndPlotCPUTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepAndPlotCPUTime.getValue()  << " \t |  " << _measureADERDGTimeStepAndPlotCalendarTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepAndPlotCalendarTime.getValue() << " \t |  " << _measureADERDGTimeStepAndPlotCPUTime.toString() << " \t |  " << _measureADERDGTimeStepAndPlotCalendarTime.toString() );
   if (logAllAdapters || _measurePredictorRerunCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictorRerun \t |  " << _measurePredictorRerunCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorRerunCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorRerunCPUTime.getValue()  << " \t |  " << _measurePredictorRerunCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorRerunCalendarTime.getValue() << " \t |  " << _measurePredictorRerunCPUTime.toString() << " \t |  " << _measurePredictorRerunCalendarTime.toString() );
   if (logAllAdapters || _measureRiemannSolverCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| RiemannSolver \t |  " << _measureRiemannSolverCPUTime.getNumberOfMeasurements() << " \t |  " << _measureRiemannSolverCPUTime.getAccumulatedValue() << " \t |  " << _measureRiemannSolverCPUTime.getValue()  << " \t |  " << _measureRiemannSolverCalendarTime.getAccumulatedValue() << " \t |  " << _measureRiemannSolverCalendarTime.getValue() << " \t |  " << _measureRiemannSolverCPUTime.toString() << " \t |  " << _measureRiemannSolverCalendarTime.toString() );
   if (logAllAdapters || _measurePredictorCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Predictor \t |  " << _measurePredictorCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorCPUTime.getValue()  << " \t |  " << _measurePredictorCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorCalendarTime.getValue() << " \t |  " << _measurePredictorCPUTime.toString() << " \t |  " << _measurePredictorCalendarTime.toString() );
   if (logAllAdapters || _measureCorrectorCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Corrector \t |  " << _measureCorrectorCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorCPUTime.getValue()  << " \t |  " << _measureCorrectorCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorCalendarTime.getValue() << " \t |  " << _measureCorrectorCPUTime.toString() << " \t |  " << _measureCorrectorCalendarTime.toString() );
   if (logAllAdapters || _measureCorrectorAndPlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| CorrectorAndPlot \t |  " << _measureCorrectorAndPlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorAndPlotCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPlotCPUTime.getValue()  << " \t |  " << _measureCorrectorAndPlotCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPlotCalendarTime.getValue() << " \t |  " << _measureCorrectorAndPlotCPUTime.toString() << " \t |  " << _measureCorrectorAndPlotCalendarTime.toString() );
   if (logAllAdapters || _measurePlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Plot \t |  " << _measurePlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotCPUTime.getValue()  << " \t |  " << _measurePlotCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotCalendarTime.getValue() << " \t |  " << _measurePlotCPUTime.toString() << " \t |  " << _measurePlotCalendarTime.toString() );

}


void exahype::repositories::RepositoryArrayStack::clearIterationStatistics() {
   _measureAugmentedAMRGridCPUTime.erase();
   _measurePlotAugmentedAMRGridCPUTime.erase();
   _measureSolutionAdjustmentAndGlobalTimeStepComputationCPUTime.erase();
   _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.erase();
   _measurePredictorAndGlobalTimeStepComputationCPUTime.erase();
   _measureADERDGTimeStepCPUTime.erase();
   _measureADERDGTimeStepAndPlotCPUTime.erase();
   _measurePredictorRerunCPUTime.erase();
   _measureRiemannSolverCPUTime.erase();
   _measurePredictorCPUTime.erase();
   _measureCorrectorCPUTime.erase();
   _measureCorrectorAndPlotCPUTime.erase();
   _measurePlotCPUTime.erase();

   _measureAugmentedAMRGridCalendarTime.erase();
   _measurePlotAugmentedAMRGridCalendarTime.erase();
   _measureSolutionAdjustmentAndGlobalTimeStepComputationCalendarTime.erase();
   _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.erase();
   _measurePredictorAndGlobalTimeStepComputationCalendarTime.erase();
   _measureADERDGTimeStepCalendarTime.erase();
   _measureADERDGTimeStepAndPlotCalendarTime.erase();
   _measurePredictorRerunCalendarTime.erase();
   _measureRiemannSolverCalendarTime.erase();
   _measurePredictorCalendarTime.erase();
   _measureCorrectorCalendarTime.erase();
   _measureCorrectorAndPlotCalendarTime.erase();
   _measurePlotCalendarTime.erase();

}
