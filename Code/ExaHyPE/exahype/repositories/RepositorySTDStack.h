// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STD_H_ 
#define _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STD_H_ 


#include "exahype/repositories/Repository.h"
#include "exahype/records/RepositoryState.h"

#include "exahype/State.h"
#include "exahype/Vertex.h"
#include "exahype/Cell.h"

#include "peano/grid/Grid.h"
#include "peano/stacks/CellSTDStack.h"
#include "peano/stacks/VertexSTDStack.h"


 #include "exahype/adapters/AugmentedAMRGrid.h" 
 #include "exahype/adapters/PlotAugmentedAMRGrid.h" 
 #include "exahype/adapters/InitialConditionAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/PredictorAndPlotAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/PredictorAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/GridErasing.h" 
 #include "exahype/adapters/ADERDGTimeStep.h" 
 #include "exahype/adapters/ADERDGTimeStepAndPlot.h" 
 #include "exahype/adapters/PredictorRerun.h" 
 #include "exahype/adapters/RiemannSolver.h" 
 #include "exahype/adapters/Predictor.h" 
 #include "exahype/adapters/Corrector.h" 
 #include "exahype/adapters/CorrectorAndPlot.h" 
 #include "exahype/adapters/Plot.h" 



namespace exahype {
      namespace repositories {
        class RepositorySTDStack;  
      }
}


class exahype::repositories::RepositorySTDStack: public exahype::repositories::Repository {
  private:
    static tarch::logging::Log _log;
  
    peano::geometry::Geometry& _geometry;
    
    typedef peano::stacks::CellSTDStack<exahype::Cell>       CellStack;
    typedef peano::stacks::VertexSTDStack<exahype::Vertex>   VertexStack;

    CellStack    _cellStack;
    VertexStack  _vertexStack;
    exahype::State          _solverState;
    peano::grid::RegularGridContainer<exahype::Vertex,exahype::Cell>  _regularGridContainer;
    peano::grid::TraversalOrderOnTopLevel                                         _traversalOrderOnTopLevel;

    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::AugmentedAMRGrid> _gridWithAugmentedAMRGrid;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PlotAugmentedAMRGrid> _gridWithPlotAugmentedAMRGrid;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::InitialConditionAndGlobalTimeStepComputation> _gridWithInitialConditionAndGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictorAndPlotAndGlobalTimeStepComputation> _gridWithPredictorAndPlotAndGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictorAndGlobalTimeStepComputation> _gridWithPredictorAndGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::GridErasing> _gridWithGridErasing;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::ADERDGTimeStep> _gridWithADERDGTimeStep;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::ADERDGTimeStepAndPlot> _gridWithADERDGTimeStepAndPlot;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictorRerun> _gridWithPredictorRerun;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::RiemannSolver> _gridWithRiemannSolver;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Predictor> _gridWithPredictor;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Corrector> _gridWithCorrector;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::CorrectorAndPlot> _gridWithCorrectorAndPlot;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Plot> _gridWithPlot;

     
   exahype::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureAugmentedAMRGridCPUTime;
    tarch::timing::Measurement _measurePlotAugmentedAMRGridCPUTime;
    tarch::timing::Measurement _measureInitialConditionAndGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measurePredictorAndGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measureGridErasingCPUTime;
    tarch::timing::Measurement _measureADERDGTimeStepCPUTime;
    tarch::timing::Measurement _measureADERDGTimeStepAndPlotCPUTime;
    tarch::timing::Measurement _measurePredictorRerunCPUTime;
    tarch::timing::Measurement _measureRiemannSolverCPUTime;
    tarch::timing::Measurement _measurePredictorCPUTime;
    tarch::timing::Measurement _measureCorrectorCPUTime;
    tarch::timing::Measurement _measureCorrectorAndPlotCPUTime;
    tarch::timing::Measurement _measurePlotCPUTime;

    tarch::timing::Measurement _measureAugmentedAMRGridCalendarTime;
    tarch::timing::Measurement _measurePlotAugmentedAMRGridCalendarTime;
    tarch::timing::Measurement _measureInitialConditionAndGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measurePredictorAndGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measureGridErasingCalendarTime;
    tarch::timing::Measurement _measureADERDGTimeStepCalendarTime;
    tarch::timing::Measurement _measureADERDGTimeStepAndPlotCalendarTime;
    tarch::timing::Measurement _measurePredictorRerunCalendarTime;
    tarch::timing::Measurement _measureRiemannSolverCalendarTime;
    tarch::timing::Measurement _measurePredictorCalendarTime;
    tarch::timing::Measurement _measureCorrectorCalendarTime;
    tarch::timing::Measurement _measureCorrectorAndPlotCalendarTime;
    tarch::timing::Measurement _measurePlotCalendarTime;

   
  public:
    RepositorySTDStack(
      peano::geometry::Geometry&                   geometry,
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset
    );
    
    /**
     * Parallel Constructor
     *
     * Used in parallel mode only where the size of the domain is not known 
     * when the type of repository is determined.  
     */
    RepositorySTDStack(
      peano::geometry::Geometry&                   geometry
    );
    
    virtual ~RepositorySTDStack();

    virtual void restart(
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
      int                                          domainLevel,
      const tarch::la::Vector<DIMENSIONS,int>&     positionOfCentralElementWithRespectToCoarserRemoteLevel
    );
         
    virtual void terminate();
        
    virtual exahype::State& getState();
    virtual const exahype::State& getState() const;
	
    virtual void iterate(int numberOfIterations=1, bool exchangeBoundaryVertices=true);

    virtual void writeCheckpoint(peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> * const checkpoint); 
    virtual void readCheckpoint( peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> const * const checkpoint );
    virtual peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>* createEmptyCheckpoint(); 

    virtual void switchToAugmentedAMRGrid();    
    virtual void switchToPlotAugmentedAMRGrid();    
    virtual void switchToInitialConditionAndGlobalTimeStepComputation();    
    virtual void switchToPredictorAndPlotAndGlobalTimeStepComputation();    
    virtual void switchToPredictorAndGlobalTimeStepComputation();    
    virtual void switchToGridErasing();    
    virtual void switchToADERDGTimeStep();    
    virtual void switchToADERDGTimeStepAndPlot();    
    virtual void switchToPredictorRerun();    
    virtual void switchToRiemannSolver();    
    virtual void switchToPredictor();    
    virtual void switchToCorrector();    
    virtual void switchToCorrectorAndPlot();    
    virtual void switchToPlot();    

    virtual bool isActiveAdapterAugmentedAMRGrid() const;
    virtual bool isActiveAdapterPlotAugmentedAMRGrid() const;
    virtual bool isActiveAdapterInitialConditionAndGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterPredictorAndPlotAndGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterPredictorAndGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterGridErasing() const;
    virtual bool isActiveAdapterADERDGTimeStep() const;
    virtual bool isActiveAdapterADERDGTimeStepAndPlot() const;
    virtual bool isActiveAdapterPredictorRerun() const;
    virtual bool isActiveAdapterRiemannSolver() const;
    virtual bool isActiveAdapterPredictor() const;
    virtual bool isActiveAdapterCorrector() const;
    virtual bool isActiveAdapterCorrectorAndPlot() const;
    virtual bool isActiveAdapterPlot() const;

   
    #ifdef Parallel
    virtual ContinueCommand continueToIterate();
    virtual void runGlobalStep();
    #endif

    virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value);
    virtual void logIterationStatistics(bool logAllAdapters) const;
    virtual void clearIterationStatistics();
};


#endif
