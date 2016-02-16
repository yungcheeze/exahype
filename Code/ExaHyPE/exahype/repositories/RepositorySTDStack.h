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


 #include "exahype/adapters/InitialGrid.h" 
 #include "exahype/adapters/PatchInitialisation.h" 
 #include "exahype/adapters/FaceDataExchange.h" 
 #include "exahype/adapters/InitialConditionAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/PredictorAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/CorrectorAndPredictorAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/Plot.h" 
 #include "exahype/adapters/GlobalTimeStepComputation.h" 
 #include "exahype/adapters/Predictor.h" 
 #include "exahype/adapters/Corrector.h" 



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

    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::InitialGrid> _gridWithInitialGrid;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PatchInitialisation> _gridWithPatchInitialisation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::FaceDataExchange> _gridWithFaceDataExchange;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::InitialConditionAndGlobalTimeStepComputation> _gridWithInitialConditionAndGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictorAndGlobalTimeStepComputation> _gridWithPredictorAndGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::CorrectorAndPredictorAndGlobalTimeStepComputation> _gridWithCorrectorAndPredictorAndGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Plot> _gridWithPlot;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::GlobalTimeStepComputation> _gridWithGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Predictor> _gridWithPredictor;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Corrector> _gridWithCorrector;

     
   exahype::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureInitialGridCPUTime;
    tarch::timing::Measurement _measurePatchInitialisationCPUTime;
    tarch::timing::Measurement _measureFaceDataExchangeCPUTime;
    tarch::timing::Measurement _measureInitialConditionAndGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measurePredictorAndGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measurePlotCPUTime;
    tarch::timing::Measurement _measureGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measurePredictorCPUTime;
    tarch::timing::Measurement _measureCorrectorCPUTime;

    tarch::timing::Measurement _measureInitialGridCalendarTime;
    tarch::timing::Measurement _measurePatchInitialisationCalendarTime;
    tarch::timing::Measurement _measureFaceDataExchangeCalendarTime;
    tarch::timing::Measurement _measureInitialConditionAndGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measurePredictorAndGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measureCorrectorAndPredictorAndGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measurePlotCalendarTime;
    tarch::timing::Measurement _measureGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measurePredictorCalendarTime;
    tarch::timing::Measurement _measureCorrectorCalendarTime;

   
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

    virtual void switchToInitialGrid();    
    virtual void switchToPatchInitialisation();    
    virtual void switchToFaceDataExchange();    
    virtual void switchToInitialConditionAndGlobalTimeStepComputation();    
    virtual void switchToPredictorAndGlobalTimeStepComputation();    
    virtual void switchToCorrectorAndPredictorAndGlobalTimeStepComputation();    
    virtual void switchToPlot();    
    virtual void switchToGlobalTimeStepComputation();    
    virtual void switchToPredictor();    
    virtual void switchToCorrector();    

    virtual bool isActiveAdapterInitialGrid() const;
    virtual bool isActiveAdapterPatchInitialisation() const;
    virtual bool isActiveAdapterFaceDataExchange() const;
    virtual bool isActiveAdapterInitialConditionAndGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterPredictorAndGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterCorrectorAndPredictorAndGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterPlot() const;
    virtual bool isActiveAdapterGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterPredictor() const;
    virtual bool isActiveAdapterCorrector() const;

   
    #ifdef Parallel
    virtual ContinueCommand continueToIterate();
    virtual void runGlobalStep();
    #endif

    virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value);
    virtual void logIterationStatistics() const;
    virtual void clearIterationStatistics();
};


#endif
