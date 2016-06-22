/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
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
 #include "exahype/adapters/SolutionUpdateAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/PredictorAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/ADERDGTimeStep.h" 
 #include "exahype/adapters/ADERDGTimeStepAndPlot.h" 
 #include "exahype/adapters/GlobalTimeStepComputation.h" 
 #include "exahype/adapters/GlobalTimeStepComputationAndPlot.h" 
 #include "exahype/adapters/FaceDataExchange.h" 
 #include "exahype/adapters/Predictor.h" 
 #include "exahype/adapters/PredictorRerun.h" 
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

    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::AugmentedAMRGrid> _gridWithAugmentedAMRGrid;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PlotAugmentedAMRGrid> _gridWithPlotAugmentedAMRGrid;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::SolutionUpdateAndGlobalTimeStepComputation> _gridWithSolutionUpdateAndGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictorAndGlobalTimeStepComputation> _gridWithPredictorAndGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::ADERDGTimeStep> _gridWithADERDGTimeStep;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::ADERDGTimeStepAndPlot> _gridWithADERDGTimeStepAndPlot;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::GlobalTimeStepComputation> _gridWithGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::GlobalTimeStepComputationAndPlot> _gridWithGlobalTimeStepComputationAndPlot;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::FaceDataExchange> _gridWithFaceDataExchange;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Predictor> _gridWithPredictor;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictorRerun> _gridWithPredictorRerun;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Corrector> _gridWithCorrector;

     
   exahype::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureAugmentedAMRGridCPUTime;
    tarch::timing::Measurement _measurePlotAugmentedAMRGridCPUTime;
    tarch::timing::Measurement _measureSolutionUpdateAndGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measurePredictorAndGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measureADERDGTimeStepCPUTime;
    tarch::timing::Measurement _measureADERDGTimeStepAndPlotCPUTime;
    tarch::timing::Measurement _measureGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measureGlobalTimeStepComputationAndPlotCPUTime;
    tarch::timing::Measurement _measureFaceDataExchangeCPUTime;
    tarch::timing::Measurement _measurePredictorCPUTime;
    tarch::timing::Measurement _measurePredictorRerunCPUTime;
    tarch::timing::Measurement _measureCorrectorCPUTime;

    tarch::timing::Measurement _measureAugmentedAMRGridCalendarTime;
    tarch::timing::Measurement _measurePlotAugmentedAMRGridCalendarTime;
    tarch::timing::Measurement _measureSolutionUpdateAndGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measurePredictorAndGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measureADERDGTimeStepCalendarTime;
    tarch::timing::Measurement _measureADERDGTimeStepAndPlotCalendarTime;
    tarch::timing::Measurement _measureGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measureGlobalTimeStepComputationAndPlotCalendarTime;
    tarch::timing::Measurement _measureFaceDataExchangeCalendarTime;
    tarch::timing::Measurement _measurePredictorCalendarTime;
    tarch::timing::Measurement _measurePredictorRerunCalendarTime;
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

    virtual void switchToAugmentedAMRGrid();    
    virtual void switchToPlotAugmentedAMRGrid();    
    virtual void switchToSolutionUpdateAndGlobalTimeStepComputation();    
    virtual void switchToPredictorAndGlobalTimeStepComputation();    
    virtual void switchToADERDGTimeStep();    
    virtual void switchToADERDGTimeStepAndPlot();    
    virtual void switchToGlobalTimeStepComputation();    
    virtual void switchToGlobalTimeStepComputationAndPlot();    
    virtual void switchToFaceDataExchange();    
    virtual void switchToPredictor();    
    virtual void switchToPredictorRerun();    
    virtual void switchToCorrector();    

    virtual bool isActiveAdapterAugmentedAMRGrid() const;
    virtual bool isActiveAdapterPlotAugmentedAMRGrid() const;
    virtual bool isActiveAdapterSolutionUpdateAndGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterPredictorAndGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterADERDGTimeStep() const;
    virtual bool isActiveAdapterADERDGTimeStepAndPlot() const;
    virtual bool isActiveAdapterGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterGlobalTimeStepComputationAndPlot() const;
    virtual bool isActiveAdapterFaceDataExchange() const;
    virtual bool isActiveAdapterPredictor() const;
    virtual bool isActiveAdapterPredictorRerun() const;
    virtual bool isActiveAdapterCorrector() const;

   
    #ifdef Parallel
    virtual ContinueCommand continueToIterate();
    virtual void runGlobalStep();
    #endif

    virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value);
    virtual void logIterationStatistics(bool logAllAdapters) const;
    virtual void clearIterationStatistics();
};


#endif
