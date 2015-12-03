// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STACK_H_ 
#define _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STACK_H_ 


#include "EulerFlow3d/repositories/Repository.h"
#include "EulerFlow3d/records/RepositoryState.h"

#include "EulerFlow3d/State.h"
#include "EulerFlow3d/Vertex.h"
#include "EulerFlow3d/Cell.h"

#include "peano/grid/Grid.h"
#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/VertexArrayStack.h"


 #include "EulerFlow3d/adapters/InitialGrid.h" 
 #include "EulerFlow3d/adapters/GridExport.h" 
 #include "EulerFlow3d/adapters/PatchInit.h" 
 #include "EulerFlow3d/adapters/InitialCondition.h" 
 #include "EulerFlow3d/adapters/Predictor.h" 
 #include "EulerFlow3d/adapters/Corrector.h" 
 #include "EulerFlow3d/adapters/SolutionExport.h" 



namespace exahype {
      namespace repositories {
        class RepositoryArrayStack;  
      }
}


class exahype::repositories::RepositoryArrayStack: public exahype::repositories::Repository {
  private:
    static tarch::logging::Log _log;
  
    peano::geometry::Geometry& _geometry;
    
    typedef peano::stacks::CellArrayStack<exahype::Cell>       CellStack;
    typedef peano::stacks::VertexArrayStack<exahype::Vertex>   VertexStack;

    CellStack    _cellStack;
    VertexStack  _vertexStack;
    exahype::State          _solverState;
    peano::grid::RegularGridContainer<exahype::Vertex,exahype::Cell>  _regularGridContainer;
    peano::grid::TraversalOrderOnTopLevel                                         _traversalOrderOnTopLevel;

    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::InitialGrid> _gridWithInitialGrid;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::GridExport> _gridWithGridExport;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PatchInit> _gridWithPatchInit;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::InitialCondition> _gridWithInitialCondition;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Predictor> _gridWithPredictor;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Corrector> _gridWithCorrector;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::SolutionExport> _gridWithSolutionExport;

  
   exahype::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureInitialGridCPUTime;
    tarch::timing::Measurement _measureGridExportCPUTime;
    tarch::timing::Measurement _measurePatchInitCPUTime;
    tarch::timing::Measurement _measureInitialConditionCPUTime;
    tarch::timing::Measurement _measurePredictorCPUTime;
    tarch::timing::Measurement _measureCorrectorCPUTime;
    tarch::timing::Measurement _measureSolutionExportCPUTime;

    tarch::timing::Measurement _measureInitialGridCalendarTime;
    tarch::timing::Measurement _measureGridExportCalendarTime;
    tarch::timing::Measurement _measurePatchInitCalendarTime;
    tarch::timing::Measurement _measureInitialConditionCalendarTime;
    tarch::timing::Measurement _measurePredictorCalendarTime;
    tarch::timing::Measurement _measureCorrectorCalendarTime;
    tarch::timing::Measurement _measureSolutionExportCalendarTime;


  public:
    RepositoryArrayStack(
      peano::geometry::Geometry&                   geometry,
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset,
      int                                          maximumSizeOfCellInOutStack,
      int                                          maximumSizeOfVertexInOutStack,
      int                                          maximumSizeOfVertexTemporaryStack
    );
    
    /**
     * Parallel Constructor
     *
     * Used in parallel mode only where the size of the domain is not known 
     * when the type of repository is determined.  
     */
    RepositoryArrayStack(
      peano::geometry::Geometry&                   geometry,
      int                                          maximumSizeOfCellInOutStack,
      int                                          maximumSizeOfVertexInOutStack,
      int                                          maximumSizeOfVertexTemporaryStack
    );
    
    virtual ~RepositoryArrayStack();

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
    virtual void switchToGridExport();    
    virtual void switchToPatchInit();    
    virtual void switchToInitialCondition();    
    virtual void switchToPredictor();    
    virtual void switchToCorrector();    
    virtual void switchToSolutionExport();    

    virtual bool isActiveAdapterInitialGrid() const;
    virtual bool isActiveAdapterGridExport() const;
    virtual bool isActiveAdapterPatchInit() const;
    virtual bool isActiveAdapterInitialCondition() const;
    virtual bool isActiveAdapterPredictor() const;
    virtual bool isActiveAdapterCorrector() const;
    virtual bool isActiveAdapterSolutionExport() const;

     
    #ifdef Parallel
    virtual ContinueCommand continueToIterate();
    virtual void runGlobalStep();
    #endif

    virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value);
    virtual void logIterationStatistics() const;
    virtual void clearIterationStatistics();
};


#endif
