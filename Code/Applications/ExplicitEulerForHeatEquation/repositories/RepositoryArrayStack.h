// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _MYPROJECT_REPOSITORIES_REPOSITORY_ARRAY_STACK_H_ 
#define _MYPROJECT_REPOSITORIES_REPOSITORY_ARRAY_STACK_H_ 


#include "ExplicitEulerForHeatEquation/repositories/Repository.h"
#include "ExplicitEulerForHeatEquation/records/RepositoryState.h"

#include "ExplicitEulerForHeatEquation/State.h"
#include "ExplicitEulerForHeatEquation/Vertex.h"
#include "ExplicitEulerForHeatEquation/Cell.h"

#include "peano/grid/Grid.h"
#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/VertexArrayStack.h"


 #include "ExplicitEulerForHeatEquation/adapters/CreateGrid.h" 
 #include "ExplicitEulerForHeatEquation/adapters/TimeStep.h" 
 #include "ExplicitEulerForHeatEquation/adapters/CreateGridAndPlot.h" 
 #include "ExplicitEulerForHeatEquation/adapters/TimeStepAndPlot.h" 



namespace myproject {
      namespace repositories {
        class RepositoryArrayStack;  
      }
}


class myproject::repositories::RepositoryArrayStack: public myproject::repositories::Repository {
  private:
    static tarch::logging::Log _log;
  
    peano::geometry::Geometry& _geometry;
    
    typedef peano::stacks::CellArrayStack<myproject::Cell>       CellStack;
    typedef peano::stacks::VertexArrayStack<myproject::Vertex>   VertexStack;

    CellStack    _cellStack;
    VertexStack  _vertexStack;
    myproject::State          _solverState;
    peano::grid::RegularGridContainer<myproject::Vertex,myproject::Cell>  _regularGridContainer;
    peano::grid::TraversalOrderOnTopLevel                                         _traversalOrderOnTopLevel;

    peano::grid::Grid<myproject::Vertex,myproject::Cell,myproject::State,VertexStack,CellStack,myproject::adapters::CreateGrid> _gridWithCreateGrid;
    peano::grid::Grid<myproject::Vertex,myproject::Cell,myproject::State,VertexStack,CellStack,myproject::adapters::TimeStep> _gridWithTimeStep;
    peano::grid::Grid<myproject::Vertex,myproject::Cell,myproject::State,VertexStack,CellStack,myproject::adapters::CreateGridAndPlot> _gridWithCreateGridAndPlot;
    peano::grid::Grid<myproject::Vertex,myproject::Cell,myproject::State,VertexStack,CellStack,myproject::adapters::TimeStepAndPlot> _gridWithTimeStepAndPlot;

  
   myproject::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureCreateGridCPUTime;
    tarch::timing::Measurement _measureTimeStepCPUTime;
    tarch::timing::Measurement _measureCreateGridAndPlotCPUTime;
    tarch::timing::Measurement _measureTimeStepAndPlotCPUTime;

    tarch::timing::Measurement _measureCreateGridCalendarTime;
    tarch::timing::Measurement _measureTimeStepCalendarTime;
    tarch::timing::Measurement _measureCreateGridAndPlotCalendarTime;
    tarch::timing::Measurement _measureTimeStepAndPlotCalendarTime;


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
        
    virtual myproject::State& getState();
    virtual const myproject::State& getState() const;

    virtual void iterate(int numberOfIterations=1, bool exchangeBoundaryVertices=true);
    
    virtual void writeCheckpoint(peano::grid::Checkpoint<myproject::Vertex, myproject::Cell> * const checkpoint); 
    virtual void readCheckpoint( peano::grid::Checkpoint<myproject::Vertex, myproject::Cell> const * const checkpoint );
    virtual peano::grid::Checkpoint<myproject::Vertex, myproject::Cell>* createEmptyCheckpoint(); 

    virtual void switchToCreateGrid();    
    virtual void switchToTimeStep();    
    virtual void switchToCreateGridAndPlot();    
    virtual void switchToTimeStepAndPlot();    

    virtual bool isActiveAdapterCreateGrid() const;
    virtual bool isActiveAdapterTimeStep() const;
    virtual bool isActiveAdapterCreateGridAndPlot() const;
    virtual bool isActiveAdapterTimeStepAndPlot() const;

     
    #ifdef Parallel
    virtual ContinueCommand continueToIterate();
    virtual void runGlobalStep();
    #endif

    virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value);
    virtual void logIterationStatistics() const;
    virtual void clearIterationStatistics();
};


#endif
