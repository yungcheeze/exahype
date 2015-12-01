// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STD_H_ 
#define _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STD_H_ 


#include "EulerFlow3d/repositories/Repository.h"
#include "EulerFlow3d/records/RepositoryState.h"

#include "EulerFlow3d/State.h"
#include "EulerFlow3d/Vertex.h"
#include "EulerFlow3d/Cell.h"

#include "peano/grid/Grid.h"
#include "peano/stacks/CellSTDStack.h"
#include "peano/stacks/VertexSTDStack.h"


 #include "EulerFlow3d/adapters/CreateGrid.h" 
 #include "EulerFlow3d/adapters/PlotGrid.h" 
 #include "EulerFlow3d/adapters/InitCells.h" 
 #include "EulerFlow3d/adapters/InitCellData.h" 
 #include "EulerFlow3d/adapters/TimeStep.h" 
 #include "EulerFlow3d/adapters/TimeStepAndPlot.h" 



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

    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::CreateGrid> _gridWithCreateGrid;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PlotGrid> _gridWithPlotGrid;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::InitCells> _gridWithInitCells;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::InitCellData> _gridWithInitCellData;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::TimeStep> _gridWithTimeStep;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::TimeStepAndPlot> _gridWithTimeStepAndPlot;

     
   exahype::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureCreateGridCPUTime;
    tarch::timing::Measurement _measurePlotGridCPUTime;
    tarch::timing::Measurement _measureInitCellsCPUTime;
    tarch::timing::Measurement _measureInitCellDataCPUTime;
    tarch::timing::Measurement _measureTimeStepCPUTime;
    tarch::timing::Measurement _measureTimeStepAndPlotCPUTime;

    tarch::timing::Measurement _measureCreateGridCalendarTime;
    tarch::timing::Measurement _measurePlotGridCalendarTime;
    tarch::timing::Measurement _measureInitCellsCalendarTime;
    tarch::timing::Measurement _measureInitCellDataCalendarTime;
    tarch::timing::Measurement _measureTimeStepCalendarTime;
    tarch::timing::Measurement _measureTimeStepAndPlotCalendarTime;

   
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

    virtual void switchToCreateGrid();    
    virtual void switchToPlotGrid();    
    virtual void switchToInitCells();    
    virtual void switchToInitCellData();    
    virtual void switchToTimeStep();    
    virtual void switchToTimeStepAndPlot();    

    virtual bool isActiveAdapterCreateGrid() const;
    virtual bool isActiveAdapterPlotGrid() const;
    virtual bool isActiveAdapterInitCells() const;
    virtual bool isActiveAdapterInitCellData() const;
    virtual bool isActiveAdapterTimeStep() const;
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
