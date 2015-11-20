// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STD_H_ 
#define _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STD_H_ 


#include "ExaHyPE/repositories/Repository.h"
#include "ExaHyPE/records/RepositoryState.h"

#include "ExaHyPE/State.h"
#include "ExaHyPE/Vertex.h"
#include "ExaHyPE/Cell.h"

#include "peano/grid/Grid.h"
#include "peano/stacks/CellSTDStack.h"
#include "peano/stacks/VertexSTDStack.h"


 #include "ExaHyPE/adapters/CreateGrid.h" 



namespace ExaHyPE {
      namespace repositories {
        class RepositorySTDStack;  
      }
}


class ExaHyPE::repositories::RepositorySTDStack: public ExaHyPE::repositories::Repository {
  private:
    static tarch::logging::Log _log;
  
    peano::geometry::Geometry& _geometry;
    
    typedef peano::stacks::CellSTDStack<ExaHyPE::Cell>       CellStack;
    typedef peano::stacks::VertexSTDStack<ExaHyPE::Vertex>   VertexStack;

    CellStack    _cellStack;
    VertexStack  _vertexStack;
    ExaHyPE::State          _solverState;
    peano::grid::RegularGridContainer<ExaHyPE::Vertex,ExaHyPE::Cell>  _regularGridContainer;
    peano::grid::TraversalOrderOnTopLevel                                         _traversalOrderOnTopLevel;

    peano::grid::Grid<ExaHyPE::Vertex,ExaHyPE::Cell,ExaHyPE::State,VertexStack,CellStack,ExaHyPE::adapters::CreateGrid> _gridWithCreateGrid;

     
   ExaHyPE::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureCreateGridCPUTime;

    tarch::timing::Measurement _measureCreateGridCalendarTime;

   
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
        
    virtual ExaHyPE::State& getState();
    virtual const ExaHyPE::State& getState() const;
	
    virtual void iterate(int numberOfIterations=1, bool exchangeBoundaryVertices=true);

    virtual void writeCheckpoint(peano::grid::Checkpoint<ExaHyPE::Vertex, ExaHyPE::Cell> * const checkpoint); 
    virtual void readCheckpoint( peano::grid::Checkpoint<ExaHyPE::Vertex, ExaHyPE::Cell> const * const checkpoint );
    virtual peano::grid::Checkpoint<ExaHyPE::Vertex, ExaHyPE::Cell>* createEmptyCheckpoint(); 

    virtual void switchToCreateGrid();    

    virtual bool isActiveAdapterCreateGrid() const;

   
    #ifdef Parallel
    virtual ContinueCommand continueToIterate();
    virtual void runGlobalStep();
    #endif

    virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value);
    virtual void logIterationStatistics() const;
    virtual void clearIterationStatistics();
};


#endif
