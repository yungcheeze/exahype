// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef MYPROJECT_ADAPTERS_TimeStepAndPlot_H_
#define MYPROJECT_ADAPTERS_TimeStepAndPlot_H_


#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"

#include "peano/grid/VertexEnumerator.h"
#include "peano/MappingSpecification.h"
#include "peano/CommunicationSpecification.h"

#include "tarch/multicore/MulticoreDefinitions.h"

#include "ExplicitEulerForHeatEquation/Vertex.h"
#include "ExplicitEulerForHeatEquation/Cell.h"
#include "ExplicitEulerForHeatEquation/State.h"

 #include "ExplicitEulerForHeatEquation/mappings/TimeStep.h"
 #include "ExplicitEulerForHeatEquation/mappings/InterpolateHangingNodes.h"
 #include "ExplicitEulerForHeatEquation/mappings/Inject.h"

 #include "ExplicitEulerForHeatEquation/adapters/TimeStepAndPlot2VTKPlotVertexValue_0.h"



namespace myproject {
      namespace adapters {
        class TimeStepAndPlot;
      } 
}


/**
 * This is a mapping from the spacetree traversal events to your user-defined activities.
 * The latter are realised within the mappings. 
 * 
 * @author Peano Development Toolkit (PDT) by  Tobias Weinzierl
 * @version $Revision: 1.10 $
 */
class myproject::adapters::TimeStepAndPlot {
  private:
    myproject::mappings::TimeStep _map2TimeStep;
    myproject::mappings::InterpolateHangingNodes _map2InterpolateHangingNodes;
    myproject::mappings::Inject _map2Inject;

    myproject::adapters::TimeStepAndPlot2VTKPlotVertexValue_0 _map2TimeStepAndPlot2VTKPlotVertexValue_0;

  public:
    static peano::MappingSpecification         touchVertexLastTimeSpecification();
    static peano::MappingSpecification         touchVertexFirstTimeSpecification();
    static peano::MappingSpecification         enterCellSpecification();
    static peano::MappingSpecification         leaveCellSpecification();
    static peano::MappingSpecification         ascendSpecification();
    static peano::MappingSpecification         descendSpecification();
    static peano::CommunicationSpecification   communicationSpecification();

    TimeStepAndPlot();

    #if defined(SharedMemoryParallelisation)
    TimeStepAndPlot(const TimeStepAndPlot& masterThread);
    #endif

    virtual ~TimeStepAndPlot();
  
    #if defined(SharedMemoryParallelisation)
    void mergeWithWorkerThread(const TimeStepAndPlot& workerThread);
    #endif

    void createInnerVertex(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void createBoundaryVertex(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void createHangingVertex(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void destroyHangingVertex(
      const myproject::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );


    void destroyVertex(
      const myproject::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );


    void createCell(
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const         fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );


    void destroyCell(
      const myproject::Cell&           fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );
        
    #ifdef Parallel
    void mergeWithNeighbour(
      myproject::Vertex&  vertex,
      const myproject::Vertex&  neighbour,
      int                                           fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareSendToNeighbour(
      myproject::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareCopyToRemoteNode(
      myproject::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareCopyToRemoteNode(
      myproject::Cell&  localCell,
      int  toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
    );

    void mergeWithRemoteDataDueToForkOrJoin(
      myproject::Vertex&  localVertex,
      const myproject::Vertex&  masterOrWorkerVertex,
      int                                       fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&  x,
      const tarch::la::Vector<DIMENSIONS,double>&  h,
      int                                       level
    );

    void mergeWithRemoteDataDueToForkOrJoin(
      myproject::Cell&  localCell,
      const myproject::Cell&  masterOrWorkerCell,
      int                                       fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                       level
    );

    bool prepareSendToWorker(
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
      int                                                                  worker
    );

    void prepareSendToMaster(
      myproject::Cell&                       localCell,
      myproject::Vertex *                    vertices,
      const peano::grid::VertexEnumerator&       verticesEnumerator, 
      const myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      const myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
    );

    void mergeWithMaster(
      const myproject::Cell&           workerGridCell,
      myproject::Vertex * const        workerGridVertices,
      const peano::grid::VertexEnumerator& workerEnumerator,
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
      int                                                                  worker,
      const myproject::State&           workerState,
      myproject::State&                 masterState
    );


    void receiveDataFromMaster(
      myproject::Cell&                        receivedCell, 
      myproject::Vertex *                     receivedVertices,
      const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
      myproject::Vertex * const               receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
      myproject::Cell&                        receivedCoarseGridCell,
      myproject::Vertex * const               workersCoarseGridVertices,
      const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
      myproject::Cell&                        workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
    );


    void mergeWithWorker(
      myproject::Cell&           localCell, 
      const myproject::Cell&     receivedMasterCell,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
    );


    void mergeWithWorker(
      myproject::Vertex&        localVertex,
      const myproject::Vertex&  receivedMasterVertex,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );
    #endif


    void touchVertexFirstTime(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void touchVertexLastTime(
      myproject::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );
    

    void enterCell(
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );


    void leaveCell(
      myproject::Cell&                          fineGridCell,
      myproject::Vertex * const                 fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      myproject::Vertex * const                 coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&                          coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&      fineGridPositionOfCell
    );


    void beginIteration(
      myproject::State&  solverState
    );


    void endIteration(
      myproject::State&  solverState
    );

    void descend(
      myproject::Cell * const          fineGridCells,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell
    );


    void ascend(
      myproject::Cell * const    fineGridCells,
      myproject::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell
    );    
};


#endif
