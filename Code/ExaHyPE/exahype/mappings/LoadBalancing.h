// This file originally was created by pdt (Peano Development Toolkit) as part
// of a code based upon the Peano project by Tobias Weinzierl. For conditions 
// of distribution and use of this project, please see the copyright notice at
// www.peano-framework.org. Feel free to adopt the license and authorship of 
// this file and your project to your needs as long as the license is in 
// agreement with the original Peano user constraints. A reference to/citation  
// of  Peano and its author is highly appreciated.
#ifndef EXAHYPE_MAPPINGS_LoadBalancing_H_
#define EXAHYPE_MAPPINGS_LoadBalancing_H_


#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"

#include "peano/grid/VertexEnumerator.h"
#include "peano/MappingSpecification.h"
#include "peano/CommunicationSpecification.h"

#include "tarch/multicore/MulticoreDefinitions.h"

#include "exahype/Vertex.h"
#include "exahype/Cell.h"
#include "exahype/State.h"


namespace exahype {
namespace mappings {
class LoadBalancing;
}
}


/**
 * Compute the load balancing metrics bottom-up
 *
 * The mapping plugs into leaveCell() only and basically realises the ideas
 * from the class documentation of mpibalancing::HotSpotBalancing.
 * 
 * @author Tobias Weinzierl
 * @version $Revision: 1.10 $
 */
class exahype::mappings::LoadBalancing {
private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log  _log;
public:
  /**
   * Nothing to be done
   */
  static peano::MappingSpecification   touchVertexLastTimeSpecification();
  /**
   * Nothing to be done
   */
  static peano::MappingSpecification   touchVertexFirstTimeSpecification();

  /**
   * Operation degenerates to nop if we translate without MPI.
   *
   * @see enterCell() where we clear the workload flags
   */
  static peano::MappingSpecification   enterCellSpecification();

  /**
   * Operation degenerates to nop if we translate without MPI.
   *
   * @see leaveCell() for a discussion of the restriction mechanism.
   */
  static peano::MappingSpecification   leaveCellSpecification();

  /**
   * Nothing to be done
   */
  static peano::MappingSpecification   ascendSpecification();

  /**
   * Nothing to be done
   */
  static peano::MappingSpecification   descendSpecification();

  /**
   * The load balancing does rely on an analysed tree grammar, i.e.
   * information propagates from fine grid levels to coarser levels.
   * However, we do not enforce any worker-master data exchange in this
   * mapping. If another mapping does require worker-master data exchange
   * (as a global time step size or time stamp has to be computed for
   * example), we do plug into the event and also update the load balancing.
   * But in general, we do not enforce and worker-master synchronisation at
   * all.
   */
  static peano::CommunicationSpecification   communicationSpecification();


  /**
   * Nop.
   */
  LoadBalancing();

#if defined(SharedMemoryParallelisation)
  /**
   * Nop.
   */
  LoadBalancing(const LoadBalancing& masterThread);
#endif

  /**
   * Nop.
   */
  virtual ~LoadBalancing();

#if defined(SharedMemoryParallelisation)
  /**
   * Nop.
   */
  void mergeWithWorkerThread(const LoadBalancing& workerThread);
#endif

  /**
   * Nop.
   */
  void createInnerVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
  );


  /**
   * Nop.
   */
  void createBoundaryVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
  );


  /**
   * Nop.
   */
  void createHangingVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
  );


  /**
   * Nop.
   */
  void destroyHangingVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
  );


  /**
   * Nop.
   */
  void destroyVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
  );


  /**
   * Nop.
   */
  void createCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const         fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
  );


  /**
   * Nop.
   */
  void destroyCell(
      const exahype::Cell&           fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
  );


#ifdef Parallel
  /**
   * Nop.
   */
  void mergeWithNeighbour(
      exahype::Vertex&  vertex,
      const exahype::Vertex&  neighbour,
      int                                           fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
  );


  /**
   * Nop.
   */
  void prepareSendToNeighbour(
      exahype::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
  );


  /**
   * Nop.
   */
  void prepareCopyToRemoteNode(
      exahype::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
  );


  /**
   * Nop.
   */
  void prepareCopyToRemoteNode(
      exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
  );


  /**
   * Nop.
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Vertex&  localVertex,
      const exahype::Vertex&  masterOrWorkerVertex,
      int                                          fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&  x,
      const tarch::la::Vector<DIMENSIONS,double>&  h,
      int                                          level
  );


  /**
   * Nop.
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Cell&  localCell,
      const exahype::Cell&  masterOrWorkerCell,
      int                                          fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
  );


  /**
   * Prepare startup send to worker
   *
   * This operation is called always when we send data to a worker. It is not
   * called when we are right in a join or fork. The operation is kind of the
   * replacement of enterCell() on the master, i.e. called for this one instead
   * of.
   *
   * !!! Reduction
   *
   * With the result, you can control whether the worker shall send back its
   * data from the master rank or not. In accordance, prepareSendToMaster() and
   * corresponding receive are then not called if reduction is switched off.
   * However, the result is only a recommendation: On the one hand, the
   * results of all the active mappings is combined. If one of them requires
   * reduction, Peano does reduce data to the master. On the other hand, Peano
   * itself might decide that it reduces nevertheless. The latter case
   * happens if the master decides that load balancing should be made.
   *
   *
   * @see peano::kernel::spacetreegrid::nodes::Node::updateCellsParallelStateAfterLoad()
   * @return Whether this node needs to send back data to its master.
   */
  bool prepareSendToWorker(
      exahype::Cell&                       fineGridCell,
      exahype::Vertex * const              fineGridVertices,
      const peano::grid::VertexEnumerator&       fineGridVerticesEnumerator,
      exahype::Vertex * const              coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      exahype::Cell&                       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell,
      int                                        worker
  );


  /**
   * TODO(Tobias): Add docu.
   */
  void mergeWithMaster(
      const exahype::Cell&                       workerGridCell,
      exahype::Vertex * const                    workerGridVertices,
      const peano::grid::VertexEnumerator&             workerEnumerator,
      exahype::Cell&                             fineGridCell,
      exahype::Vertex * const                    fineGridVertices,
      const peano::grid::VertexEnumerator&             fineGridVerticesEnumerator,
      exahype::Vertex * const                    coarseGridVertices,
      const peano::grid::VertexEnumerator&             coarseGridVerticesEnumerator,
      exahype::Cell&                             coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&         fineGridPositionOfCell,
      int                                              worker,
      const exahype::State&                      workerState,
      exahype::State&                            masterState
  );


  /**
   * Nop.
   */
  void prepareSendToMaster(
      exahype::Cell&                       localCell,
      exahype::Vertex *                    vertices,
      const peano::grid::VertexEnumerator&       verticesEnumerator, 
      const exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      const exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
  );


  /**
   * Nop.
   */
  void receiveDataFromMaster(
      exahype::Cell&                        receivedCell, 
      exahype::Vertex *                     receivedVertices,
      const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
      exahype::Vertex * const               receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
      exahype::Cell&                        receivedCoarseGridCell,
      exahype::Vertex * const               workersCoarseGridVertices,
      const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
      exahype::Cell&                        workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
  );


  /**
   * TODO(Tobias): Add docu.
   */
  void mergeWithWorker(
      exahype::Cell&           localCell, 
      const exahype::Cell&     receivedMasterCell,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
  );


  /**
   * Nop.
   */
  void mergeWithWorker(
      exahype::Vertex&        localVertex,
      const exahype::Vertex&  receivedMasterVertex,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
  );
#endif


  /**
   * Nop.
   */
  void touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
  );



  /**
   * Nop.
   */
  void touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
  );


  /**
   * TODO(Tobias): Add docu.
   */
  void enterCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
  );


  /**
   * TODO(Tobias): Add docu.
   */
  void leaveCell(
      exahype::Cell&                          fineGridCell,
      exahype::Vertex * const                 fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const                 coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&                          coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&      fineGridPositionOfCell
  );


  /**
   * Nop.
   */
  void beginIteration(
      exahype::State&  solverState
  );


  /**
   * Nop.
   */
  void endIteration(
      exahype::State&  solverState
  );


  /**
   * Nop.
   */
  void descend(
      exahype::Cell * const          fineGridCells,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell
  );


  /**
   * Nop.
   */
  void ascend(
      exahype::Cell * const    fineGridCells,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell
  );
};


#endif
