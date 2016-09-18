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

#ifndef EXAHYPE_MAPPINGS_Sending_H_
#define EXAHYPE_MAPPINGS_Sending_H_

#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include "peano/CommunicationSpecification.h"
#include "peano/MappingSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/multicore/MulticoreDefinitions.h"

#include "exahype/Cell.h"
#include "exahype/State.h"
#include "exahype/Vertex.h"

namespace exahype {
  namespace mappings {
    class Sending;
  }
}

/**
 * Determine a global time step size
 *
 * The global time step computation runs through all the cells. Per cell, it
 * runs through all involved solvers and determines the corresponding minimal
 * time step sizes. Once the traversal terminates, all solvers thus know what
 * the minimal permitted time step size is. We now take the minimal time step
 * sizes and inform the solvers about them through
 * updateMinNextPredictorTimeStepSize().
 *
 * In the subsequent time step, these minimal time step sizes then are used by
 * synchroniseTimeStepping() (see notably the mapping NewTimeStep) to move the
 * patch forward in time. There is no need to take extra care for the
 * optimistic time step choice - we can determine from outside whether we tend
 * to overshoot and thus have to rerun the predictor. This is done in the
 * runner.
 *
 * <h2>Multicore parallelisation</h2>
 * See documentation of _minTimeStepSizes/
 *
 * <h2>MPI parallelisation</h2>
 *
 * @author Dominic Charrier, Tobias Weinzierl
 */
class exahype::mappings::Sending {
 public:
  #ifdef Parallel
  static bool SkipReductionInBatchedTimeSteps;
  #endif
 private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log _log;

  /**
   * Local copy of the state.
   */
  exahype::State _localState;

  /**
   * A semaphore that is locked if a thread performs a restriction
   * operation.
   */
  static tarch::multicore::BooleanSemaphore _semaphoreForRestriction;

  /**
   * A minimum solver time step size for each thread in a multicore run.
   *
   * We could directly compute the minimal time step sizes and the minimum time
   * stamp: Run per cell through all patch descriptions and update the global
   * solver objects immediately. This is not very clever in a multicore
   * environment as we then have to protect the global data access with a
   * semaphore which serialises the computation (that, in the worst case, is
   * done in each time step). So we do compute the time step size and the stamp
   * locally in a vector, and we project it back to the global data in
   * endIteration().
   */
  std::vector<double> _minTimeStepSizes;

  /**
   * Prepare a appropriately sized vector _minTimeStepSizes
   * with elements initiliased to MAX_DOUBLE.
   */
  void prepareEmptyLocalTimeStepData();

#ifdef Parallel
  /**
   * Tag that is used to exchange all the solver instances in MPI
   */
  static int _mpiTag;

/**
 * Loop over all the solvers and check
 * if a cell description (ADERDGCellDescription,
 * FiniteVolumesCellDescription,...) is registered
 * for the solver type. If so, send
 * out data or empty messages to the rank \p toRank that
 * owns the neighbouring domain.
 *
 * If not so, send out empty messages for the particular
 * solver.
 *
 * TODO(Dominic): Make messaging solver functionality?
 *
 * \note Not thread-safe.
 */
static void sendSolverDataToNeighbour(
    const int                                    toRank,
    const tarch::la::Vector<DIMENSIONS,int>&     src,
    const tarch::la::Vector<DIMENSIONS,int>&     dest,
    const int                                    srcCellDescriptionIndex,
    const int                                    destCellDescriptionIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level);

/**
 * Loop over all the solvers and check
 * send out empty messages for the particular
 * solver.
 *
 * \note Not thread-safe.
 */
static void sendEmptySolverDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS,int>&      src,
    const tarch::la::Vector<DIMENSIONS,int>&      dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level);
#endif

 public:
  /**
   * Run through whole tree. Run concurrently on fine grid.
   */
  static peano::MappingSpecification enterCellSpecification();

  /**
   * Nop.
   */
  static peano::MappingSpecification touchVertexLastTimeSpecification();
  static peano::MappingSpecification touchVertexFirstTimeSpecification();
  static peano::MappingSpecification leaveCellSpecification();
  static peano::MappingSpecification ascendSpecification();
  static peano::MappingSpecification descendSpecification();

  /**
   * The global time step computation does synchronise the individual cells
   * with the solver instances from the master. The actual
   * synchronisation/consistency routines
   * for the solvers are done in NewTimeStep and
   * the SpaceTimePredictor.
   *
   * The fundamental job of the global time step mapping is to report back to
   * the master what time step is permitted. As all those operations are
   * actually done in enterCell---also the veto of a global time step is done
   * in the Riemann solver, i.e. before enterCell---we can send back data as
   * soon as the traversal operation leaves the local subtree.
   */
  static peano::CommunicationSpecification communicationSpecification();

  /**
    * Nop.
    */
  Sending();
  #if defined(SharedMemoryParallelisation)
   /**
    * Nop.
    */
  Sending(const Sending& masterThread);
  #endif

  /**
   * If the fine grid cell functions as compute cell for a solver,
   * compute a stable time step size.
   * Further update the time stamp of the compute cell
   * and update the minimum solver time stamp and
   * time step size.
   */
  void enterCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Restricts face data from \p cellDescription to
   * a parent cell description if the fine grid cell associated with
   * cell description is adjacent to a boundary of the
   * coarse grid cell associated with the parent cell description.
   *
   * \note We use a semaphore to make this operation thread-safe.
   */
  void leaveCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Start iteration/grid sweep.
   * Make the state clear its accumulated values.
   *
   * \note Is called once per rank.
   */
  void beginIteration(exahype::State& solverState);

  /**
   * Runs over all the registered solvers and sets the
   * reduced minimum time step sizes. Then updates the minimum time stamp
   * of the solvers.
   *
   * \note Is called once per rank.
   */
  void endIteration(exahype::State& solverState);
#ifdef Parallel
  /**
   * Prepare a vertex that is sent to the neighbour
   *
   * We run the $2^d$ adjacent cells and for each cell that is local, we do
   * send its d boundary values that are adjacent to the vertex away.
   *
   * This algorithm translates into a loop very similar to
   * RiemannSolver::touchVertexFirstTime():
   *
   * - Run through all the $2^d$ adjacent cells. Only those that belong to
   *   toRank are of interest. Skip the others. See below for remarks on
   *   interest.
   * - For any cell assigned to toRank, there are d faces that are adjacent to
   *   vertex.
   * - Get the heap indices of all the surrounding cells. Not that some of
   *   them, by definition, are remote.
   *
   * When we run through all the cells adjacent to a vertex, we may communicate
   * only local cells to other ranks. This defines the first two entries in the
   * corresponding if statement. Furthermore, only those cell pairs sharing a
   * face do exchange data. This is done in the third line. The Manhattan
   * distance of the two entries has to be exactly one. Finally (notably on rank
   * 0), we may only send out data if the corresponding cell is inside the
   * domain.
   *
   * <h2>Enumeration</h2>
   *
   * The faces are enumerated: left, right, bottom, top, front, back. Let n be
   * the normal of the face starting with 0. Then, the face index is 2i+f where
   * f means whether it is the face that runs through the left bottom corner (0)
   * or not (1).
   *
   *
   * <h2>MPI administration</h2>
   *
   * Please note that the communication specification deploys the whole heap
   * management to the Peano kernel. There is thus no need to invoke any
   * MPI-specific heap communication operation.
   */
  void prepareSendToNeighbour(exahype::Vertex& vertex, int toRank,
                              const tarch::la::Vector<DIMENSIONS, double>& x,
                              const tarch::la::Vector<DIMENSIONS, double>& h,
                              int level);


  ///////////////////////////////////////
  // WORKER->MASTER->WORKER (Reduction+Broadcast)
  ///////////////////////////////////////
  /**
   * This routine is called on the worker.
   *
   * Send the local array of minimal time step sizes up to the master. This is
   * one MPI_Send on the whole array.
   */
  void prepareSendToMaster(
      exahype::Cell& localCell, exahype::Vertex* vertices,
      const peano::grid::VertexEnumerator& verticesEnumerator,
      const exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * The reduction-broadcast process consists of a send of time step
   * data from the workers to the master which merges=reduces the time step data and then
   * sends the reduced time step data back to workers which merge it finally.
   * Here, we adopt a worker centric point of view, i.e., the reducing of
   * time step data by the master belongs to the sending process
   * of its workers.
   */
  void mergeWithMaster(
      const exahype::Cell& workerGridCell,
      exahype::Vertex* const workerGridVertices,
      const peano::grid::VertexEnumerator& workerEnumerator,
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      int worker, const exahype::State& workerState,
      exahype::State& masterState);

  /**
   * TODO(Dominic): Revise documentation.
   *
   * We use this hook to send solver information to a worker. It is important
   * that we return false here. Whether or not we need data back from the worker
   * is solely up to prepareSendToWorker(). This
   * mapping determines whether we do sync ranks or not.
   *
   * Exchange all the global states with the worker
   *
   * This ensures that all workers have correct states with correct time step
   * sizes.
   *
   * Through the result of this routine, we can skip worker-master data
   * transfer as all other mappings return false. Such a skip is advantageous
   * if the runner has decided to trigger multiple grid traversals in one
   * batch. This in turn automatically disables the load balancing.
   *
   * Our strategy thus is as follows: If we may skip the Sending, i.e. the
   * user has enabled this optimisation in the ExaHyPE spec file, then we
   * return false if load balancing is disabled.
   */
  bool prepareSendToWorker(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      int worker);

  /**
   * Nop.
   */
  void receiveDataFromMaster(
      exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
      const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
      exahype::Vertex* const receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
      exahype::Cell& receivedCoarseGridCell,
      exahype::Vertex* const workersCoarseGridVertices,
      const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
      exahype::Cell& workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Nop.
   */
  void mergeWithNeighbour(exahype::Vertex& vertex,
                          const exahype::Vertex& neighbour, int fromRank,
                          const tarch::la::Vector<DIMENSIONS, double>& x,
                          const tarch::la::Vector<DIMENSIONS, double>& h,
                          int level);

  /**
   * Nop.
   */
  void prepareCopyToRemoteNode(
      exahype::Cell& localCell, int toRank,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);

  /**
   * Nop.
   */
  void prepareCopyToRemoteNode(exahype::Vertex& localVertex, int toRank,
                               const tarch::la::Vector<DIMENSIONS, double>& x,
                               const tarch::la::Vector<DIMENSIONS, double>& h,
                               int level);

  /**
   * Nop.
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
      int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
      const tarch::la::Vector<DIMENSIONS, double>& h, int level);
  /**
   * Nop.
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
      int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);




  /**
   * Nop.
   */
  void mergeWithWorker(exahype::Cell& localCell,
                       const exahype::Cell& receivedMasterCell,
                       const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
                       const tarch::la::Vector<DIMENSIONS, double>& cellSize,
                       int level);
  /**
   * Nop.
   */
  void mergeWithWorker(exahype::Vertex& localVertex,
                       const exahype::Vertex& receivedMasterVertex,
                       const tarch::la::Vector<DIMENSIONS, double>& x,
                       const tarch::la::Vector<DIMENSIONS, double>& h,
                       int level);
#endif
   /**
    * Nop.
    */
   virtual ~Sending();
 #if defined(SharedMemoryParallelisation)
   /**
    * Nop.
    */
   void mergeWithWorkerThread(const Sending& workerThread);
 #endif
   /**
    * Nop.
    */
   void createInnerVertex(
       exahype::Vertex& fineGridVertex,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
   /**
    * Nop.
    */
   void createBoundaryVertex(
       exahype::Vertex& fineGridVertex,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
   /**
    * Nop.
    */
   void createHangingVertex(
       exahype::Vertex& fineGridVertex,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
   /**
    * Nop.
    */
   void destroyHangingVertex(
       const exahype::Vertex& fineGridVertex,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
   /**
    * Nop.
    */
   void destroyVertex(
       const exahype::Vertex& fineGridVertex,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
   /**
    * Nop.
    */
   void createCell(
       exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
       const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);
   /**
    * Nop.
    */
   void destroyCell(
       const exahype::Cell& fineGridCell,
       exahype::Vertex* const fineGridVertices,
       const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);
  /**
   * Nop.
   */
  void touchVertexFirstTime(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
  /**
   * Nop.
   */
  void touchVertexLastTime(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

  /**
   * Nop.
   */
  void descend(
      exahype::Cell* const fineGridCells,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell);
  /**
   * Nop.
   */
  void ascend(exahype::Cell* const fineGridCells,
              exahype::Vertex* const fineGridVertices,
              const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
              exahype::Vertex* const coarseGridVertices,
              const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
              exahype::Cell& coarseGridCell);
};
#endif
