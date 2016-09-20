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
 
#ifndef EXAHYPE_MAPPINGS_Prediction_H_
#define EXAHYPE_MAPPINGS_Prediction_H_

#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include "peano/CommunicationSpecification.h"
#include "peano/MappingSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "tarch/multicore/MulticoreDefinitions.h"

#include "tarch/multicore/BooleanSemaphore.h"

#include "exahype/Cell.h"
#include "exahype/State.h"
#include "exahype/Vertex.h"

namespace exahype {
  namespace mappings {
    class Prediction;
  }
}

/**
 * This mapping realises the space-time predictor
 *
 * We do run through all the cells and, per solver, run the space-time
 * prediction and the volume integral.
 * Most of the interesting stuff is done in enterCell().
 *
 * The other methods despite leaveCell(...) as well as MPI- and Shared-Mem-specific methods
 * perform no operation.
 *
 * This mapping's enterCell(...) and leaveCell(...) methods are further used to
 * prolongate coarse grid face unknowns down to fine grid cells and to
 * restrict coarse grid face unknowns up to coarse grid cells.
 *
 * As all state data is encoded in the global solver states and as all data
 * accesses are read-only, the mapping has no object attributes.
 *
 * <h2>MPI</h2>
 *
 * For a valid computation of the space-time prediction, we need to know the
 * correct time step on each rank. We thus distribute the all solvers among the
 * ranks when we start them up. See prepareSendToWorker() and the corresponding
 * receiveDataFromMaster().
 *
 * <h2>Shared Memory</h2>
 *
 * As the mapping accesses the state data in a read-only fashion, no special
 * attention is required here.
 *
 * <h2>Optimisations</h2>
 * We dedicate each thread a fixed size space-time predictor,
 * space-time volume flux, predictor, and volume flux
 * field. There is no need to store these massive
 * quantities on the heap for each cell as it was done
 * in the baseline code.
 * This massively reduces the memory footprint of
 * the method and might lead to a more
 * cache-friendly code since we reuse
 * the temporary variables multiple times
 * per grid traversal (unverified).
 *
 * @author Dominic E. Charrier and Tobias Weinzierl
 *
 * @developers:
 * TODO(Dominic): Need to propagate face data
 * from master to worker (prolongation).
 * Correct face data on Ancestor and Descendant is
 * made available by the functions provided in this mapping.
 * This mapping should therefore be merged with the
 * space-time predictor mapping.
 * Need to propagate face data from worker to master
 * (erasing,joins).
 */
class exahype::mappings::Prediction {
 private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log _log;

  /**
   * Temporary variable: Degrees of freedom of the space-time predictor
   * for all solvers in case we parallelise over the cell descriptions associated with a cell.
   */
  double** _lQi=nullptr;

  /**
   * Temporary variable: Old degrees of freedom of the old Picard loop space-time predictor
   * for all solvers in case we parallelise over the cell descriptions associated with a cell.
   * Theses values are only used within the Picard loop of the space-time predictor computation.
   */
  double** _lQi_old=nullptr;

  /**
   * Temporary variable: Degrees of freedom of the Picard loop right-hand side
   * for all solvers in case we parallelise over the cell descriptions associated with a cell.
   * Theses values are only used within the Picard loop of the space-time predictor computation.
   */
  double** _rhs=nullptr;

  /**
   * Temporary variable: Degrees of freedom of the initial Picard loop right-hand side
   * for all solvers in case we parallelise over the cell descriptions associated with a cell.
   * Theses values are only used within the Picard loop of the space-time predictor computation.
   */
  double** _rhs_0=nullptr;

  /**
   * Temporary variable: Degrees of freedom of the space-time predictor volume flux
   * for all solvers in case we parallelise over the cell descriptions associated with a cell.
   */
  double** _lFi=nullptr;
  /**
   * Temporary variable: Degrees of freedom of time averaged space-time predictor
   * for all solvers in case we parallelise over the cell descriptions associated with a cell.
   */
  double** _lQhi=nullptr;
  /**
   * Temporary variable: Degrees of freedom of the time averaged space-time predictor volume flux
   * for all solvers in case we parallelise over the cell descriptions associated with a cell.
   */
  double** _lFhi=nullptr;

  /**
   * Initialises the temporary variables.
   *
   * \note We parallelise over the domain
   * (mapping is copied for each thread) and
   * over the solvers registered on a cell.
   *
   * \note We need to initialise the temporary variables
   * in this mapping and not in the solvers since the
   * solvers in exahype::solvers::RegisteredSolvers
   * are not copied for every thread.
   */
  void initTemporaryVariables();

  /**
   * Deletes the temporary variables.
   *
   * \note We need to initialise the temporary variables
   * in this mapping and not in the solvers since the
   * solvers in exahype::solvers::RegisteredSolvers
   * are not copied for every thread.
   */
  void deleteTemporaryVariables();

 public:
  static peano::MappingSpecification touchVertexLastTimeSpecification();
  static peano::MappingSpecification touchVertexFirstTimeSpecification();
  static peano::MappingSpecification enterCellSpecification();
  static peano::MappingSpecification leaveCellSpecification();
  static peano::MappingSpecification ascendSpecification();
  static peano::MappingSpecification descendSpecification();

  /**
   * Please consult the specification's documentation in NewTimeStep.
   */
  static peano::CommunicationSpecification communicationSpecification();

  /**
   * Initialise the temporary variables (of the master thread
   * in a shared memory build).
   */
  Prediction();

  /**
   * Delete the temporary variables (for
   * master and worker threads in a shared memory build).
   */
  ~Prediction();

  #if defined(SharedMemoryParallelisation)
  /**
   * Initialise the temporary variables of a
   * worker thread.
   */
  Prediction(const Prediction& masterThread);
  /**
   * Nop
   */
  void mergeWithWorkerThread(const Prediction& workerThread);
  #endif

  /**
   * This method first synchronises the time step sizes and time stamps, and
   * then resets the Riemann solve flags and the face data exchange counter for all
   * solvers for which a valid cell description was registered on this cell.
   *
   * Directly after, it runs through all solvers assigned to a cell and invoke the solver's
   * spaceTimePredictor(...) as well as the solver's volumeIntegral(...) routines if
   * the the fine grid cell functions as a compute cell (Cell) for the solver.
   * Please see the discussion in the class header.
   *
   * @see enterCellSpecification()
   */
  void enterCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Nop.
   */
  void leaveCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Nop
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
   * Nop
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
   * Nop
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
   * Nop
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
   * Nop
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
   * Nop
   */
  void createCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Nop
   */
  void destroyCell(
      const exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  //
  // Below all methods are nop.
  //
  //===================================
#ifdef Parallel
  /*
   * Nop
   */
  void prepareSendToNeighbour(exahype::Vertex& vertex, int toRank,
                              const tarch::la::Vector<DIMENSIONS, double>& x,
                              const tarch::la::Vector<DIMENSIONS, double>& h,
                              int level);

  /*
   * Nop
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
   * Nop
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
   * Nop
   */
  void mergeWithNeighbour(exahype::Vertex& vertex,
                          const exahype::Vertex& neighbour, int fromRank,
                          const tarch::la::Vector<DIMENSIONS, double>& x,
                          const tarch::la::Vector<DIMENSIONS, double>& h,
                          int level);

  /**
   * Nop
   */
  void prepareCopyToRemoteNode(exahype::Vertex& localVertex, int toRank,
                               const tarch::la::Vector<DIMENSIONS, double>& x,
                               const tarch::la::Vector<DIMENSIONS, double>& h,
                               int level);

  /**
   * Nop
   */
  void prepareCopyToRemoteNode(
      exahype::Cell& localCell, int toRank,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);

  /**
   * Nop
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
      int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
      const tarch::la::Vector<DIMENSIONS, double>& h, int level);

  /**
   * Nop
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
      int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);

  /**
   * Nop
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
   * Nop
   */
  void prepareSendToMaster(
      exahype::Cell& localCell, exahype::Vertex* vertices,
      const peano::grid::VertexEnumerator& verticesEnumerator,
      const exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Nop
   */
  void mergeWithWorker(exahype::Cell& localCell,
                       const exahype::Cell& receivedMasterCell,
                       const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
                       const tarch::la::Vector<DIMENSIONS, double>& cellSize,
                       int level);

  /**
   * Nop
   */
  void mergeWithWorker(exahype::Vertex& localVertex,
                       const exahype::Vertex& receivedMasterVertex,
                       const tarch::la::Vector<DIMENSIONS, double>& x,
                       const tarch::la::Vector<DIMENSIONS, double>& h,
                       int level);
#endif
  /**
   * Nop
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
   * Nop
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
   * Nop
   */
  void beginIteration(exahype::State& solverState);

  /**
   * Nop
   */
  void endIteration(exahype::State& solverState);

  /**
   * Nop
   */
  void descend(
      exahype::Cell* const fineGridCells,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell);

  /**
   * Nop
   */
  void ascend(exahype::Cell* const fineGridCells,
              exahype::Vertex* const fineGridVertices,
              const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
              exahype::Vertex* const coarseGridVertices,
              const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
              exahype::Cell& coarseGridCell);
};

#endif
