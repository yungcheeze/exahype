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
 * TODO(Comment):
 *
 * <h2>Shared Memory</h2>
 *
 * As the mapping accesses the state data in a read-only fashion, no special
 * attention is required here.
 *
 * <h2>Bug fixes and optimisations we have not yet performed<h2>
 *
 * 1. Solved: We currently send out and receive 2^{d-1} messages per face. This will become an even bigger issue
 *    when we introduce space-time face data.
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
   * A semaphore that is locked if a thread performs a restriction
   * operation.
   */
  static tarch::multicore::BooleanSemaphore _semaphoreForRestriction;

#if defined(Debug)
  /**
   * Some counters for debugging purposes.
   */
  static int _parentOfCellOrAncestorNotFound;
  static int _parentOfCellOrAncestorFound;
  static int _parentOfDescendantFound;
#endif

  /**
   * Computes the space-time predictor quantities, extrapolates fluxes
   * and (space-time) predictor values to the boundary and
   * computes the volume integral.
   *
   * TODO(Dominic): Dedicate each thread a fixed size Prediction and
   * PredictionVolumeFlux field. No need to store these massive
   * quantities on the heap for each cell.
   */
  static void computePredictionAndVolumeIntegral(exahype::records::ADERDGCellDescription& p);

  /**
   * Sets the extrapolated predictor and fluctuations of an Ancestor to zero.
   */
  static void prepareAncestor(
      exahype::records::ADERDGCellDescription& p);


  /**
   * Prolongates face data from a parent ADERDGCellDescription to
   * \p cellDescription if the fine grid cell associated with
   * \p cellDescription is adjacent to a boundary of the
   * coarse grid cell associated with the parent cell description.
   *
   * \note This function assumes a top-down traversal of the grid and must thus
   * be called from the enterCell(...) or descend(...) functions.
   *
   * \note This method makes only sense if \p cellDescription is of type
   * Descendant and EmptyDescendant.
   */
  static void prolongateADERDGFaceData(
      const exahype::records::ADERDGCellDescription& cellDescription,
      const int parentIndex,
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex);

  /**
   * Prolongates the boundary layer (or the complete solution) from a parent FiniteVolumesCellDescription to
   * \p cellDescription if the fine grid cell associated with
   * \p cellDescription is adjacent to a boundary of the
   * coarse grid cell associated with the parent cell description.
   *
   * \note This function assumes a top-down traversal of the grid and must thus
   * be called from the enterCell(...) or descend(...) functions.
   *
   * \note This method makes only sense if \p cellDescription is of type
   * Descendant and EmptyDescendant.
   */
  static void prolongateFiniteVolumesFaceData(
      const exahype::records::FiniteVolumesCellDescription& cellDescription,
      const int parentIndex,
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex);

  /**
   * Restricts face data from \p cellDescription to
   * a parent ADERDGCellDescription if the fine grid cell associated with
   * \p cellDescription is adjacent to a boundary of the
   * coarse grid cell associated with the parent cell description.
   *
   * \note This function assumes a bottom-up traversal of the grid and must thus
   * be called from the leaveCell(...) or ascend(...) functions.
   *
   * \note This method makes only sense if \p cellDescription is of type Cell.
   *
   * \note We use a semaphore to make this operation thread-safe.
   */
  static void restrictADERDGFaceData(
      const exahype::records::ADERDGCellDescription& cellDescription,
      const int parentIndex,
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex);

  /**
   * Restricts the finite volume solution values (volume data) from \p cellDescription to
   * a parent FiniteVolumesCellDescription if the fine grid cell associated with
   * \p cellDescription is adjacent to a boundary of the
   * coarse grid cell associated with the parent cell description.
   *
   * \note This function assumes a bottom-up traversal of the grid and must thus
   * be called from the leaveCell(...) or ascend(...) functions.
   *
   * \note This method makes only sense if \p cellDescription is of type Cell.
   */
  static void restrictFiniteVolumesSolution(
      const exahype::records::FiniteVolumesCellDescription& cellDescription,
      const int parentIndex,
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex);

  /**
   * Picks out the subcell indices that are not at position \p d.
   */
  static tarch::la::Vector<DIMENSIONS - 1, int> getSubfaceIndex(
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
      const int d);


  #ifdef Parallel
  /**
   * We need read access to the state.
   *
   * \note While currently no forking and joining will happen after the initial
   * grid setup, this will happen in the dynamic AMR case. It thus doesn't hurt to
   * use work with the _state already in the static AMR case.
   * TODO(Dominic): Delete this comment as soon as we have dynamic AMR+MPI.
   *
   * \note Do not set values on the state. Write access to the state is not thread-safe.
   */
  exahype::State* _state;
  /**
   * Tag that is used to exchange all the solver instances in MPI
   */
  static int _mpiTag;

  /**
   * Count the listings of remote ranks that share a vertex
   * adjacent to the face \p faceIndex of a cell.
   * This value is either 0 or 2^{d-1}.
   *
   * If we count 2^{d-1} listings, this implies that this rank
   * shares a whole face with a remote rank.
   *
   * More interestingly, we know from this number how
   * many vertices will try to exchange neighbour information
   * that is related to this face.
   *
   * @developers:
   * TODO(Dominic): We currently check for uniqueness of the
   * remote rank. This might however not be necessary.
   */
  static int countListingsOfRemoteRankAtFace(
      const int faceIndex,
      exahype::Vertex* const verticesAroundCell,
      const peano::grid::VertexEnumerator& verticesEnumerator);

  /**
   * Checks for all cell descriptions (ADER-DG, FV, ...)
   * corresponding to the heap index \p cellDescriptionsIndex
   * if now is the time to send out face data to a
   * neighbouring rank.
   *
   * <h2>Details<\h2>
   * On every cell description, we hold a field of 2*d
   * counters. If a face is part of the MPI boundary,
   * we initialise the corresponding counter with
   * value 2^{d-1}.
   *
   * In the mergeWithNeighbour(...) routine,
   * we then decrement the counters for the face
   * every time one of the 2^{d-1}
   * adjacent vertices touches the face.
   *
   * If the counters hold the value zero, this function returns true.
   * Otherwise it returns false
   *
   * @see decrementCounters
   */
  static bool needToSendFaceData(
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      int cellDescriptionsIndex);


  /**
   * Every call of this function decrements the
   * faceDataExchangeCounter for the face corresponding
   * to the source and destination position pair \p src and \p dest
   * for all cell descriptions corresponding to \p cellDescriptionsIndex.
   *
   * @see hasToSendFace
   */
  static void decrementCounters(
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      int cellDescriptionsIndex);

  /**
   * TODO(Dominic): Docu.
   *
   * TODO(Dominic): Make messaging solver functionality?
   *
   * \note Not thread-safe.
   */
  static void sendADERDGFaceData(
      int toRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      int level,
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      int srcCellDescriptionIndex,
      int destCellDescriptionIndex);
  #endif


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
   * Enter a cell
   *
   * This method first synchronises the time step sizes and time stamps, and
   * then resets the Riemann solve flags and the face data exchange counter for all
   * solvers for which a valid cell description was registered on this cell.
   *
   * Directly after, it runs through all solvers assigned to a cell and invoke the solver's
   * spaceTimePredictor(...) as well as the solver's volumeIntegral(...) routines if
   * the the fine grid cell functions as a compute cell (Cell) for the solver.
   * Please see the discussion in the class header.
   *
   * <h2>MPI</h2>
   *
   * The predictor writes something into the boundary arrays lQhbnd and lFhbnd.
   * These values are write-only and have to exchanged with neighbouring ranks.
   *
   * <h2>AMR</h2>
   *
   * If the fine grid cell functions as a helper cell of type Descendant for a solver, this method invokes
   * prolongateFaceData(...).
   * Further clears the face data of helper cells of type Ancestor.
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
   * If the fine grid cell functions as Cell or Ancestor for a solver, this function method
   * restrictFaceData(...).
   */
  void leaveCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);


#if defined(SharedMemoryParallelisation)
  /**
   * Access the states mapping in a read-only fashion.
   */
  Prediction(const Prediction& masterThread);
#endif

#if defined(SharedMemoryParallelisation)
  /**
   * Nop
   */
  void mergeWithWorkerThread(const Prediction& workerThread);
#endif

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

  /**
   * Exchange all the global states with the worker
   *
   * This ensures that all workers have correct states with correct time step
   * sizes.
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
   * See prepareSendToWorker(). This is the counterpart operation.
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



  //
  // Below all methods are nop.
  //
  //===================================



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
  Prediction();

  /**
   * Nop
   */
  virtual ~Prediction();

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
