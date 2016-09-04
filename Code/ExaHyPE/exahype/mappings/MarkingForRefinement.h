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
 
#ifndef EXAHYPE_MAPPINGS_MarkingForRefinement_H_
#define EXAHYPE_MAPPINGS_MarkingForRefinement_H_

#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include "peano/CommunicationSpecification.h"
#include "peano/MappingSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "tarch/multicore/MulticoreDefinitions.h"

#include "exahype/Cell.h"
#include "exahype/State.h"
#include "exahype/Vertex.h"

#include "peano/utils/Globals.h"

namespace exahype {
namespace mappings {
class MarkingForRefinement;
}
}

/**
 * This mapping is one of the two mappings dealing with
 * the refinement and erasing of compute cells.
 *
 * This mapping is only used for marking compute cells as candidates
 * for refinement or erasing. The newly created cells
 * are initialised in mapping Refinement.
 *
 * <h2>MPI<h2>
 * TODO(Dominic): Add docu.
 *
 * @developers:
 * 2. TODO(Dominic): Need to propagate refinement events, cell type
 * changes, and volume data from master to worker (refinement,forks).
 * Need to propagate refinement events, cell type changes,
 * and volume data from worker to master (erasing,joins).
 * 3. TODO(Dominic): The mapping in the current form does not
 * require communication between neighbouring cells. However
 * if we consider non-local refinement criteria
 * as in M. Dumbser's papers, this will be necessary.
 *
 * @author Dominic Etienne Charrier
 */
class exahype::mappings::MarkingForRefinement {
 private:
  /**
   * The logging device of this mapping.
   */
  static tarch::logging::Log _log;

  /**
   * Local copy of the state.
   */
  exahype::State _state;

#ifdef Parallel

  /**
   * Send out ADER-DG solution values to the master or worker
   * if a cell description is of type Cell.
   *
   * If the cell description is of type Ancestor, Descendant
   * EmptyAncestor, EmptyDescendant change the type to Ancestor and Descendant and
   * thus ensure that the cell description holds
   * face data on the boundary. Later during the traversal,
   * the MarkingForAugmentation::enterCell(...) method checks if the cell description is part
   * of the master-worker MPI boundary and does not change the type anymore.
   */
  static void sendADERDGDataToMasterOrWorker(
      int cellDescriptionsIndex,
      int toRank,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);

  /**
   * Send out Finite Volumes solution values to the master or worker
   * if a cell description is of type Cell.
   *
   * If the cell description is of type Ancestor, Descendant
   * EmptyAncestor, EmptyDescendant
   * change the type to Ancestor and Descendant and
   */
  static void sendFiniteVolumesDataToMasterOrWorker(
      int cellDescriptionsIndex,
      int toRank,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);


  /**
   * Sets heap indices of all cell descriptions (ADER-DG, FV, ...) that were received due to
   * a fork or join event to multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex,
   * and the parent index of the cell descriptions to the specified \p parentIndex.
   */
  static void resetHeapIndices(
      const int cellDescriptionsIndex,
      const int parentIndex);

  /**
   * Receives ADER-DG solution values from the master or worker
   * if a cell description is of type Cell.
   *
   * If the cell description is of type Ancestor, Descendant
   * EmptyAncestor, EmptyDescendant
   * change the type to Ancestor and Descendant and
   */
  static void receiveADERDGDataFromMasterOrWorker(
      const int cellDescriptionsIndex,
      const int fromRank,
      const peano::heap::MessageType& messageType,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize,
      const int level,
      const int receivedMetadataIndex);

  /**
   * Returns false if the \p cellDescriptionsIndex is invalid,
   * or if no cell descriptions is registered for this cellDescriptionsIndex,
   * i.e., the vector is empty.
   * Further returns false if the geometry information on the cell descriptions
   * the \p cellDescriptionsIndex is pointing at does not match
   * with \p cellCentre and \p cellSize.
   *
   * Returns true otherwise.
   */
  static bool geometryInfoDoesMatch(
      const int cellDescriptionsIndex,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const int level);

  /**
   * Receives Finite Volumes solution values from the master or worker
   * if a cell description is of type Cell.
   *
   * If the cell description is of type Ancestor, Descendant
   * EmptyAncestor, EmptyDescendant
   * change the type to Ancestor and Descendant and
   */
  static void receiveFiniteVolumesDataFromMasterOrWorker(
      const int cellDescriptionsIndex,
      const int fromRank,
      const peano::heap::MessageType& messageType,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize,
      const int level,
      const int receivedMetadataIndex);
#endif

 public:

  /**
    * Run through the whole grid. Run concurrently on the fine grid.
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
   * Receive data and state before first touch first time.
   * Send data and state fater last touch vertex first time.
   * Further let Peano handle heap data exchange internally.
   */
  static peano::CommunicationSpecification communicationSpecification();


  /**
   * Copy the state and start to send synchronous data.
   */
  void beginIteration(exahype::State& solverState);

  /**
   * Finish sending synchronous data.
   */
  void endIteration(exahype::State& solverState);

  /**
   * If the fine grid cell functions as compute cell of a solver:
   * Mark the compute cell as candidates for refinement or erasing
   * if no other refinement event is set.
   * If the compute cell is marked for refinement, refinement on all
   * adjacent vertices of the fine grid cell is triggered.
   */
  void enterCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

#ifdef Parallel
  /**
   * TODO(Dominic): Implement & add docu.
   *
   * Move data to neighbour
   *
   * Throughout the joins or forks of subdomains, data has to be moved
   * between the nodes. This operation allows the user to plug into these
   * movements. Different to the neighbour communciation, a move always
   * implies that Peano is copying the data structure bit-wise to the
   * remote node. However, if you have heap data and pointers, e.g.,
   * associated to your vertices/cells you have to take care yourself that
   * this data is moved as well.
   *
   * If data is inside your computational domain, you fork, and this data
   * now is located at the new boundary, Peano moves the data as well, i.e.
   * the move also can be a (global) copy operation.
   *
   * @param localCell The local cell. This is not a copy, i.e. you may
   *                    modify the cell before a copy of it is sent away.
   */
  void prepareCopyToRemoteNode(
      exahype::Cell& localCell, int toRank,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);

  /**
   * TODO(Dominic): Implement & add docu.
   *
   * Merge with remote data due to fork or join
   *
   * This operation takes remote data and merges it into the local copy, as
   * data is moved from one rank to another, e.g. Do not use an assignment
   * operator on the whole record, as you may overwrite only PDE-specific
   * fields in localVertex.
   *
   * @param localCell    Local cell data. Some information here is already
   *                     set: Adjacency information from the master, e.g.,
   *                     already is merged. Thus, do not overwrite
   *                     non-PDE-specific data.
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
      int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);

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
  void prepareSendToNeighbour(exahype::Vertex& vertex, int toRank,
                              const tarch::la::Vector<DIMENSIONS, double>& x,
                              const tarch::la::Vector<DIMENSIONS, double>& h,
                              int level);

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
   * Nop.
   */
  void prepareSendToMaster(
      exahype::Cell& localCell, exahype::Vertex* vertices,
      const peano::grid::VertexEnumerator& verticesEnumerator,
      const exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

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
   MarkingForRefinement();

 #if defined(SharedMemoryParallelisation)
   /**
    * Nop.
    */
   MarkingForRefinement(const MarkingForRefinement& masterThread);
 #endif

   /**
    * Nop.
    */
   virtual ~MarkingForRefinement();

 #if defined(SharedMemoryParallelisation)
   /**
    * Nop.
    */
   void mergeWithWorkerThread(const MarkingForRefinement& workerThread);
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
  void leaveCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

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
