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
 
#ifndef EXAHYPE_MAPPINGS_RiemannSolver_H_
#define EXAHYPE_MAPPINGS_RiemannSolver_H_

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
class RiemannSolver;
}
}

/**
 * Run the Riemann solves on the faces
 *
 *
 * @todo Dominic, bitte ordentlich dokumentieren, was hier wann, wo und warum passiert.
 * Bitte auch dokumentieren, falls Du was mal probiert hast und es nicht funktioniert hat.
 *
 *
 * <h2>Synchronisation in multi-rank environment</h2>
 *
 * The individual Riemann solves have to have knowledge about
 * admissible time step sizes. This information is read from the solver
 * instances. In return, they update the solver instances with the maximal
 * permitted time step size for the subsequent iteration. This synchronisation
 * is realised in the mapping NewTimeStep, i.e. no solver synchronisation is to
 * be done in this mapping. It is also realised in SpaceTimePredictor. So
 * please consult the documentation there.
 *
 * @author Dominic E. Charrier and Tobias Weinzierl
 */
class exahype::mappings::RiemannSolver {
 private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log _log;

#ifdef Debug
  /*
   *  Counter for the interior face solves for debugging purposes.
   */
  int _interiorFaceSolves;
  /*
   *  Counter for the boundary face solves for debugging purposes.
   */
  int _boundaryFaceSolves;
#endif

  /**
   * Solve the Riemann problem at the interface between two cells ("left" and
   * "right"). This method only performs a Riemann solve if at least one of the
   * cell descriptions (per solver) associated with the two cells is of type
   * ::Cell and none of the two cells belongs to the boundary.
   * In case a Riemann problem is solved,
   * the method further sets the ::riemannSolvePerformed
   * flags for the particular faces on both cell descriptions (per solver).
   *
   * This method further synchronises the ADERDGCellDescription
   * with the corresponding solver if this is required by the time stepping scheme.
   * This operation must be performed in mergeWithNeighbour(...) and
   * touchVertexFirstTime(...) since both callbacks touch the
   * ADERDGCellDescriptions before the other callbacks.
   *
   * <h2>Rationale</h2>
   *
   * We did originally split up the boundary condition handling and the Riemann
   * updates into two mappings. This offers a functional decomposition. However,
   * both mappings then need a significiant number of technical administrative
   * code (cmp all the loops in touchVertexFirstTime and the redundant code to
   * manage the semaphores). We thus decided to merge both aspects. This also
   * should make sense from a performance point of view.
   *
   * We could potentially remove the face indices here if we had normals that
   * point outwards. However, we don't evaluate the direction of the normal and
   * thus need these counters as a Riemann problem on a face either could be
   * triggered by the left cell or by the right cell.
   *
   * \note The current implementation might classify cells with vertices that are part of the
   * boundary of the domain or outside to be classified as inside of the domain (volume-ratio based).
   *
   * \note We cannot solely check for indices of value
   * multiscalelinked::HangingVertexBookkepper::DomainBoundaryAdjacencyIndex
   * in vertex.getCellDescriptions() to determine if we are on the boundary of the domain
   * since these values are overwritten by multiscalelinked::HangingVertexBookkepper::RemoteAdjacencyIndex
   * if the domain boundary aligns with an MPI boundary
   * (see multiscalelinkedcell::HangingVertexBookkeeper::updateCellIndicesInMergeWithNeighbour(...)).
   *
   * \note Not thread-safe.
   *
   * @param[in] cellDescriptionIndexOfLeftCell
   * @param[in] cellDescriptionIndexOfRightCell
   * @param[in] faceIndexForLeftCell    The index of the interface
   *                                    from the perspective of the "left" cell.
   *                                     One out of
   *                                    (EXAHYPE_FACE_LEFT=0,EXAHYPE_FACE_RIGHT=1,...,EXAHYPE_FACE_TOP=5).
   * @param[in] faceIndexForRightCell   The index of the interface from the
   *                                    perspective of the "right" cell.
   * @param[in] normalNonZero           Non zero component of the
   *                                    normal vector orthogonal to the
   *                                    interface.
   */
  void solveRiemannProblemAtInterface(const int cellDescriptionIndexOfLeftCell,
                                      const int cellDescriptionIndexOfRightCell,
                                      const int faceIndexForLeftCell,
                                      const int faceIndexForRightCell,
                                      const int normalNonZero);
  /**
   * Apply the boundary conditions at the face with index \p faceIndex.
   *
   * This method further synchronises the ADERDGCellDescription
   * with the corresponding solver if this is required by the time stepping scheme.
   * This operation must be performed in mergeWithNeighbour(...) and
   * touchVertexFirstTime(...) since both callbacks touch the
   * ADERDGCellDescriptions before the other callbacks.
   *
   * \note Not thread-safe.
   *
   * @param[in] cellDescription         The cell description
   * @param[in] faceIndex               The index of the interface
   *                                    from the perspective of the cell/cell
   * description. One out of
   *                                    (EXAHYPE_FACE_LEFT=0,EXAHYPE_FACE_RIGHT=1,...,EXAHYPE_FACE_TOP=5).
   * @param[in] normalNonZero           Non zero component of the
   *                                    normal vector orthogonal to the
   *                                    interface.
   * \note Not thread-safe.
   */
  void applyBoundaryConditions(records::ADERDGCellDescription& cellDescription,
                               const int faceIndex, const int normalNonZero);

#ifdef Parallel
  /**
   * Single-sided version of the other solveRiemannProblemAtInterface(). It
   * works only on one cell and one solver within this cell and in return
   * hands in the F and Q values explicitly through  indexOfQValues and
   * indexOfFValues. The Riemann solver is invoked and the bits are set
   * accordingly no matter of what they did hold before, i.e. different to
   * the standard solveRiemannProblemAtInterface() operation, we do not
   * check whether we shall run a Riemann solver or not.
   *
   * This method further synchronises the ADERDGCellDescription
   * with the corresponding solver if this is required by the time stepping scheme.
   * This operation must be performed in mergeWithNeighbour(...) and
   * touchVertexFirstTime(...) since both callbacks touch the
   * ADERDGCellDescriptions before the other callbacks.
   *
   * \note Not thread-safe.
   */
//  static void solveRiemannProblemAtInterface(
//      records::ADERDGCellDescription& cellDescription,
//      const int faceIndexForCell,
//      const int normalNonZero,  // TODO(Tobias): is redundant. We should be able to // derive this from faceIndexForCell
//      const int indexOfQValues, const int indexOfFValues);

  /**
   * Iterates over the received metadata and every time
   * we find a valid entry, we call mergeWithNeighbourData
   * on the solver corresponding to the metadata.
   * if we want to receive the neighbour data
   * or if we just want to drop it.
   * \note Not thread-safe.
   */
  static void mergeWithNeighbourData(
      const int fromRank,
      const int srcCellDescriptionIndex,
      const int destCellDescriptionIndex,
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int level,
      const exahype::MetadataHeap::HeapEntries& receivedMetadata);

  /**
   * Iterates over the received metadata and
   * drop the received neighbour data.
   *
   * \note Not thread-safe.
   */
  static void dropNeighbourData(
      const int fromRank,
      const int srcCellDescriptionIndex,
      const int destCellDescriptionIndex,
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int level,
      const exahype::MetadataHeap::HeapEntries& receivedMetadata);
#endif

 public:
  /**
   * Call the touch vertex first time event on every vertex of
   * the grid. Run in parallel but avoid fine grid races.
   *
   * <h2>Shared Memory</h2>
   * The AvoidFineGridRaces multithreading specification prevents that
   * more than one threads write data for the same face of the grid
   * at the same time.
   *
   * The specification is realised by touching the vertices in
   * a red-black (X and O in the figure below) manner:
   * We might first process the X-vertices in parallel
   * and then the O-vertices.
   * In each step, faces adjacent to the vertices,
   * do not overlap and race conditions can thus
   * not occur.
   *
   *     |    |
   *     |    |
   * ----O----X-----
   *     |    |
   *     |    |
   * ----X----O-----
   *     |    |
   *     |    |
   *
   * TODO(Dominic): It might be useful to introduce a multithreading specification
   * "AvoidFineGridRacesOnlyRed" that processes only the red
   * vertices and not the black ones. De facto, the second sweep only
   * finds the riemannSolvePerfomed flags set and does nothing in
   * our current implementation.
   */
  static peano::MappingSpecification touchVertexFirstTimeSpecification();
  /**
   * @todo Dominic, warum darf das Zeugs so stehen?
   */
  static peano::MappingSpecification enterCellSpecification();

  /**
   * Nop
   */
  static peano::MappingSpecification touchVertexLastTimeSpecification();
  /**
   * Nop
   */
  static peano::MappingSpecification leaveCellSpecification();

  /**
   * Nop
   */
  static peano::MappingSpecification ascendSpecification();

  /**
   * Nop
   */
  static peano::MappingSpecification descendSpecification();

  /**
   * The mapping does synchronise through synchroniseTimeStepping() invoked
   * on the solvers. Yet, no data is transported through the vertices, the
   * cell or the state object.
   *
   * Though we need valid solvers for the actual Riemann solves, we do not
   * do any solver exchange in this mapping. The appropriate data exchange
   * is done in NewTimeStep() and SolutionUpdate().
   */
  static peano::CommunicationSpecification communicationSpecification();

  /**
   * Solve Riemann problems on all interior faces that are adjacent
   * to this vertex and impose boundary conditions on faces that
   * belong to the boundary. This is done for all cell descriptions
   * belonging to the cells that are an interior face.
   *
   * The routine itself runs the loop over the faces. The actual
   * functionality is outsourced to solveRiemannProblemAtInterface().
   *
   * The function ensures implicitly that interior faces
   * do not align with MPI boundaries. In this case, no operation
   * is performed.
   *
   * This method sets the riemannSolvePerformed flag on a cell description
   * if boundary conditions have been imposed for this cell description.
   * This method sets the riemannSolvePerformed flag on both cell descriptions
   * (per solver) for interior faces if a Riemann solve has been performed for
   * both cell descriptions.
   *
   * \note The function itself is not thread-safe.
   * Thread-safety of this function is ensured by setting
   * RiemannSolver::touchVertexFirstTimeSpecification()
   * to peano::MappingSpecification::AvoidFineGridRaces.
   *
   * <h2>Limiter identification</h2>
   * Each ADER-DG solver analyses the local min and max values within a cell.
   * This information however is not stored in the cell but on the 2d faces
   * of a cell. See Cell::setSolutionMinMaxAndAnalyseValidity() for details.
   * The face information then in the subsequent step has to be merged which
   * is done when we trigger the Riemann solve.
   *
   * @see Cell::mergeSolutionMinMaxOnFace()
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
   * Resets counters that are solely used for debugging.
   */
  void beginIteration(exahype::State& solverState);

  /**
   * Prints the output of counters that are solely used for debugging.
   */
  void endIteration(exahype::State& solverState);

#ifdef Parallel

  /**
   * Merge vertex with the incoming vertex from a neighbouring computation node.
   *
   * When Peano is running in parallel the data exchange is done vertex-wise
   * between two grid iterations, i.e. the predictor sends out data in one step
   * and in the following step we receive this data and merge it into the local
   * Riemann integrals. The routine is thus very similar to
   * touchVertexFirstTime():
   *
   * - We identify incoming faces.
   * - We create (temporary) indices on the heap.
   * - We receive the data into these indices.
   *
   * As we always receive data in the iteration following the sends and as
   * Peano inverts the traversal direction after each grid sweep, we have to
   * invert the order in which data is received, too.
   *
   * \note It happens or is possible that this operation is performed
   * after touchVertexLastTime(...) was invoked.
   *
   * <h2>Send</h2>
   * Sending out data corresponds logically to a projection of cell data onto
   * the faces. Therefore, I realise it within
   * SpaceTimePredictor::prepareSendToNeighbour().
   *
   *
   * <h2>Min max analysis</h2>
   * The min/max analysis runs analogously. We do send out min and max from
   * either side to the other rank and then merge min and max on both sides
   * into the local data.
   *
   * <h2>Face data exchange counters</h2>
   * Cell::isInside() does not imply that all adjacent vertices are
   * inside. If we count down the counter only on
   * vertices that are inside we might not send out all faces
   * of a cell that is close to the boundary.
   * @see Prediction::countListingsOfRemoteRankByInsideVerticesAtFace.
   */
  void mergeWithNeighbour(exahype::Vertex& vertex,
                          const exahype::Vertex& neighbour, int fromRank,
                          const tarch::la::Vector<DIMENSIONS, double>& x,
                          const tarch::la::Vector<DIMENSIONS, double>& h,
                          int level);
  /**
   * Receive kick-off message from master
   *
   * Counterpart of prepareSendToWorker(). This operation is called once when
   * we receive data from the master node.
   *
   * To do a proper Riemann solve, it is important that all the solvers have
   * the right state. We therefore receive for each individual solver a couple
   * of messages and merge them into the global instance before we continue
   * with the actual iteration. However, this data exchange in the mapping
   * SpaceTimePredictor. See the documentation there for rationale.
   *
   * @see SpaceTimePredictor::receiveDataFromMaster()
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
  void prepareCopyToRemoteNode(
      exahype::Cell& localCell, int toRank,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);
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
  RiemannSolver();
#if defined(SharedMemoryParallelisation)
  /**
   * Nop.
   */
  RiemannSolver(const RiemannSolver& masterThread);
#endif
  /**
   * Nop.
   */
  virtual ~RiemannSolver();
#if defined(SharedMemoryParallelisation)
  /**
   * Nop.
   */
  void mergeWithWorkerThread(const RiemannSolver& workerThread);
#endif
  /**
   * Nop
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
