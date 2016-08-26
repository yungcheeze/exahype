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
 
#ifndef EXAHYPE_MAPPINGS_MarkingForAugmentation_H_
#define EXAHYPE_MAPPINGS_MarkingForAugmentation_H_

#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include "peano/CommunicationSpecification.h"
#include "peano/MappingSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "tarch/multicore/MulticoreDefinitions.h"

#include "exahype/Cell.h"
#include "exahype/State.h"
#include "exahype/Vertex.h"

#include "exahype/solvers/Solver.h"

#include "peano/utils/Globals.h"

namespace exahype {
namespace mappings {
class MarkingForAugmentation;
}
}

/**
 * TODO(Dominic): Docu.
 *
 * @developers:
 * 1. TODO(Dominic): We need to add a veto that the cell type
 * of an ancestor/descendant can only be Ancestor/Descendant
 * if the cell is a worker root or a master leaf at the fork boundary.
 */
class exahype::mappings::MarkingForAugmentation {
 private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log _log;

//  #define logInfoM(a,b) std::cout << "[dominic-acer],rank:" << tarch::parallel::Node::getInstance().getRank() << " info         exahype::mappings::MarkingForAugmentation::" << a << " " << b << std::endl;

  /**
   * Pointer to the state.
   * We need direct read access to the state to identify during the
   * traversal for which neighbouring ranks forking was triggered
   * and which ranks are forking.
   */
  const exahype::State* _state;

  /**
   * The augmentation control states.
   */
  enum class AugmentationControl {
    /**
     * Indicates that a spacetree cell is next to another spacetree cell
     * of type exahype::records::ADERDGCellDescription::Cell.
     */
    NextToCell = 0,
    /**
     * Indicates that a spacetree cell is next to another spacetree cell
     * of type exahype::records::ADERDGCellDescription::Ancestor or
     * exahype::records::ADERDGCellDescription::EmptyAncestor.
     */
    NextToAncestor = 1,
    /**
     * Indicates that a spacetree cell is both, next to another spacetree cell
     * of type exahype::records::ADERDGCellDescription::Ancestor or
     * exahype::records::ADERDGCellDescription::EmptyAncestor, and
     * a spacetree cell of type
     * exahype::records::ADERDGCellDescription::Cell.
     */
    NextToCellAndAncestor = 2,
    /**
     * Indicates that a spacetree cell is neither, next to another spacetree cell
     * of type exahype::records::ADERDGCellDescription::Ancestor or
     * exahype::records::ADERDGCellDescription::EmptyAncestor, nor
     * next to a spacetree cell of type exahype::records::ADERDGCellDescription::Cell.
     *
     * A cell of type exahype::records::ADERDGCellDescription::Descendant can then request erasing.
     * A cell of type exahype::records::ADERDGCellDescription::Cell does then not need
     * to request augmenting.
     */
    Default = 3
  };

  /**
   * Returns  AugmentationControl::::NextToCell if the cell has a neighbour of
   * type exahype::records::ADERDGCellDescription:Cell,
   * exahype::solvers::Solver::NextToAncestor if the cell has a neighbour of
   * type exahype::records::ADERDGCellDescription:Ancestor or
   * exahype::records::ADERDGCellDescription::EmptyAncestor.
   * Is both the case, this function returns
   * AugmentationControl::::NextToCellOrAncestor.
   * If none of the previous is the case, this function returns
   * AugmentationControl::Default.
   */
  AugmentationControl augmentationCriterion(
      const int solverNumber,
      const exahype::records::ADERDGCellDescription::Type type, const int level,
      const tarch::la::Vector<THREE_POWER_D, int>&
          neighbourCellDescriptionIndices) const;


#ifdef Parallel
  /**
   * TODO(Dominic): Add docu.
   */
  void decodeADERDGMetadataInMergeWithNeigbour(const int destCellDescriptionIndex,const int receivedMetadataIndex) const;

  /**
   * TODO(Dominic): Add docu.
   */
  void decodeFiniteVolumesMetadataInMergeWithNeigbour(const int destCellDescriptionIndex,const int receivedMetadataIndex) const;
#endif

 public:
  /**
   * These flags are used to inform Peano about your operation. It tells the
   * framework whether the operation is empty, whether it works only on the
   * spacetree leaves, whether the operation can restart if the thread
   * crashes (resiliency), and so forth. This information allows Peano to
   * optimise the code.
   *
   * @see peano::MappingSpecification for information on thread safety.
   */
  static peano::MappingSpecification touchVertexLastTimeSpecification();
  static peano::MappingSpecification touchVertexFirstTimeSpecification();
  static peano::MappingSpecification enterCellSpecification();
  static peano::MappingSpecification leaveCellSpecification();
  static peano::MappingSpecification ascendSpecification();
  static peano::MappingSpecification descendSpecification();

  static peano::CommunicationSpecification communicationSpecification();

  /**
   * TODO(Dominic): Add docu.
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
   * TODO(Dominic): Add docu.
   */
  void enterCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * TODO(Dominic): Add docu.
   */
  void leaveCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

#ifdef Parallel
  /**
   * TODO(Dominic): Add docu.
   */
  void mergeWithNeighbour(exahype::Vertex& vertex,
                          const exahype::Vertex& neighbour, int fromRank,
                          const tarch::la::Vector<DIMENSIONS, double>& x,
                          const tarch::la::Vector<DIMENSIONS, double>& h,
                          int level);

  /**
   * TODO(Dominic): Add docu.
   */
  void prepareSendToNeighbour(exahype::Vertex& vertex, int toRank,
                              const tarch::la::Vector<DIMENSIONS, double>& x,
                              const tarch::la::Vector<DIMENSIONS, double>& h,
                              int level);



  //
  // Below all methods are nop.
  //
  // ==================================



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
    MarkingForAugmentation();

  #if defined(SharedMemoryParallelisation)
    /**
     * Nop.
     */
    MarkingForAugmentation(const MarkingForAugmentation& masterThread);
  #endif

    /**
     * Nop.
     */
    virtual ~MarkingForAugmentation();

  #if defined(SharedMemoryParallelisation)
    /**
     * Nop.
     */
    void mergeWithWorkerThread(const MarkingForAugmentation& workerThread);
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
  void beginIteration(exahype::State& solverState);

  /**
   * Nop.
   */
  void endIteration(exahype::State& solverState);

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
