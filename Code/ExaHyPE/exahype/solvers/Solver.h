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
 
#ifndef _EXAHYPE_SOLVERS_SOLVER_H_
#define _EXAHYPE_SOLVERS_SOLVER_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>

#include "peano/utils/Globals.h"

#include "tarch/la/Vector.h"

#include "exahype/profilers/Profiler.h"
#include "exahype/profilers/simple/NoOpProfiler.h"

#include "peano/heap/DoubleHeap.h"

// TODO(Dominic): Change to enum
#define EXAHYPE_FACE_LEFT 0
#define EXAHYPE_FACE_RIGHT 1
#define EXAHYPE_FACE_FRONT 2
#define EXAHYPE_FACE_BACK 3
#define EXAHYPE_FACE_BOTTOM 4
#define EXAHYPE_FACE_TOP 5

// todo 08/02/16:Dominic Etienne Charrier
// move somewhere else
// is evaluated at compile time
constexpr int power(int basis, int exp) {
  return (exp == 0) ? 1 : basis * power(basis, exp - 1);
}

#ifdef ALIGNMENT
constexpr int addPadding(const int originalSize) {
  return ALIGNMENT/8 * static_cast<int>((originalSize+(ALIGNMENT/8-1))/(ALIGNMENT/8));
}
#else
constexpr int addPadding(const int originalSize) {
  return originalSize;
}
#endif

namespace exahype {
  /**
   * Forward declaration of exahype::Cell.
   */
  class Cell;
  /**
   * We store the degrees of freedom associated with the ADERDGCellDescription and FiniteVolumesCellDescription
   * instances on this heap.
   * We further use this heap to send and receive face data from one MPI rank to the other.
   */
  typedef peano::heap::PlainDoubleHeap DataHeap;

  namespace solvers {
    class Solver;

    bool FuseAlgorithmicTimeSteps = false;

    typedef std::vector<Solver*> RegisteredSolversEntries;
    /**
     * All the registered solvers. Has to be declared extern in C++ standard as
     * it is instantiated in the corresponding cpp file.
     */
    // TODO: std::vector<std::unique_ptr<Solver>> ?!
    extern std::vector<Solver*> RegisteredSolvers;
}  // namespace solvers
}  // namespace exahype

/**
 * Describes one solver.
 */
class exahype::solvers::Solver {
 public:
  /**
   * The type of a solver.
   */
  enum class Type { ADER_DG, FiniteVolumes };

  /**
   * The time stepping mode.
   */
  enum class TimeStepping {
    /**
     * In the global time stepping mode, every cells works with the same time step.
     */
    Global,
    /**
     * In the fixed time stepping mode, we assume that each cell advanced in
     * time with the prescribed time step size. No CFL condition is checked.
     */
    GlobalFixed
    // Local, Anarchic
  };

  /**
   * The refinement control.
   */
  enum class RefinementControl { Keep = 0, Refine = 1, Erase = 2 };

  /**
   * TODO(Dominic): Docu.
   */
  typedef struct SubcellPosition {
    int parentIndex;
    tarch::la::Vector<DIMENSIONS, int> subcellIndex;

    SubcellPosition() : parentIndex(-1), subcellIndex(-1) {}
    ~SubcellPosition() {}
  } SubcellPosition;

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
   * Default return value of function getElement(...)
   * If we do not find the element in a vector
   * stored at a heap address.
   */
  static const int NotFound;

 protected:
  /**
   * Each solver has an identifier/name. It is used for debug purposes only.
   */
  const std::string _identifier;

  const Type _type;

  /**
   * The number of state variables of the conservation or balance law.
   */
  const int _numberOfVariables;

  /**
   * The number of parameters, e.g, material parameters.
   */
  const int _numberOfParameters;

  /**
   * The number of nodal basis functions that are employed in each
   * coordinate direction.
   */
  const int _nodesPerCoordinateAxis;

  /**
   * The maximum extent a cell is allowed to have in each coordinate direction.
   */
  const double _maximumMeshSize;

  /**
   * The time stepping mode of this solver.
   */
  const TimeStepping _timeStepping;

  /**
   * A profiler for this solver.
   */
  std::unique_ptr<profilers::Profiler> _profiler;

 public:
  Solver(const std::string& identifier, exahype::solvers::Solver::Type type,
         int numberOfVariables, int numberOfParameters,
         int nodesPerCoordinateAxis, double maximumMeshSize,
         exahype::solvers::Solver::TimeStepping timeStepping,
         std::unique_ptr<profilers::Profiler> profiler =
             std::unique_ptr<profilers::Profiler>(
                 new profilers::simple::NoOpProfiler("")));

  virtual ~Solver() { _profiler->writeToConfiguredOutput(); }

  // Disallow copy and assignment
  Solver(const Solver& other) = delete;
  Solver& operator=(const Solver& other) = delete;

  /**
   * Return a string representation for the type \p param.
   */
  static std::string toString(const exahype::solvers::Solver::Type& param);

  /**
   * Return a string representation for the time stepping mode \p param.
   */
  static std::string toString(const exahype::solvers::Solver::TimeStepping& param);

  /**
   * Returns the maximum extent a mesh cell is allowed to have
   * in all coordinate directions.
   * This maximum mesh size is used both as a
   * constraint on the AMR as well as to set up the initial
   * grid. If you return the extent of the computational domain in
   * each coordinate direction or larger values,
   * you indicate that this solver is not active in the domain.
   */
  double getMaximumMeshSize() const;

  /**
   * Returns the identifier of this solver.
   */
  std::string getIdentifier() const;

  /**
   * Returns the type of this solver.
   */
  Type getType() const;

  /**
   * Returns the time stepping algorithm this solver is using.
   */
  TimeStepping getTimeStepping() const;

  /**
   * Returns the number of state variables.
   */
  int getNumberOfVariables() const;

  /**
   * Returns the number of parameters, e.g.,material constants etc.
   */
  int getNumberOfParameters() const;

  /**
   * If you use a higher order method, then this operation returns the
   * polynomial degree plus one. If you use a Finite Volume method, it
   * returns the number of cells within a patch per coordinate axis.
   */
  int getNodesPerCoordinateAxis() const;

  /**
   * This operation allows you to impose time-dependent solution values
   * as well as to add contributions of source terms.
   * Please be aware that this operation is called per time step if
   * the corresponding predicate hasToUpdateSolution() yields true for the
   * region.
   */
  virtual void solutionAdjustment(
      double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,
      const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) = 0;

  virtual bool hasToAdjustSolution(
      const tarch::la::Vector<DIMENSIONS, double>& center,
      const tarch::la::Vector<DIMENSIONS, double>& dx, double t) = 0;

  /**
   * @defgroup AMR Solver routines for adaptive mesh refinement
   */
  ///@{
  /**
   * The refinement criterion that must be defined by the user.
   *
   */
  // @todo: 16/04/06:Dominic Etienne Charrier Consider to correct the level in
  // the invoking code, i.e., level-> level-1
  // since this is was the user expects.
  virtual exahype::solvers::Solver::RefinementControl refinementCriterion(
      const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,
      const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
      const int level) = 0;

  virtual std::string toString() const;

  virtual void toString(std::ostream& out) const;

  /**
   * Run over all solvers and identify the minimal time stamp.
   */
  virtual double getMinTimeStamp() const = 0;

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  virtual double getMinTimeStepSize() const = 0;

  virtual void updateNextTimeStepSize(double value) = 0;

  virtual void initInitialTimeStamp(double value) = 0;

  /**
   * Copies the time stepping data from the global solver onto the patch's time
   * meta data.
   *
   * \param[in] element Index of the cell description in
   *                    the array at address \p cellDescriptionsIndex.
   */
  virtual void synchroniseTimeStepping(
      const int cellDescriptionsIndex,
      const int element) = 0;

  virtual void startNewTimeStep() = 0;

  virtual double getNextMinTimeStepSize() const =0;

  /**
   * Run over all solvers and identify the minimal time stamp.
   */
  static double getMinSolverTimeStampOfAllSolvers();

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  static double getMinSolverTimeStepSizeOfAllSolvers();

  /**
   * Run over all solvers and identify the maximal time stamp.
   *
   * On the individual patches, we do use the min time stamp so
   * far, so the routine returns the maximum over all min solver
   * time stamps.
   */
  static double getMaxSolverTimeStampOfAllSolvers();

  static bool allSolversUseTimeSteppingScheme(solvers::Solver::TimeStepping scheme);

  static double getCoarsestMeshSizeOfAllSolvers();
  static double getFinestMaximumMeshSizeOfAllSolvers();

  /**
   * Returns true if the index \p cellDescriptionsIndex
   * is a valid heap index.
   */
  virtual bool isValidCellDescriptionIndex(
      const int cellDescriptionsIndex) const = 0;

  /**
   * If an entry for this solver exists,
   * return the element index of the cell description
   * in the array at address \p cellDescriptionsIndex.
   * Otherwise and if \p cellDescriptionsIndex is an
   * invalid index, return Solver::NotFound.
   */
  virtual int tryGetElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const = 0;

  /**
   * Modify a cell description in enter cell event.
   * This event should be used for single cell operations
   * like marking for refinement, erasing, augmenting,
   * or deaugmenting.
   *
   * Returns true if a performed action requires to
   * refine the mesh.
   *
   * \note We use this at the moment only
   * for refinement events. We can consider later
   * on to merge the time stepping functionality
   * (solution update, predictor comp.) into
   * this hook.
   */
  virtual bool enterCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) = 0;

  /**
   * Refinement routine that should be used for
   * collective children-parent operations.
   *
   * Returns true if the grid cell associated with
   * the cell description can be erased due to
   * an erasing/deaugmenting process that was finished
   * for the cell description.
   * Returns false otherwise.
   *
   * \note We use this at the moment only
   * for refinement events. We can consider later
   * on to merge the time stepping functionality
   * (solution update, predictor comp.) into
   * this hook.
   */
  virtual bool leaveCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) = 0;

  /**
   * Update the solution of a cell description.
   */
  virtual void updateSolution(
      const int cellDescriptionsIndex,
      const int element) = 0;

  /**
   * Receive solver data from neighbour rank and write
   * it on the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array at address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * \see tryGetElement
   */
  virtual void mergeNeighbours(
        const int                                 cellDescriptionsIndex1,
        const int                                 element1,
        const int                                 cellDescriptionsIndex2,
        const int                                 element2,
        const tarch::la::Vector<DIMENSIONS, int>& pos1,
        const tarch::la::Vector<DIMENSIONS, int>& pos2) = 0;

  /**
   * Take the cell descriptions \p element
   * from array at address \p cellDescriptionsIndex
   * and merge it with boundary data.
   *
   * \param[in] element Index of the cell description
   *                    at address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * \see tryGetElement
   */
  virtual void mergeWithBoundaryData(
        const int                                 cellDescriptionsIndex,
        const int                                 element,
        const tarch::la::Vector<DIMENSIONS, int>& posCell,
        const tarch::la::Vector<DIMENSIONS, int>& posBoundary) =0;

  #ifdef Parallel
  /**
   * Send solver copy to remote node
   *
   * <h2>Time step restriction</h2>
   *
   * The restrictions of global solver data is not done through the solver
   * objects directly via MPI. See GlobalTimeStepComputation for example.
   *
   * @deprecated
   */
  virtual void sendToRank(int rank, int tag) = 0;

  /**
   * Receive solver copy from remote node
   *
   * @deprecated
   */
  virtual void receiveFromMasterRank(int rank, int tag) = 0;

  /**
   * Send solver data to neighbour rank. Read the data from
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the ADERDGCellDescription
   *                    holding the data to send out in
   *                    the heap vector at \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * \see tryGetElement
   */
  virtual void sendDataToNeighbour(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send empty solver data to neighbour rank.
   */
  virtual void sendEmptyDataToNeighbour(
      const int                                    toRank,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Merge cell description \p element in
   * the cell descriptions array stored at \p
   * cellDescriptionsIndex with metadata.
   *
   * Currently, the neighbour metadata is only the neighbour
   * type as int \p neighbourTypeAsInt.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithNeighbourMetadata(
      const int neighbourTypeAsInt,
      const int cellDescriptionsIndex,
      const int element) = 0;

  /**
   * Receive solver data from neighbour rank and write
   * it on the cell description \p element in
   * the cell descriptions array stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithNeighbourData(
      const int                                    fromRank,
      const int                                    neighbourTypeAsInt,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Drop solver data from neighbour rank.
   */
  virtual void dropNeighbourData(
      const int                                    fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send solver data to master or worker rank. Read the data from
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void sendDataToWorkerOrMasterDueToForkOrJoin(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send empty solver data to master or worker rank
   * due to fork or join.
   */
  virtual void sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
      const int                                    toRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Merge with solver data from master or worker rank
   * that was sent out due to a fork or join. Wrote the data to
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithWorkerOrMasterDataDueToForkOrJoin(
      const int                                    fromRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Drop solver data from master or worker rank
   * that was sent out due to a fork or join.
   */
  virtual void dropWorkerOrMasterDataDueToForkOrJoin(
      const int                                    fromRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  ///////////////////////////////////
  // WORKER->MASTER
  ///////////////////////////////////

  /**
   * Send solver data to master rank. Read the data from
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the heap vector at \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void sendDataToMaster(
      const int                                    masterRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send empty solver data to master rank.
   */
  virtual void sendEmptyDataToMaster(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Merge with solver data from worker rank.
   * Write the data to the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithWorkerData(
      const int                                     workerRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) = 0;

  /**
   * Drop solver data from worker rank.
   */
  virtual void dropWorkerData(
      const int                                     workerRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) = 0;

  ///////////////////////////////////
  // MASTER->WORKER
  ///////////////////////////////////

  /**
   * Send solver data to worker rank. Read the data from
   * the cell description \p element in the cell descriptions
   * vector stored at \p cellDescriptionsIndex.
   *
   * \param[in] element Index of the ADERDGCellDescription
   *                    holding the data to send out in
   *                    the heap vector at \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void sendDataToWorker(
      const int                                     workerRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) = 0;

  /**
   * Send empty solver data to worker rank.
   */
  virtual void sendEmptyDataToWorker(
      const int                                     workerRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) = 0;

  /**
   * Merge with solver data from master rank
   * that was sent out due to a fork or join. Write the data to
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithMasterData(
      const int                                     masterRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) = 0;

  /**
   * Drop solver data from master rank.
   */
  virtual void dropMasterData(
      const int                                     masterRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) = 0;
  #endif
};

#endif
