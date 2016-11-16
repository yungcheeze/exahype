/*
 * LimitingADERDGSolver.h
 *
 *  Created on: 26 Oct 2016
 *      Author: dominic
 */

#ifndef LIMITEDADERDGSOLVER_H_
#define LIMITEDADERDGSOLVER_H_

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "tarch/multicore/BooleanSemaphore.h"

#include "exahype/plotters/Plotter.h"

#include "exahype/mappings/Merging.h"
#include "exahype/mappings/Prediction.h"
#include "exahype/mappings/LimiterStatusSpreading.h"
#include "exahype/mappings/SolutionRecomputation.h"

#include "exahype/plotters/Plotter.h"

#include "exahype/profilers/simple/NoOpProfiler.h"

namespace exahype {
namespace solvers {

class LimitingADERDGSolver;

} /* namespace solvers */
} /* namespace exahype */

class exahype::solvers::LimitingADERDGSolver : public exahype::solvers::Solver {

  /**
   * These mappings need access to the _solver and _limiter fields of LimitingADERDGSolver.
   */
  friend class exahype::mappings::Merging;
  friend class exahype::mappings::Prediction;
  friend class exahype::mappings::LimiterStatusSpreading;
  friend class exahype::mappings::SolutionRecomputation;
private:
  /**
   * A flag indicating that the limiter domain has changed.
   * This might be the case if either a cell has been
   * newly marked as troubled or healed.
   */
  bool _limiterDomainChanged;

  typedef exahype::records::ADERDGCellDescription SolverPatch;
  typedef peano::heap::PlainHeap<SolverPatch> SolverHeap;

  typedef exahype::records::FiniteVolumesCellDescription LimiterPatch;
  typedef peano::heap::PlainHeap<LimiterPatch> LimiterHeap;


  /**
   * Log device.
   */
  static tarch::logging::Log _log;

  /**
   * The ADERDG solver.
   */
  std::unique_ptr<exahype::solvers::ADERDGSolver> _solver;

  /**
   * The finite volumes solver used for the a posteriori subcell limiting.
   */
  std::unique_ptr<exahype::solvers::FiniteVolumesSolver> _limiter;

  /**
   * TODO(Dominc): Remove after docu is recycled.
   *
   * This operation sets the solutions' minimum and maximum value on a cell.
   * The routine is to be invoked after the code has determined the new minimum
   * and maximum value within a cell. In turn, it evaluates whether the new
   * minimum and maximum value have decreased or grown, respectively.
   *
   * If the new min/max values indicate that the new solution comprises
   * oscillations, the routine returns false. This is an indicator that the
   * solution should be limited.
   *
   * If the new min/max values fit, the routine returns true.
   *
   * <h2>Implementation</h2>
   * We hold the min/max information exclusively on the faces. The first thing
   * the routine does is to project the min/max values into the cell. For this
   * it evaluates the 2d faces. The projected value then is compared to the
   * arguments. Once the results of the operation is determined, the routine
   * writes the new arguments onto the 2d face entries. This, on the one hand,
   * stores the data for the subsequent time step, but it also propagates the
   * min/max information into the face-connected neighbours.
   *
   * @param  min          New minimum values within the cell. Array of length
   *                      _numberOfUnknowns.
   * @param  max          New maximum values within the cell
   * @param  solverIndex  Number of the solver within the cell. Please ensure
   *                      that solverIndex refers to an ADER-DG solver.
   * @return True if the new min and max values fit into the restricted min
   *   max solutions. Return false if we seem to run into oscillations.
   */
  //  void setSolutionMinMax(double* min, double* max) const;

  /**
   * Merge the solution min and max values on a face between two cell
   * descriptions. Signature is similar to that of the solver of a Riemann problem.
   */
  void mergeSolutionMinMaxOnFace(
      SolverPatch& pLeft,
      SolverPatch& pRight,
      const int faceIndexLeft,
      const int faceIndexRight) const;

  /**
   * Determine a new merged limiter status based on the neighbours merged
   * limiter status.
   *
   * This method is used during the limiter status spreading.
   * It assumes that the merged limiter statuses are initialised with
   * Ok or Troubled before the spreading.
   * It further assumes that a unique value is written again to
   * the merged limiter status fields after the limiter status of all neighbours
   * of a solver patch have been merged with the limiter status of the patch.
   * (see exahype::solvers::LimitingADERDG::updateLimiterStatus).
   *
   * Determining the merged limiter status:
   * | Status     | Neighbour's Status | Merged Status |
   * --------------------------------------------------|
   * | O          | O/NNT              | O
   * | ~          | T                  | NT
   * | ~          | NT                 | NNT
   * | NNT        | T                  | NT
   */
  void mergeWithNeighbourLimiterStatus(
      SolverPatch& solverPatch,
      const int faceIndex,
      const SolverPatch::LimiterStatus& neighbourLimiterStatus) const;

  void mergeWithNeighbourLimiterStatus(
      SolverPatch& solverPatch,
      const int faceIndex,
      const SolverPatch::LimiterStatus& neighbourLimiterStatus,
      const SolverPatch::LimiterStatus& neighbourOfNeighbourLimiterStatus) const;

  /**
   * Single-sided variant of mergeSolutionMinMaxOnFace() that is required
   * for MPI where min and max value are explicitly exchanged through messages.
   */
  void mergeSolutionMinMaxOnFace(
      SolverPatch&  cellDescription,
      int               faceIndex,
      double* min, double* max) const;

  /**
   * Determine the limiter status after a limiter status spreading
   * iteration.
   */
  exahype::solvers::LimitingADERDGSolver::SolverPatch::LimiterStatus determineLimiterStatus(SolverPatch& solverPatch) const;

  /**
   * Checks if updated solution
   * of the ADER-DG solver is valid
   * or if it contains unphysical oscillations.
   */
  bool solutionIsTroubled(SolverPatch& solverPatch);

  /**
   * Computes the cell-local minimum and maximum
   * of the solution (per variable)
   * and makes them accessible per face.
   */
  void determineMinAndMax(SolverPatch& solverPatch);
public:

  LimitingADERDGSolver(
      std::unique_ptr<exahype::solvers::ADERDGSolver> solver,
      std::unique_ptr<exahype::solvers::FiniteVolumesSolver> limiter,
      std::unique_ptr<profilers::Profiler> profilerSolver =
          std::unique_ptr<profilers::Profiler>(
              new profilers::simple::NoOpProfiler("")),
      std::unique_ptr<profilers::Profiler> profilerLimiter =
          std::unique_ptr<profilers::Profiler>(
              new profilers::simple::NoOpProfiler("")));

  virtual ~LimitingADERDGSolver() {
    _solver.reset();
    _limiter.reset();
  }

  // Disallow copy and assignment
  LimitingADERDGSolver(const ADERDGSolver& other) = delete;
  LimitingADERDGSolver& operator=(const ADERDGSolver& other) = delete;

  /*
   * A time stamp minimised over all the ADERDG and FV solver
   * patches.
   */
  double getMinTimeStamp() const override;

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  double getMinTimeStepSize() const override;

  double getMinNextTimeStepSize() const override;

  void updateMinNextTimeStepSize(double value) override;

  void initInitialTimeStamp(double value) override;

  void synchroniseTimeStepping(
          const int cellDescriptionsIndex,
          const int element) override;

  void startNewTimeStep() override;

  /**
   * Roll back the time step data to the
   * ones of the previous time step.
   */
  void rollbackToPreviousTimeStep();

  void reinitTimeStepData() override;

  bool isValidCellDescriptionIndex(
      const int cellDescriptionsIndex) const override ;

  int tryGetElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const override;

  int tryGetLimiterElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const;

  ///////////////////////////////////
  // MODIFY CELL DESCRIPTION
  ///////////////////////////////////
  bool enterCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) override;

  bool leaveCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) override;

  ///////////////////////////////////
  // CELL-LOCAL
  //////////////////////////////////
  double startNewTimeStep(
      const int cellDescriptionsIndex,
      const int element,
      double*   tempEigenvalues) override;

 /**
   * Rollback to the previous time step, i.e,
   * overwrite the time step size and time stamp
   * fields of the solver and limiter patches
   * by the values used in the previous iteration.
   */
  void rollbackToPreviousTimeStep(
      const int cellDescriptionsIndex,
      const int solverElement);

  /**
   * TODO(Dominic): I need the whole limiter recomputation
   * procedure also for the initial conditions.
   *
   * This includes computing, sending, and merging
   * of the min/max values.
   */
  void setInitialConditions(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override;

  /**
   * Adds a limiter patch to the initially troubled cells, their
   * neighbours, and their neighbour's neighbours.
   * Further imposes finite volumes boundary conditions onto the limiter
   * solution in cells with limiter status Troubled and NeighbourIsTroubled and
   * then projects this solution onto the DG solver space for those cells.
   * Finally, projects the ADER-DG initial conditions onto
   * the FV limiter space for cells with limiter status
   * NeighbourIsNeighbourOfTroubledCell.
   */
  void initialiseLimiter(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const;

  /**
   * This method assumes the ADERDG solver's cell-local limiter status has
   * already been determined.
   *
   * \see determineLimiterStatusAfterLimiterStatusSpreading(...)
   */
  void updateSolution(
      const int cellDescriptionsIndex,
      const int element,
      double** tempStateSizedArrays,
      double** tempUnknowns,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override;

  /**
   * We use this method to determine the limiter status of
   * a cell description after a solution update. This methods checks
   * which cell is troubled and which cell holds a valid solution.
   * If the limiter subdomain changes, i.e., if a cell changes from holding a
   * valid solution (Ok) to troubled (Troubled) or vice versa,
   * this function returns true.
   *
   * The function further sets the merged limiter statuses per face
   * to the value Ok or Troubled.
   * It does not change the value of the cell-based limiter status field.
   * This value is necessary to keep track of the local history of the
   * limiter status.
   */
  bool determineLimiterStatusAfterSolutionUpdate(
      const int cellDescriptionsIndex,
      const int element);

  /**
   * Update the limiter status with a unified value based on
   * the merged limiter statuses per face of the cell.
   *
   * After the status has been determined, it
   * is written to the merged limiter status fields per face.
   *
   * <h2>Determining the unified value</h2>
   * If all of the merged limiter status fields
   * are set to Troubled, the limiter status is changed to Troubled.
   * (There is either all or none of the statuses set to Troubled.)
   *
   * Otherwise and if at least one of the merged statuses is set to NeighbourOfTroubledCell,
   * the status is set to NeighbourOfTroubledCell.
   *
   * Otherwise and if at least one of the merged statuses is set to NeighbourIsNeighbourOfTroubledCell,
   * the status is set to NeighbourIsNeighbourOfTroubledCell.
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   */
  void updateLimiterStatus(int cellDescriptionsIndex, int element);

  /**
   * Reinitialises cells that have been subject to a limiter status change.
   * This method is invoked (during and??) after the limiter status spreading.
   *
   * The method has to take into account which solution, the solver's
   * or the limiter's, was populated with valid solution values
   * in the last iteration. The action of this method is
   * thus based on the new and old limiter status.
   *
   * We perform the following actions based on the
   * old and new limiter status:
   *
   * | New Status | Old Status | Action                                                                                                                                        |
   * ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
   * | O          | O          | Do nothing.                                                                                                                                   |
   * | ~          | T/NT/NNT   | Remove the limiter patch.                                                                                                                     |
   * | T/NT/NNT   | T/NT       | Roll back the limiter solution.                                                                                                               |
   * | ~          | O/NNT      | Roll back the solver solution. Initialise limiter patch if necessary. Project (old, valid) solver solution onto the limiter's solution space. |
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   *
   * We do not overwrite the old limiter status set in this method.
   * We compute the new limiter status based on the merged limiter statuses associated
   * with the faces.
   *
   * TODO(Dominic)
   * Adapters:
   * LimitingADERDGSolver LimiterStatusSpreading
   * LimitingADERDGSolver Reinitialisation
   * LimitingADERDGSolver Recomputation
   */
  void reinitialiseSolvers(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const;

  /**
   * Recompute the solution in cells that have been subject to a limiter status change
   * This method is invoked after the solver reinitialisation
   * (see exahype::solvers::LimitingADERDGSolver::reinitialiseSolvers).
   *
   * It evolves the solution of the solver and limiter in the reinitialised cells to the
   * correct time stamp.
   *
   * We perform the following actions based on the
   * new limiter status:
   *
   * |New Status | Action                                                                                                                                      |
   * ----------------------------------------------------------------------------------------------------------------------------------------------------------|
   * |O          | Do nothing. Solver solution has been evolved correctly before.                                                                              |
   * |T/NT       | Evolve FV solver project result onto the ADER-DG space.                                                                                     |
   * |NNT        | Evolve solver and project its solution onto the limiter solution space. (We had to do a rollback beforehand in the reinitialisation phase.) |
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   *
   * We do not overwrite the old limiter status set in this method.
   * We compute the new limiter status based on the merged limiter statuses associated
   * with the faces.
   *
   * <h2>Overlapping status spreading and reinitialisation with solution reconputation</h2>
   * We can recompute the new solution in cells with status Troubled after one iteration
   * since old solution values from direct neighbours are available then.
   *
   * We can recompute the
   *
   * TODO(Dominic)
   * Adapters:
   * LimitingADERDGSolver LimiterStatusSpreading
   * LimitingADERDGSolver Reinitialisation
   * LimitingADERDGSolver Recomputation
   */
  void recomputeSolution(
      const int cellDescriptionsIndex,
      const int element,
      double** tempStateSizedArrays,
      double** tempUnknowns,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator);

  void preProcess(
      const int cellDescriptionsIndex,
      const int element) override;

  void postProcess(
      const int cellDescriptionsIndex,
      const int element) override;

  void prolongateDataAndPrepareDataRestriction(
      const int cellDescriptionsIndex,
      const int element) override;

  void restrictData(
      const int cellDescriptionsIndex,
      const int element,
      const int parentCellDescriptionsIndex,
      const int parentElement,
      const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) override;

  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  /**
   *
   */
  void mergeLimiterStatusOfNeighbours(
        const int                                 cellDescriptionsIndex1,
        const int                                 element1,
        const int                                 cellDescriptionsIndex2,
        const int                                 element2,
        const tarch::la::Vector<DIMENSIONS, int>& pos1,
        const tarch::la::Vector<DIMENSIONS, int>& pos2);

  /**
   * TODO(Dominic): Not sure if this is necessary.
   */
  void mergeLimiterStatusOfNeighboursOfNeighbours(
          const int                                 cellDescriptionsIndex1,
          const int                                 element1,
          const int                                 cellDescriptionsIndex2,
          const int                                 element2,
          const tarch::la::Vector<DIMENSIONS, int>& pos1,
          const tarch::la::Vector<DIMENSIONS, int>& pos2);

  /**
   * Merge solver boundary data (and other values) of two adjacent
   * cells.
   *
   * The solver involved in the neighbour merge
   * is selected according to the following scheme:
   *
   * | Status A | Status B | Solver to Merge
   * ---------------------------------------
   * | O        | O        | ADER-DG       |
   * | O        | NNT      | ADER-DG       |// O|NNT x O|NNT
   * | NNT      | O        | ADER-DG       |
   * | NNT      | NNT      | ADER-DG       |
   *
   * | NNT      | NT       | FV            |
   * | NT       | NNT      | FV            | // NT&NNT | N&NNT
   *
   * | NT       | NT       | FV            |
   * | NT       | T        | FV            |
   * | T        | NT       | FV            | // T|NT x T|NT
   * | T        | T        | FV            |
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   */
  void mergeNeighbours(
      const int                                 cellDescriptionsIndex1,
      const int                                 element1,
      const int                                 cellDescriptionsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2,
      double**                                  tempFaceUnknownsArrays,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) override;

  /**
   * Merges only the min max of two neighbours.
   *
   * This method is used to detect cells that are
   * troubled after the imposition of initial conditions.
   */
  void mergeSolutionMinMax(
      const int                                 cellDescriptionsIndex1,
      const int                                 element1,
      const int                                 cellDescriptionsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2,
      double**                                  tempFaceUnknownsArrays,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices);

  /**
   * Depending on the limiter status, we impose boundary conditions
   * onto the solution of the solver or of the limiter.
   */
  void mergeWithBoundaryData(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      double**                                  tempFaceUnknownsArrays,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) override;


#ifdef Parallel
  static void sendCellDescriptions(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const peano::heap::MessageType&               messageType,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) {
    assertionMsg(false, "Please implement!");
  }

  static void sendEmptyCellDescriptions(
      const int                                     toRank,
      const peano::heap::MessageType&               messageType,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) {
    assertionMsg(false, "Please implement!");
  }

  /**
   * Receives cell descriptions from rank \p fromRank
   * and resets the data heap indices to -1.
   *
   * If a received cell description has the same
   * solver number as a cell description in the
   * array at address \p cellDescriptionsIndex,
   * we merge the metadata (time stamps, time step size)
   * of both cell descriptions.
   *
   * If no cell description in the array at address
   * \p cellDescriptionsIndex can be found with the
   * same solver number than a received cell description,
   * we push the received cell description to
   * the back of the array at address \p cellDescriptions
   * Index.
   *
   * This operation is intended to be used in combination
   * with the solver method mergeWithWorkerOrMasterDataDueToForkOrJoin(...).
   * Here, we would merge first the cell descriptions sent by the master and worker
   * and then merge the data that is sent out right after.
   */
  static void mergeCellDescriptionsWithRemoteData(
      const int                                     fromRank,
      const int                                     cellDescriptionsIndex,
      const peano::heap::MessageType&               messageType,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) {
    assertionMsg(false, "Please implement!");
  }

  /**
   * Drop cell descriptions received from \p fromRank.
   */
  static void dropCellDescriptions(
      const int                                     fromRank,
      const peano::heap::MessageType&               messageType,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) {
    assertionMsg(false, "Please implement!");
  }

  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  void mergeWithNeighbourMetadata(
        const int neighbourTypeAsInt,
        const int cellDescriptionsIndex,
        const int element) override {
    assertionMsg(false,"Please implement!");
  }

  void sendDataToNeighbour(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     elementIndex,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false,"Please implement!");
  }

  void sendEmptyDataToNeighbour(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }

  void mergeWithNeighbourData(
      const int                                    fromRank,
      const int                                    neighbourTypeAsInt,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      double**                                     tempFaceUnknownsArrays,
      double**                                     tempStateSizedVectors,
      double**                                     tempStateSizedSquareMatrices,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override {
    assertionMsg(false,"Please implement!");
  }


  void dropNeighbourData(
      const int                                     fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }

  void sendDataToWorkerOrMasterDueToForkOrJoin(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }

  void sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }

  void mergeWithWorkerOrMasterDataDueToForkOrJoin(
      const int                                     fromRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }

  void dropWorkerOrMasterDataDueToForkOrJoin(
      const int                                     fromRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }

  ///////////////////////////////////
  // WORKER->MASTER
  ///////////////////////////////////
  bool hasToSendDataToMaster(
        const int cellDescriptionsIndex,
        const int element) override {
    assertionMsg(false, "Please implement!");
  }

  void sendDataToMaster(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override {
    assertionMsg(false, "Please implement!");
  }

  void mergeWithWorkerData(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override {
    assertionMsg(false, "Please implement!");
  }

  void sendDataToMaster(
      const int                                     masterRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }

  void sendEmptyDataToMaster(
      const int                                     masterRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }

  void mergeWithWorkerData(
      const int                                    workerRank,
      const int                                    workerTypeAsInt,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override {
    assertionMsg(false, "Please implement!");
  }

  void dropWorkerData(
      const int                                     workerRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }

  ///////////////////////////////////
  // MASTER->WORKER
  ///////////////////////////////////
  void sendDataToWorker(
      const                                        int workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override {
    assertionMsg(false, "Please implement!");
  }

  void mergeWithMasterData(
      const                                        int masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override {
    assertionMsg(false, "Please implement!");
  }

  void sendDataToWorker(
      const int                                     workerRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }

  void sendEmptyDataToWorker(
      const int                                     workerRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }

  void mergeWithMasterData(
      const int                                     masterRank,
      const int                                     masterTypeAsInt,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }

  void dropMasterData(
      const int                                     masterRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override {
    assertionMsg(false, "Please implement!");
  }
#endif

  std::string toString() const override {
    assertionMsg(false, "Please implement!");
  }

  void toString (std::ostream& out) const override {
    assertionMsg(false, "Please implement!");
  }
};


#endif /* LIMITEDADERDGSOLVER_H_ */
