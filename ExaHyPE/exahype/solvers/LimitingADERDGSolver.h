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

#include "exahype/plotters/Plotter.h"

#include "exahype/profilers/simple/NoOpProfiler.h"

#include "exahype/solvers/TemporaryVariables.h"

namespace exahype {
namespace solvers {

class LimitingADERDGSolver;

} /* namespace solvers */
} /* namespace exahype */

class exahype::solvers::LimitingADERDGSolver : public exahype::solvers::Solver {

private:
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
   * A flag indicating that the limiter domain has changed.
   * This might be the case if either a cell has been
   * newly marked as troubled or healed.
   */
  bool _limiterDomainHasChanged;

  /**
   * The limiterDomainHasChanged for the next
   * iteration.
   */
  bool _nextLimiterDomainHasChanged;

  /**
   * The maximum relaxation parameter
   * used for the discrete maximum principle.
   */
  double _DMPMaximumRelaxationParameter;

  /**
   * The difference scaling
   * used for the discrete maximum principle.
   */
  double _DMPDifferenceScaling;

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
   * Determine a new limiter status for the given direction based on the neighbour's
   * limiter status and the cell's reduced limiter status.
   *
   * Computes the new limiter status \f$ L_\rm{new} \f$ per direction
   * according to:
   *
   * \f[
   *  L_\rm{new} = \begin{cases}
   *  T & L = T \\
   *  \max \left( 0, \max \left( L, L_\rm{neighbour} \right) -1 \right) & \rm{else}
   *   \end{cases}
   * \f]
   *
   * with \f$ L \f$, \f$ L_\rm{neighbour} \f$, denoting the current limiter status
   * of the cell and the neighbour, respectively, and \f$  T  \f$ indicates the status
   * of a troubled cell.
   */
  void mergeWithLimiterStatus(
      SolverPatch& solverPatch,
      const int faceIndex,
      const SolverPatch::LimiterStatus& neighbourLimiterStatus) const;

  void mergeWithNeighbourLimiterStatus(
      SolverPatch& solverPatch,
      const int faceIndex,
      const SolverPatch::LimiterStatus& neighbourLimiterStatus,
      const SolverPatch::LimiterStatus& neighbourOfNeighbourLimiterStatus) const;

  /**
   * Checks if the updated solution
   * of the ADER-DG solver contains unphysical oscillations (false)
   * or not (true).
   *
   * Compute the new min and max at the same time.
   */
  bool evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(SolverPatch& solverPatch);

  /**
   * Checks if the updated solution
   * of the ADER-DG solver contains unphysical oscillations (false)
   * or not (true).
   */
  bool evaluateDiscreteMaximumPrinciple(SolverPatch& solverPatch);

  /**
   * Checks if the updated solution
   * of the ADER-DG solver is
   * a physically admissible one (true).
   */
  bool evaluatePhysicalAdmissibilityCriterion(SolverPatch& solverPatch);

  /**
   * Computes the cell-local minimum and maximum
   * of the solution (per variable)
   * and makes them accessible per face.
   */
  void determineSolverMinAndMax(SolverPatch& solverPatch);

  /**
   * Computes the cell-local minimum and maximum
   * of the limiter's solution (per variable)
   * and makes them accessible per face.
   *
   * This method is used for troubled cells that
   * do not hold a valid ADER-DG solution,
   * as well as their neighbours.
   */
  void determineLimiterMinAndMax(SolverPatch& solverPatch,LimiterPatch& limiterPatch);

  /**
   * Updates the merged limiter status based on the cell-local ADER-DG solution
   * values,
   */
  void updateLimiterStatusAfterSolutionUpdate(SolverPatch& solverPatch,const bool isTroubled);

  void mergeWithLimiterStatusOfNeighbours(
      SolverPatch& solverPatch,
      const tarch::la::Vector<THREE_POWER_D, int>& neighbourCellDescriptionsIndices) const;

#ifdef Parallel
  /**
   * Data messages per neighbour communication.
   * This information is required by the sendEmpty...(...)
   * methods.
   */
  static const int DataMessagesPerNeighbourCommunication;
  /**
   * Data messages per fork/join communication.
   * This information is required by the sendEmpty...(...)
   * methods.
   */
  static const int DataMessagesPerForkOrJoinCommunication;
  /**
   * Data messages per master worker communication.
   * This information is required by the sendEmpty...(...)
   * methods.
   */
  static const int DataMessagesPerMasterWorkerCommunication;

  /**
   * Send the solution minimum and maximum values per variable
   * and further the merged limiter status of the solver patch
   * \p element in heap array \p cellDescriptionsIndex to the
   * neighbour.
   *
   * \see exahype::solvers::Solver::sendDataToNeighbour
   * for a description of the parameters.
   */
  void sendMinAndMaxToNeighbour(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  /**
   * Receive the solution minimum and maximum values per variable
   * and further the merged limiter status for the solver patch
   * \p element in heap array \p cellDescriptionsIndex from the
   * neighbour.
   *
   * \see exahype::solvers::Solver::sendDataToNeighbour
   * for a description of the parameters.
   */
  void mergeWithNeighbourMinAndMax(
      const int                                    fromRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  /**
   * Single-sided variant of mergeSolutionMinMaxOnFace() that is required
   * for MPI where min and max value are explicitly exchanged through messages.
   */
  void mergeSolutionMinMaxOnFace(
      SolverPatch&  cellDescription,
      const int     faceIndex,
      const double* const min, const double* const  max) const;

//  /**
//   * Send the limiter status
//   * of the solver patch \p element in heap array
//   * \p cellDescriptionsIndex to the respective neighbour.
//   *
//   * \see exahype::solvers::Solver::sendDataToNeighbour
//   * for a description of the parameters.
//   *
//   * <h2>Heap</h2>
//   * We currently use the DoubleHeap to
//   * communicate the enum.
//   * While it would make sense to use the
//   * MetadataHeap for this purpose,
//   * there is a possibility of intermixing
//   * the communication of the limiter status
//   * with the communication of the solver
//   * metadata.
//   */
//  void sendLimiterStatusToNeighbour(
//      const int                                     toRank,
//      const int                                     cellDescriptionsIndex,
//      const int                                     element,
//      const tarch::la::Vector<DIMENSIONS, int>&     src,
//      const tarch::la::Vector<DIMENSIONS, int>&     dest,
//      const tarch::la::Vector<DIMENSIONS, double>&  x,
//      const int                                     level) const;
#endif

public:
  static bool irregularChangeOfLimiterDomainOfOneSolver();

  /**
   * Create a limiting ADER-DG solver.
   *
   * <h2>Discrete maximum principle</h2>
   * By default this constructor initialises the maximum relaxation
   * parameter to the value to \f$ \delta_0 = 1\cdot 10^{-4} \f$
   * and the difference scaling parameter to \f$ \epsilon = 1\cdot 10^{-3} \f$.
   * See Dumbser et al., 2014. doi:10.1016/j.jcp.2014.08.009 for more details on
   * the notation.
   */
  LimitingADERDGSolver(
      const std::string& identifier,
      std::unique_ptr<exahype::solvers::ADERDGSolver> solver,
      std::unique_ptr<exahype::solvers::FiniteVolumesSolver> limiter,
      const double DMPRelaxationParameter=1e-4,
      const double DMPDifferenceScaling=1e-3
      );

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

  void initSolver(const double timeStamp, const tarch::la::Vector<DIMENSIONS,double>& boundingBox) override;

  void synchroniseTimeStepping(
          const int cellDescriptionsIndex,
          const int element) override;

  void reconstructStandardTimeSteppingData() {
    _solver->reconstructStandardTimeSteppingData();
    _limiter->setMinTimeStamp(_solver->getMinCorrectorTimeStamp());
    _limiter->setMinTimeStepSize(_solver->getMinCorrectorTimeStepSize());
  }

  /**
   * We always override the limiter time step
   * data by the ADER-DG one before a solution update.
   */
  void startNewTimeStep() override;

  void zeroTimeStepSizes() override;

  bool getLimiterDomainHasChanged() {
    return _limiterDomainHasChanged;
  }

  bool getNextLimiterDomainHasChanged() {
    return _nextLimiterDomainHasChanged;
  }

  void updateNextLimiterDomainHasChanged(bool state) {
    _nextLimiterDomainHasChanged |= state;
  }

  /**
   * Roll back the time step data to the
   * ones of the previous time step.
   */
  void rollbackToPreviousTimeStep();

  void reconstructStandardTimeSteppingDataAfterRollback();

  void reinitialiseTimeStepData() override;

  void updateNextMinCellSize(double minCellSize) override;
  void updateNextMaxCellSize(double maxCellSize) override;
  double getNextMinCellSize() const override;
  double getNextMaxCellSize() const override;
  double getMinCellSize() const override;
  double getMaxCellSize() const override;

  bool isValidCellDescriptionIndex(
      const int cellDescriptionsIndex) const override ;

  /**
   * Returns the index of the solver patch registered for the solver with
   * index \p solverNumber in exahype::solvers::RegisteredSolvers.
   * If no patch is registered for the solver, this function
   * returns exahype::solvers::Solver::NotFound.
   */
  int tryGetElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const override {
    return _solver->tryGetElement(cellDescriptionsIndex,solverNumber);
  }

  /**
   * Returns the index of the limiter patch registered for the solver with
   * index \p solverNumber in exahype::solvers::RegisteredSolvers.
   * If no patch is registered for the solver, this function
   * returns exahype::solvers::Solver::NotFound.
   */
  int tryGetLimiterElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const {
    return _limiter->tryGetElement(cellDescriptionsIndex,solverNumber);
  }

  /**
   * Returns the index of the limiter patch corresponding to
   * the solver patch with index \p solverElement.
   * Both patches link to the same ::LimitingADERDGSolver.
   * If no limiter patch is found, this function
   * returns exahype::solvers::Solver::NotFound.
   */
  int tryGetLimiterElementFromSolverElement(
      const int cellDescriptionsIndex,
      const int solverElement) const {
    SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);
    return _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  }

  /**
    * \see exahype::amr::computeSubcellPositionOfCellOrAncestor
    */
  SubcellPosition computeSubcellPositionOfCellOrAncestor(
      const int cellDescriptionsIndex,
      const int element) override;

  ///////////////////////////////////
  // MODIFY CELL DESCRIPTION
  ///////////////////////////////////
  /**
   * Based on the limiter status of a solver patch
   * and the solver patch's type, we perform the
   * following actions:
   *
   * | New Status | Type                        | Action                                                                                           |
   * ----------------------------------------------------------------------------------------------------------------------------------------------|
   * | O/NNT      | Any                         | Do nothing.                                                                                      |
   * | T/NT       | Cell                        | Set RefinementRequested event on parent cell if its current event is None or AugmentingRequested |
   * | T/NT       | Descendant                  | Set RefinementRequested event if current event is None or AugmentingRequested                    |
   * | T/NT       | Else                        | Do nothing                                                                                       |
   *
   * \note Currently we assume that the problem and load-balancing is so well-behaved that
   * we always find a Cell as parent of a Descendant on the same MPI rank. We further do not
   * consider Master-Worker boundaries in the lookup of the parent.
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   */
  bool markForRefinementBasedOnLimiterStatus(
        exahype::Cell& fineGridCell,
        exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
        const bool initialGrid,
        const int solverNumber);

  void mergeWithLimiterStatusOfNeighbours(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const;

  bool markForRefinement(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const bool initialGrid,
      const int solverNumber) override;

  bool updateStateInEnterCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const bool initialGrid,
      const int solverNumber) override;

  bool updateStateInLeaveCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) override;

  bool attainedStableState(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      const int solverNumber) const override;

  void finaliseStateUpdates(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) override;

  ///////////////////////////////////
  // CELL-LOCAL
  //////////////////////////////////
  bool evaluateRefinementCriterionAfterSolutionUpdate(
      const int cellDescriptionsIndex,
      const int element) override;

  double startNewTimeStep(
      const int cellDescriptionsIndex,
      const int element,
      double*   tempEigenvalues) override;

  void zeroTimeStepSizes(const int cellDescriptionsIndex, const int solverElement) override;

  void reconstructStandardTimeSteppingData(const int cellDescriptionsIndex,int element) const {
    _solver->reconstructStandardTimeSteppingData(cellDescriptionsIndex,element);

    SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
    const int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    if (limiterElement!=exahype::solvers::Solver::NotFound) {
      LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
      limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
    }
  }

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
   * TODO(Dominic): Docu
   */
  void reconstructStandardTimeSteppingDataAfterRollback(
      const int cellDescriptionsIndex,
      const int element) const;

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
   * This method assumes the ADERDG solver's cell-local limiter status has
   * already been determined.
   *
   * Before performing an update with the limiter,
   * set the ADER-DG time step sizes.
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
   * Determine the new cell-local min max values.
   *
   * Must be invoked after ::determineLimiterStatusAfterSolutionUpdate.
   */
  void determineMinAndMax(
      const int cellDescriptionsIndex,
      const int element);

  /**
   * Evaluates a discrete maximum principle (DMP) and
   * the physical admissibility detection (PAD) criterion for
   * the solution values stored for a solver patch.
   *
   * This method then invokes
   * ::determinLimiterStatusAfterSolutionUpdate(SolverPatch&,const bool)
   * with the result of these checks.
   */
  bool updateLimiterStatusAndMinAndMaxAfterSolutionUpdate(
      const int cellDescriptionsIndex,
      const int element);

  /**
   * Similar to ::determineLimiterStatusAfterSolutionUpdate(const int,const int)
   * but does not evaluate the discrete maximum principle.
   */
  bool updateLimiterStatusAndMinAndMaxAfterSetInitialConditions(
      const int cellDescriptionsIndex,
      const int element);

  /**
   * Update the limiter status based on the cell-local solution values.
   * If a cell
   *
   * \param[in] isTroubled A bool indicating if the patch's solution is (still) troubled
   *
   * \return True if the limiter domain changes irregularly in the cell, i.e.,
   * if a patch with status Ok, NeighbourOfTroubled3, NeighbourOfTroubled4
   * changes its status to Troubled.
   *
   * If the limiter status changes regularly, i.e., from NeighbourOfTroubled1
   * to Troubled or from Troubled to NeighbourOfTroubled3, NeighbourOfTroubled4, this
   * methods returns false.
   */
  bool determineLimiterStatusAfterSolutionUpdate(
      SolverPatch& solverPatch,const bool isTroubled) const;

  /**
   * Update the merged limiter status with a unified value based on
   * the merged limiter statuses per face of the cell.
   *
   * \see ::determineLimiterStatus
   */
  void updateLimiterStatus(
      const int cellDescriptionsIndex,const int element) const;

  /**
   * Update the limiter status with a unified value based on
   * the merged limiter statuses per face of the cell.
   *
   * \see ::determineLimiterStatus
   */
  void updatePreviousLimiterStatus(
      const int cellDescriptionsIndex,const int element) const;


  /**
   * Deallocates the limiter patch for solver patches
   * that are either not of type Cell or flagged with limiter status
   * Ok.
   *
   * On the other hand, allocate a limiter patch for solver
   * patches of type Cell that are flagged with a limiter
   * status other than Ok.
   *
   * \return true if a limiter patch was allocated. Return false
   * otherwise.
   */
  bool allocateOrDeallocateLimiterPatch(
      const int cellDescriptionsIndex,
      const int solverNumber) const;

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
   * LimitingADERDGSolver SolutionRecomputation
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
   * |O          | Do nothing. Solver solution has been evolved correctly before. DG solution is correct.                                                                             |
   * |T/NT       | Evolve FV solver. Project result onto the ADER-DG space. Recompute the space-time predictor if not initial recomputation.                                                                                  |
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
   * LimitingADERDGSolver SolutionRecomputation
   */
  void recomputeSolution(
      const int cellDescriptionsIndex,
      const int element,
      exahype::solvers::SolutionUpdateTemporaryVariables& solutionUpdateTemporaryVariables,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator);

  /**
   * !!! Only for fused time stepping !!!
   *
   * Recompute the predictor in particular cells.
   *
   * We perform the following actions based on the
   * new and old limiter status:
   *
   * |New Status | Old Status | Action                                                                                                                                      |
   * ----------------------------------------------------------------------------------------------------------------------------------------------------------|
   * |O/NNTT     | O/NT/NNT   | Do nothing. Underlying DG solution has not changed.
   * |           | T          | Recompute predictor. Cell was skipped before in predictor computation
   * |           |            | since it is marked T. See mapping Prediction::enterCell.
   * |NT         | *          | Recompute predictor. DG solution has been recomputed.
   * |T          | *          | Not necesssary to compute predictor
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   *
   * We do not overwrite the old limiter status set in this method.
   * We compute the new limiter status based on the merged limiter statuses associated
   * with the faces.
   */
  void recomputePredictor(
      const int cellDescriptionsIndex,
        const int element,
        exahype::solvers::PredictionTemporaryVariables& predictionTemporaryVariables,
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

  /**
   * Restrict and merge the merged limiter status of a solver patch of
   * type Cell or Ancestor with the next Ancestor and
   * all the intermediate EmptyAncestors.
   * Perform the merge only if the cell associated with the solver patch is
   * adjacent to the boundary of the cell associated with the intermediate
   * EmptyAncestors or next Ancestor.
   */
  void mergeLimiterStatusWithAncestors(
      const int cellDescriptionsIndex,
      const int element);

  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  void mergeNeighbours(
      const int                                 cellDescriptionsIndex1,
      const int                                 element1,
      const int                                 cellDescriptionsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) override;

  /**
   * Merge solver boundary data (and other values) of two adjacent
   * cells based on their limiter status.
   *
   * The solver involved in the neighbour merge
   * is selected according to the following scheme:
   *
   * | Status 1 | Status 2 | Solver to Merge
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
   *
   * <h2>Solution Recomputation</h2>
   * If we perform a solution recomputation, we do not need to perform
   * a solution update in the non-troubled solver patches and the patches
   * with status neighbourIsNeighbourOfTroubledCell. Instead we
   * can simply reuse the already computed update to advance in
   * time to the desired time stamp.
   *
   * We thus do not need to merge these patches with neighbour data
   * in the recomputation phase.
   *
   * TODO(Dominic): Remove limiterstatus1 and limiterStatus2 argument.
   * They depend on the isRecomputation value
   */
  void mergeNeighboursBasedOnLimiterStatus(
      const int                                 cellDescriptionsIndex1,
      const int                                 element1,
      const int                                 cellDescriptionsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2,
      const bool                                isRecomputation,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices);

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
      double**                                  tempFaceUnknowns,
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
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) override;

  /**
   * Merge solver boundary data (and other values) of a
   * cell with the boundary conditions based on the cell's
   * limiter status.
   *
   * The solver involved in the merge
   * is selected according to the following scheme:
   *
   * | Status   | Solver to Merge |
   * ------------------------------
   * | O        | ADER-DG         |
   * | NNT      | ADER-DG         |
   *
   * | NT       | FV              |
   * | T        | FV              |
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   *
   * <h2>Solution Recomputation</h2>
   * If we perform a solution recomputation, we do not need to perform
   * a solution update in the non-troubled solver patches and the patches
   * with status neighbourIsNeighbourOfTroubledCell. Instead we
   * can simply reuse the already computed update to advance in
   * time to the desired time stamp.
   *
   * We thus do not need to merge these patches with boundary data
   * in the recomputation phase.
   *
   * \param[in] isRecomputation Flag indicating if this merge is part of a solution recomputation phase.
   */
  void mergeWithBoundaryDataBasedOnLimiterStatus(
        const int                                 cellDescriptionsIndex,
        const int                                 element,
        const SolverPatch::LimiterStatus&         limiterStatus,
        const tarch::la::Vector<DIMENSIONS, int>& posCell,
        const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
        const bool                                isRecomputation,
        double**                                  tempFaceUnknowns,
        double**                                  tempStateSizedVectors,
        double**                                  tempStateSizedSquareMatrices);

  void prepareNextNeighbourMerging(
      const int cellDescriptionsIndex,const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const override;

#ifdef Parallel
  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  void mergeWithNeighbourMetadata(
      const exahype::MetadataHeap::HeapEntries& neighbourMetadata,
      const tarch::la::Vector<DIMENSIONS, int>& src,
      const tarch::la::Vector<DIMENSIONS, int>& dest,
      const int                                 cellDescriptionsIndex,
      const int                                 element) override;

  void sendDataToNeighbour(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  /**
   * Send data or empty data to the neighbour data based
   * on the limiter status.
   *
   * \param[in] isRecomputation Indicates if this called within a solution recomputation
   *                            process.
   * \param[in] limiterStatus   This method is used in two different contexts (see \p isRecomputation)
   *                            which either make use of the unified face-wise limiter status (isRecomputation)
   *                            or make use of the cell-wise limiter status (!isRecomputation).
   *
   * <h2>SolutionRecomputation</h2>
   * In case this method is called within a solution recomputation,
   * the ADER-DG solver does only send empty messages to the neighbour.
   * Otherwise it merges the received data and adds it to the update.
   *
   * \note This method assumes that there has been a unified face-wise limiter status value
   * determined and written back to the faces a-priori.
   *
   * <h2>Possible optimisations</h2>
   * Depending on isRecomputation we do not need to send both, solver and limiter,
   * data for patches with status NeighbourIsNeighbourOfTroubledCell and NeighbourOfTroubledCell.
   */
  void sendDataToNeighbourBasedOnLimiterStatus(
        const int                                    toRank,
        const int                                    cellDescriptionsIndex,
        const int                                    element,
        const tarch::la::Vector<DIMENSIONS, int>&    src,
        const tarch::la::Vector<DIMENSIONS, int>&    dest,
        const bool                                   isRecomputation,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level) const;

  void sendEmptyDataToNeighbour(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void mergeWithNeighbourData(
      const int                                    fromRank,
      const exahype::MetadataHeap::HeapEntries&    neighbourMetadata,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      double**                                     tempFaceUnknowns,
      double**                                     tempStateSizedVectors,
      double**                                     tempStateSizedSquareMatrices,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  /**
   * Merge or drop received neighbour data based
   * on the limiter status.
   *
   * \param isRecomputation Indicates if this called within a solution recomputation
   *                        process.
   * \param limiterStatus   This method is used in two different contexts (see \p isRecomputation)
   *                        which either make use of the unified face-wise limiter status (isRecomputation)
   *                        or make use of the cell-wise limiter status (!isRecomputation).
   *
   * <h2>SolutionRecomputation</h2>
   * In case this method is called within a solution recomputation,
   * the ADER-DG solver drops the received boundary data.
   * Otherwise it merges the received data and adds it to the update.
   *
   *  \note This method assumes that there has been a unified face-wise limiter status value
   *  determined and written back to the faces.
   */
  void mergeWithNeighbourDataBasedOnLimiterStatus(
      const int                                    fromRank,
      const exahype::MetadataHeap::HeapEntries&    neighbourMetadata,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const bool                                   isRecomputation,
      double**                                     tempFaceUnknowns,
      double**                                     tempStateSizedVectors,
      double**                                     tempStateSizedSquareMatrices,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  void dropNeighbourData(
      const int                                     fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  ///////////////////////////////////////
  // NEIGHBOUR - Limiter status spreading
  ///////////////////////////////////////

  /**
   * Receive and merge the merged limiter status sent
   * by the neighbour at position \p src - \p dest.
   *
   * see exahype::solvers::Solver::mergeWithNeighbourData
   * for details on the parameters.
   *
   * @deprecated see mergeWithNeighbourMetadata
   */
  void mergeWithNeighbourLimiterStatus(
      const int                                    fromRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  /**
   * Drop the merged limiter status sent
   * by the neighbour at position \p src - \p dest.
   *
   * \see exahype::solvers::Solver::dropNeighbourData
   * for details on the parameters.
   */
  void dropNeighbourLimiterStatus(
      const int                                    fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  /**
   *  Send the merged limiter status at position \p dest- \p src to
   *  the corresponding neighbour.
   *
   * \see exahype::solvers::Solver::sendDataToNeighbour
   * for details on the parameters.
   */
  void sendLimiterStatusToNeighbour(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const;

  /**
   *  Send an empty array instead of
   *  the merged limiter status at position \p dest- \p src to
   *  the corresponding neighbour.
   *
   *  \see exahype::solvers::Solver::sendEmptyDataToNeighbour
   * for details on the parameters.
   */
  void sendEmptyDataInsteadOfLimiterStatusToNeighbour(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const;

  ///////////////////////////////////////
  // NEIGHBOUR - Solution Recomputation
  ///////////////////////////////////////
  /**
   * TODO(Dominic):
   * Add more docu.
   *
   * We do not send the minimum and maximum solution values
   * during the solution recomputation process.
   */
  void sendEmptySolverAndLimiterDataToNeighbour(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const;

  /**
   * TODO(Dominic):
   * Add more docu.
   *
   * We do not send the minimum and maximum solution values
   * during the solution recomputation process.
   */
  void dropNeighbourSolverAndLimiterData(
        const int                                     fromRank,
        const tarch::la::Vector<DIMENSIONS, int>&     src,
        const tarch::la::Vector<DIMENSIONS, int>&     dest,
        const tarch::la::Vector<DIMENSIONS, double>&  x,
        const int                                     level) const;

  /////////////////////////////////////
  // FORK OR JOIN
  /////////////////////////////////////
  void sendDataToWorkerOrMasterDueToForkOrJoin(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void mergeWithWorkerOrMasterDataDueToForkOrJoin(
      const int                                     fromRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void dropWorkerOrMasterDataDueToForkOrJoin(
      const int                                     fromRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  ///////////////////////////////////
  // WORKER->MASTER
  ///////////////////////////////////
  bool hasToSendDataToMaster(
        const int cellDescriptionsIndex,
        const int element) override;

  void sendDataToMaster(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void mergeWithWorkerData(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void sendDataToMaster(
      const int                                     masterRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void sendEmptyDataToMaster(
      const int                                     masterRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void mergeWithWorkerData(
      const int                                    workerRank,
      const exahype::MetadataHeap::HeapEntries&    workerMetadata,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void dropWorkerData(
      const int                                     workerRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  ///////////////////////////////////
  // MASTER->WORKER
  ///////////////////////////////////
  void sendDataToWorker(
      const                                        int workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void mergeWithMasterData(
      const                                        int masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void sendDataToWorker(
      const int                                     workerRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void sendEmptyDataToWorker(
      const int                                     workerRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void mergeWithMasterData(
      const int                                     masterRank,
      const exahype::MetadataHeap::HeapEntries&     masterMetadata,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void dropMasterData(
      const int                                     masterRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;
#endif

  std::string toString() const override;

  void toString (std::ostream& out) const override;

  const std::unique_ptr<exahype::solvers::FiniteVolumesSolver>&
  getLimiter () const {
    return _limiter;
  }

  const std::unique_ptr<exahype::solvers::ADERDGSolver>&
  getSolver () const {
    return _solver;
  }
};


#endif /* LIMITEDADERDGSOLVER_H_ */
