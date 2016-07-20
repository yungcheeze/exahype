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
#include <vector>

#include "peano/utils/Globals.h"

#include "tarch/la/Vector.h"

#include "exahype/profilers/Profiler.h"
#include "exahype/profilers/simple/NoOpProfiler.h"
#include "exahype/records/ADERDGCellDescription.h"

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

constexpr int addPadding(const int originalSize) {
  return ALIGNMENT/8 * static_cast<int>((originalSize+(ALIGNMENT/8-1))/(ALIGNMENT/8));
}

namespace exahype {
namespace solvers {
class Solver;

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
   * The type of this solver.
   */
  enum class Type { ADER_DG };

  /**
   * The time stepping modus.
   */
  enum class TimeStepping {
    /**
     * In the global time stepping mode, every cells works with the same time step.
     */
    Global,
    // Local, Anarchic
  };

  /**
   * The refinement controls a solver can use.
   */
  enum class RefinementControl { Keep = 0, Refine = 1, Erase = 2 };

 protected:
  /**
   * Each solver has an identifier/name. It is used for debug purposes only.
   */
  const std::string _identifier;

  const Type _type;

  /**
   * Each solver has a kernel number that says which kernel is to be
   * executed. Typically this is an ascending index starting from 0.
   */
  const int _kernelNumber;

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
   * The number of unknowns/basis functions associated with each face of an element.
   * This number includes the unknowns of all state variables.
   */
  const int _unknownsPerFace;

  /**
   * The total number of unknowns/basis functions associated with the 2^d faces of an element.
   * This number includes the unknowns of all state variables.
   */
  const int _unknownsPerCellBoundary;

  /**
   * The total number of unknowns/basis functions associated with the volume of a cell.
   * This number includes the unknowns of all state variables.
   */
  const int _unknownsPerCell;

  /**
   * The total number of volume flux unknowns/basis functions associated with the volume of a cell.
   * This number includes the unknowns of all state variables.
   */
  const int _fluxUnknownsPerCell;

  /**
   * The total number of space-time unknowns/basis functions associated with the
   * space-time volume of a cell and its time stepping interval.
   * This number includes the unknowns of all state variables.
   */
  const int _spaceTimeUnknownsPerCell;

  /**
   * The total number of space-time volume flux unknowns/basis functions associated with the
   * space-time volume of a cell and its time stepping interval.
   * This number includes the unknowns of all state variables.
   */
  const int _spaceTimeFluxUnknownsPerCell;

  /**
    * The size of data required to store cell volume based unknowns and associated parameters.
    */
  const int _dataPerCell;

  /**
   * The time stepping mode of this solver.
   */
  const TimeStepping _timeStepping;

  /**
   * A profiler for this solver.
   */
  std::unique_ptr<profilers::Profiler> _profiler;

  /**
   * Minimum corrector time stamp.
   */
  double _minCorrectorTimeStamp;

  /**
   * Minimum predictor time stamp. Always equal or larger
   * than the minimum corrector time stamp.
   */
  double _minPredictorTimeStamp;

  /**
   * Corrector time step size.
   */
  double _minCorrectorTimeStepSize;

  /**
   * Predictor time step size.
   */
  double _minPredictorTimeStepSize;

  /**
   * Predictor time step size.
   */
  double _minNextPredictorTimeStepSize;

 public:
  Solver(const std::string& identifier, exahype::solvers::Solver::Type type, int kernelNumber,
         int numberOfVariables, int numberOfParameters, int nodesPerCoordinateAxis,
         double maximumMeshSize,
         exahype::solvers::Solver::TimeStepping timeStepping,
         std::unique_ptr<profilers::Profiler> profiler =
             std::unique_ptr<profilers::Profiler>(
                 new profilers::simple::NoOpProfiler));

  virtual ~Solver() {
    // TODO(guera): Remove?!
    _profiler->writeToCout();
  }

  // Disallow copy and assignment
  Solver(const Solver& other) = delete;
  Solver& operator=(const Solver& other) = delete;

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
   * Returns the time stepping mode of this solver.
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
   * This operation returns the number of space time
   * unknowns per cell.
   *
   * Note that this operation might only have a meaning for space-time type
   * discretisation methods.
   */
  int getSpaceTimeUnknownsPerCell() const;

  /**
   * This operation returns the number of space time
   * flux unknowns per cell.
   *
   * Note that this operation might only have a meaning for space-time type
   * discretisation methods.
   */
  int getSpaceTimeFluxUnknownsPerCell() const;

  /**
   * This operation returns the number of unknowns per cell located in
   * the interior of a cell.
   */
  int getUnknownsPerCell() const;

  /**
   * This operation returns the number of flux unknowns per cell
   * located in the interior of a cell.
   */
  int getFluxUnknownsPerCell() const;

  /**
   * This operation returns the number of unknowns that are located
   * on or in the vicinity of the boundary of a cell.
   */
  int getUnknownsPerCellBoundary() const;

  /**
   * This operation returns the number of unknowns that are located
   * on or in the vicinity of each face of a cell.
   */
  int getUnknownsPerFace() const;

  /**
   * If you use a higher order method, then this operation returns the
   * polynomial degree plus one. If you use a Finite Volume method, it
   * returns the number of cells within a patch per coordinate axis.
   */
  int getNodesPerCoordinateAxis() const;

  /**
   * This operation returns the size of data required
   * to store cell volume based unknowns and associated parameters.
   */
  int getDataPerCell() const;

  /**
   * @brief Adds the solution update to the solution.
   *
   * @param[inout] luh  Cell-local solution DoF.
   * @param[in]    lduh Cell-local update DoF.
   * @param[dt]    dt   Time step size.
   */
  virtual void solutionUpdate(double* luh, const double* const lduh,
                              const double dt) = 0;

  /**
   * @brief Computes the volume flux contribution to the cell update.
   *
   * @param[inout] lduh Cell-local update DoF.
   * @param[in]    dx   Extent of the cell in each coordinate direction.
   * @param[dt]    dt   Time step size.
   */
  virtual void volumeIntegral(
      double* lduh, const double* const lFhi,
      const tarch::la::Vector<DIMENSIONS, double>& dx) = 0;

  /**
   * @brief Computes the surface integral contributions
   * to the cell update.
   *
   * @param[inout] lduh   Cell-local update DoF.
   * @param[in]    lFhbnd Cell-local DoF of the boundary extrapolated fluxes.
   * @param[in]    dx     Extent of the cell in each coordinate direction.
   */
  virtual void surfaceIntegral(
      double* lduh, const double* const lFhbnd,
      const tarch::la::Vector<DIMENSIONS, double>& dx) = 0;

  /**
   * @brief Computes the normal fluxes (or fluctuations) at the interface of two
   *cells.
   *
   * @param[inout] FL             Flux DoF belonging to the left cell.
   * @param[inout] FR             Flux DoF belonging the right cell.
   * @param[in]    QL             DoF of the boundary extrapolated predictor
   *belonging to the left cell.
   * @param[in]    QR             DoF of the boundary extrapolated predictor
   *belonging to the right cell.
   * @param[in]    normalNonZero  Index of the nonzero normal vector component,
   *i.e., 0 for e_x, 1 for e_y, and 2 for e_z.
   */
  virtual void riemannSolver(double* FL, double* FR, const double* const QL,
                             const double* const QR, const double dt,
                             const int normalNonZero) = 0;

  /**
   * @brief Computes cell-local predictor space-time, volume, and face DoF.
   *
   * Computes the cell-local space-time predictor lQi, the space-time volume
   *flux lFi,
   * the predictor lQhi, the volume flux lFhi, the boundary
   * extrapolated predictor lQhbnd and normal flux lFhbnd.
   *
   * @param[inout] lQi    Space-time predictor DoF.
   * @param[inout] lFi    Space-time flux DoF.
   * @param[inout] lQhi   Predictor DoF
   * @param[inout] lFhi   Volume flux DoF.
   * @param[inout] lQhbnd Boundary extrapolated predictor DoF.
   * @param[inout] lFhbnd Boundary extrapolated normal fluxes (or fluctuations).
   * @param[out]   luh    Solution DoF.
   * @param[in]    dx     Extent of the cell in each coordinate direction.
   * @param[in]    dt     Time step size.
   */
  virtual void spaceTimePredictor(
      double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd,
      double* lFhbnd, const double* const luh,
      const tarch::la::Vector<DIMENSIONS, double>& dx, const double dt) = 0;

  /**
   * @brief Returns a stable time step size.
   *
   * @param[in] luh       Cell-local solution DoF.
   * @param[in] dx        Extent of the cell in each coordinate direction.
   */
  virtual double stableTimeStepSize(
      const double* const luh,
      const tarch::la::Vector<DIMENSIONS, double>& dx) = 0;

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
      const int level) = 0;  // @todo make abstract

  /**
   * Project coarse grid face unknowns
   * on level \p coarseGridLevel down to level \p fineGridLevel
   * and writes them to the fine grid unknowns
   *
   * \note For the considered AMR concept, the difference in levels can
   * be larger than one. Let \f$l\f$ be the level difference. The
   * vector \p subfaceIndex does contain values in the range
   * \f$0,1,\ldots,3^l-1\f$.
   */
  virtual void faceUnknownsProlongation(
      double* lQhbndFine, double* lFhbndFine, const double* lQhbndCoarse,
      const double* lFhbndCoarse, const int coarseGridLevel,
      const int fineGridLevel,
      const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex) = 0;

  /**
   * Restricts fine grid face unknowns on level \p fineGridLevel
   * up to level \p coarseGridLevel and adds them to the coarse grid unknowns.
   *
   * \note For the considered AMR concept, the difference in levels can
   * be larger than one. Let \f$l\f$ be the level difference. The
   * vector \p subfaceIndex does contain values in the range
   * \f$0,1,\ldots,3^l-1\f$.
   */
  virtual void faceUnknownsRestriction(
      double* lQhbndCoarse, double* lFhbndCoarse, const double* lQhbndFine,
      const double* lFhbndFine, const int coarseGridLevel,
      const int fineGridLevel,
      const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex) = 0;

  /**
   * Project coarse grid face unknowns
   * on level \p coarseGridLevel down to level \p fineGridLevel
   * and writes them to the fine grid unknowns
   *
   * \note For the considered AMR concept, the difference in levels is always
   * equal to one. The vector \p subcellIndex does contain values in the range
   * \f$0,1,2\f$.
   */
  virtual void volumeUnknownsProlongation(
      double* luhFine, const double* luhCoarse, const int coarseGridLevel,
      const int fineGridLevel,
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) = 0;

  /**
   * Restricts fine grid volume unknowns on level \p fineGridLevel
   * up to level \p coarseGridLevel and adds them to the coarse grid unknowns.
   *
   * \note For the considered AMR concept, the difference in levels is always
   * equal to one. The vector \p subcellIndex does contain values in the range
   * \f$0,1,2\f$.
   */
  virtual void volumeUnknownsRestriction(
      double* luhCoarse, const double* luhFine, const int coarseGridLevel,
      const int fineGridLevel,
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) = 0;
  ///@}

  /**
   * @todo Dominic, please add a description what this routine does.
   */
  void synchroniseTimeStepping(
      exahype::records::ADERDGCellDescription& p) const;

  void startNewTimeStep();

  void updateMinNextPredictorTimeStepSize(
      const double& nextPredictorTimeStepSize);

  double getMinNextPredictorTimeStepSize() const;

  // todo 16/02/25:Dominic Etienne Charrier: It follows stuff that must be
  // revised:

  // todo 25/02/16:Dominic Etienne Charrier
  // Remove the time stamps that are not used in ExaHype.
  void setMinCorrectorTimeStamp(double minCorectorTimeStamp);

  double getMinCorrectorTimeStamp() const;

  void setMinPredictorTimeStamp(double minPredictorTimeStamp);

  double getMinPredictorTimeStamp() const;

  // @todo 25/02/16:Dominic Etienne Charrier
  // @Tobias: The time step size getters are only used for
  // debugging/assertion purposes at the moment and will
  // be removed if the time stepping works robust again.
  double getMinCorrectorTimeStepSize() const;

  void setMinPredictorTimeStepSize(double minPredictorTimeStepSize);

  double getMinPredictorTimeStepSize() const;
  /**
   * Update predictor time step size
   *
   * This operation takes the minimum of the current predictor time step size
   * and the argument handed in. The routine is used in
   * GlobalTimeStepComputation to determine the subsequent time step size.
   *
   * <h1>Thread-safety</h1>
   *
   * This operation is not thread safe.
   *
   */
  void updateNextPredictorTimeStepSize(double nextPredictorTimeStepSize);

  /**
   * This operation has to different branches: one for the master and one for
   * the worker. If we are in the master, we basically do only send out all
   * the solver data to the worker. If we are on the worker, we do overwrite
   * all solver data accordingly.
   */
  void sendToRank(int rank, int tag);
  void receiveFromRank(int rank, int tag);

  void toString();

  /**
   * Run over all solvers and identify the minimal time stamp.
   */
  static double getMinSolverTimeStamp();

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  static double getMinSolverTimeStepSize();
};

#endif
