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
  enum class Type { ADER_DG, FiniteVolumes };

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
  Solver(
    const std::string& identifier,
    exahype::solvers::Solver::Type type,
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

  /**
   * Send solver copy to remote node
   */
  virtual void sendToRank(int rank, int tag) = 0;


  /**
   * Receive solver copy from remote node
   */
  virtual void receiveFromMasterRank(int rank, int tag) = 0;
  virtual void receiveFromWorkerRank(int rank, int tag) = 0;

  void toString();

  /**
   * Run over all solvers and identify the minimal time stamp.
   */
  virtual double getMinTimeStamp() const = 0;

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  virtual double getMinTimeStepSize() const = 0;

  virtual void updateNextTimeStepSize( double value ) = 0;

  virtual void initInitialTimeStamp(double value) = 0;

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
};

#endif
