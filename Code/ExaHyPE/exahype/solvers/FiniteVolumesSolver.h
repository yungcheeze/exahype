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

#ifndef _EXAHYPE_SOLVERS_FINITE_VOLUMES_SOLVER_H_
#define _EXAHYPE_SOLVERS_FINITE_VOLUMES_SOLVER_H_


#include "exahype/solvers/Solver.h"

namespace exahype {
  namespace solvers {
    class FiniteVolumesSolver;
  }  // namespace solvers
}  // namespace exahype



class exahype::solvers::FiniteVolumesSolver: public exahype::solvers::Solver {
  private:
    /**
     * Total number of unknowns in a cell.
     */
    int _unknownsPerCell;
  public:
    FiniteVolumesSolver(
      const std::string& identifier,
      int numberOfVariables,
      int numberOfParameters,
      int nodesPerCoordinateAxis,
      double maximumMeshSize,
      exahype::solvers::Solver::TimeStepping timeStepping,
      std::unique_ptr<profilers::Profiler> profiler =
        std::unique_ptr<profilers::Profiler>( new profilers::simple::NoOpProfiler)
    );

    virtual ~FiniteVolumesSolver() {}

    // Disallow copy and assignment
    FiniteVolumesSolver(const FiniteVolumesSolver& other) = delete;
    FiniteVolumesSolver& operator=(const FiniteVolumesSolver& other) = delete;

    virtual double stableTimeStepSize(
        const double* const luh,
        const tarch::la::Vector<DIMENSIONS, double>& dx) override;

    virtual void solutionAdjustment(
        double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,
        const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) override;

    virtual bool hasToAdjustSolution(
        const tarch::la::Vector<DIMENSIONS, double>& center,
        const tarch::la::Vector<DIMENSIONS, double>& dx, double t) override;

    virtual exahype::solvers::Solver::RefinementControl refinementCriterion(
        const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,
        const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
        const int level) override;

    virtual double getMinTimeStamp() const override;

    /**
     * This operation returns the number of unknowns per cell located in
     * the interior of a cell.
     */
    int getUnknownsPerCell() const;

    /**
     * Run over all solvers and identify the minimal time step size.
     */
    virtual double getMinTimeStepSize() const override;

    virtual void updateNextTimeStepSize( double value ) override;

    virtual void initInitialTimeStamp(double value) override;

    virtual void startNewTimeStep() override;

    virtual double getNextMinTimeStepSize() const override;
};


#endif
