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

    /**
     * Total number of unknowns per cell face.
     */
    int _unknownsPerFace;

    /**
     * Total number of unknowns per cell boundary.
     */
    int _unknownsPerCellBoundary;

    /**
     * Minimum time stamps of all patches.
     */
    double _minTimeStamp;

    /**
     * Minimum time step size of all patches.
     */
    double _minTimeStepSize;

    /**
     * Next minimum step size of all patches.
     */
    double _nextMinTimeStepSize;

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

    /**
     * @param luh is a pointer to 3^d pointers to doubles
     */
    virtual double stableTimeStepSize(
        double* luh[THREE_POWER_D],
        const tarch::la::Vector<DIMENSIONS, double>& dx) = 0;

    /**
     *
     */
    virtual void solutionAdjustment( double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) = 0;

    virtual bool hasToAdjustSolution(
        const tarch::la::Vector<DIMENSIONS, double>& center,
        const tarch::la::Vector<DIMENSIONS, double>& dx, double t) = 0;

    virtual exahype::solvers::Solver::RefinementControl refinementCriterion(
        const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,
        const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
        const int level) = 0;

    /**
     * @param luh is a pointer to 3^d pointers to doubles
     * @param dt Time step size that is to be used.
     * @param maxAdmissibleDt Maximum time step size that would have been
     *        possible. If maxAdmissibleDt<dt, then we know that no time
     *        step has been done.
     */
    virtual void solutionUpdate(double* luh[THREE_POWER_D], const tarch::la::Vector<DIMENSIONS, double>& dx, const double dt, double& maxAdmissibleDt) = 0;

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

    void sendToRank(int rank, int tag) override;

    void receiveFromMasterRank(int rank, int tag) override;
};


#endif
