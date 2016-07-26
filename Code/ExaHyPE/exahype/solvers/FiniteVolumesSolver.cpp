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

#include "exahype/solvers/FiniteVolumesSolver.h"


exahype::solvers::FiniteVolumesSolver::FiniteVolumesSolver(
  const std::string& identifier,
  int numberOfVariables,
  int numberOfParameters,
  int nodesPerCoordinateAxis,
  double maximumMeshSize,
  exahype::solvers::Solver::TimeStepping timeStepping,
  std::unique_ptr<profilers::Profiler> profiler):
  Solver(
    identifier,
    exahype::solvers::Solver::Type::FiniteVolumes,
    numberOfVariables,
    numberOfParameters,
    nodesPerCoordinateAxis,
    maximumMeshSize,
    timeStepping,
    std::move(profiler)
  ),
  _unknownsPerCell( (numberOfVariables+numberOfParameters) * power(nodesPerCoordinateAxis, DIMENSIONS + 0))
  {
  assertion3(_unknownsPerCell>0, numberOfVariables, numberOfParameters, nodesPerCoordinateAxis);

}


double exahype::solvers::FiniteVolumesSolver::stableTimeStepSize(
    const double* const luh,
    const tarch::la::Vector<DIMENSIONS, double>& dx)  {

}


void exahype::solvers::FiniteVolumesSolver::solutionAdjustment(
    double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt)  {

}


bool exahype::solvers::FiniteVolumesSolver::hasToAdjustSolution(
    const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t)  {

}


exahype::solvers::Solver::RefinementControl exahype::solvers::FiniteVolumesSolver::refinementCriterion(
    const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level)  {

}


int exahype::solvers::FiniteVolumesSolver::getUnknownsPerCell() const {
  return _unknownsPerCell;
}


double exahype::solvers::FiniteVolumesSolver::getMinTimeStamp() const {

}


double exahype::solvers::FiniteVolumesSolver::getMinTimeStepSize() const {

}


void exahype::solvers::FiniteVolumesSolver::updateNextTimeStepSize( double value ) {

}


void exahype::solvers::FiniteVolumesSolver::initInitialTimeStamp(double value) {

}


void exahype::solvers::FiniteVolumesSolver::startNewTimeStep() {

}


double exahype::solvers::FiniteVolumesSolver::getNextMinTimeStepSize() const {

}
