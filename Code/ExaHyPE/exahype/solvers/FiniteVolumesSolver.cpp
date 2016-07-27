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
  _unknownsPerCell( (numberOfVariables+numberOfParameters) * power(nodesPerCoordinateAxis, DIMENSIONS + 0)),
  _minTimeStamp(std::numeric_limits<double>::max()),
  _minTimeStepSize(std::numeric_limits<double>::max()),
  _nextMinTimeStepSize(0.0) {
  assertion3(_unknownsPerCell>0, numberOfVariables, numberOfParameters, nodesPerCoordinateAxis);

}


int exahype::solvers::FiniteVolumesSolver::getUnknownsPerCell() const {
  return _unknownsPerCell;
}


double exahype::solvers::FiniteVolumesSolver::getMinTimeStamp() const {
  return _minTimeStamp;
}


double exahype::solvers::FiniteVolumesSolver::getMinTimeStepSize() const {
  return _minTimeStepSize;
}


void exahype::solvers::FiniteVolumesSolver::updateNextTimeStepSize( double value ) {
  _nextMinTimeStepSize = std::min(_nextMinTimeStepSize,value);
}


void exahype::solvers::FiniteVolumesSolver::initInitialTimeStamp(double value) {
  _minTimeStamp = value;
}


void exahype::solvers::FiniteVolumesSolver::startNewTimeStep() {
  _minTimeStamp       += _minTimeStepSize;
  _minTimeStepSize     = _nextMinTimeStepSize;
  _nextMinTimeStepSize = std::numeric_limits<double>::max();
}


double exahype::solvers::FiniteVolumesSolver::getNextMinTimeStepSize() const {
  return _nextMinTimeStepSize;
}
