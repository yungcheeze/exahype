#include "exahype/solvers/Solver.h"


std::vector<exahype::solvers::Solver*>  exahype::solvers::RegisteredSolvers;


exahype::solvers::Solver::Solver(const std::string& identifier, Type type, int kernelNumber, int numberOfVariables, int nodesPerCoordinateAxis, TimeStepping timeStepping):
  _identifier                  (identifier),
  _type                        (type),
  _kernelNumber                (kernelNumber),
  _numberOfVariables           (numberOfVariables),
  _nodesPerCoordinateAxis      (nodesPerCoordinateAxis),
  _unknownsPerFace             (numberOfVariables * power(nodesPerCoordinateAxis,DIMENSIONS-1)),
  _unknownsPerCellBoundary     (DIMENSIONS_TIMES_TWO * _unknownsPerFace),
  _unknownsPerCell             (numberOfVariables * power(nodesPerCoordinateAxis,DIMENSIONS+0)),
  _fluxUnknownsPerCell         (_unknownsPerCell*DIMENSIONS),
  _spaceTimeUnknownsPerCell    (numberOfVariables * power(nodesPerCoordinateAxis,DIMENSIONS+1)),
  _spaceTimeFluxUnknownsPerCell(_spaceTimeUnknownsPerCell*DIMENSIONS),
  _timeStepping                (timeStepping),
  _minTimeStamp       (std::numeric_limits<double>::max()),
  _minTimeStepSize    (std::numeric_limits<double>::max()),
  _minNextTimeStepSize(std::numeric_limits<double>::max())
{
  // do nothing
}


std::string exahype::solvers::Solver::getIdentifier() const {
  return _identifier;
}


exahype::solvers::Solver::Type exahype::solvers::Solver::getType() const {
  return _type;
}


int exahype::solvers::Solver::getNumberOfVariables() const {
  return _numberOfVariables;
}


int exahype::solvers::Solver::getNodesPerCoordinateAxis() const {
  return _nodesPerCoordinateAxis;
}


int exahype::solvers::Solver::getUnknownsPerFace() const {
  return _unknownsPerFace;
}


int exahype::solvers::Solver::getUnknownsPerCellBoundary() const {
  return _unknownsPerCellBoundary;
}


int exahype::solvers::Solver::getUnknownsPerCell() const {
  return _unknownsPerCell;
}


int exahype::solvers::Solver::getFluxUnknownsPerCell() const {
  return _fluxUnknownsPerCell;
}


int exahype::solvers::Solver::getSpaceTimeUnknownsPerCell() const {
  return _spaceTimeUnknownsPerCell;
}


int exahype::solvers::Solver::getSpaceTimeFluxUnknownsPerCell() const {
  return _spaceTimeFluxUnknownsPerCell;
}


void exahype::solvers::Solver::synchroniseTimeStepping(exahype::records::ADERDGCellDescription& p) const {
  if (_timeStepping==GlobalTimeStepping) {
    p.setTimeStamp   (_minTimeStamp);
    p.setTimeStepSize(_minTimeStepSize);
  }
}


void exahype::solvers::Solver::startNewTimeStep() {
  _minTimeStepSize     = _minNextTimeStepSize;
  _minTimeStamp        = _minTimeStamp+_minNextTimeStepSize;

  _minNextTimeStepSize = std::numeric_limits<double>::max();
}

void exahype::solvers::Solver::updateMinNextTimeStepSize (const double& minNextTimeStepSize) {
  _minNextTimeStepSize = std::min( _minNextTimeStepSize, minNextTimeStepSize );
}

double exahype::solvers::Solver::getMinNextTimeStepSize () const {
  return _minNextTimeStepSize;
}

void exahype::solvers::Solver::setMinTimeStamp(double minTimeStamp) {
  _minTimeStamp = minTimeStamp;
}

double exahype::solvers::Solver::getMinTimeStamp() const {
  return _minTimeStamp;
}

double exahype::solvers::Solver::getMinTimeStepSize() const {
  return _minTimeStepSize;
}
