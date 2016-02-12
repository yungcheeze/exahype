//$<EXAHYPE_SOURCE_FILE_COPYRIGHT_NOTE>$
#include "exahype/solvers/Solve.h"

#include <limits>

exahype::solvers::Solve::Solve (
    int solverNumber,
    int parentSolve,
    exahype::solvers::Solve::Type type,
    exahype::solvers::Solve::TimeStepping timeStepping,
    bool sameTimeStepSize,
    bool active
):
  _solverNumber    (solverNumber),
  _parentSolve     (parentSolve),
  _timeStepping    (timeStepping),
  _sameTimeStepSize(sameTimeStepSize)
{
  _active                    = active;
  _type                      = type;
  _predictorMinTimeStamp     = std::numeric_limits<double>::max();
  _correctorMinTimeStamp     = std::numeric_limits<double>::max();
  _correctorTimeStepSize     = std::numeric_limits<double>::max();
  _predictorTimeStepSize     = std::numeric_limits<double>::max();
  _nextPredictorTimeStepSize = std::numeric_limits<double>::max();
}

exahype::solvers::Solve::Solve (const exahype::solvers::Solve& anotherSolve):
  _solverNumber    (anotherSolve._solverNumber),
  _timeStepping    (anotherSolve._timeStepping),
  _sameTimeStepSize(anotherSolve._sameTimeStepSize)
{
  _parentSolve               = anotherSolve._parentSolve;
  _active                    = anotherSolve._active;
  _type                      = anotherSolve._type;
  _timeStepping              = anotherSolve._timeStepping;
  _correctorMinTimeStamp     = anotherSolve._correctorMinTimeStamp;
  _correctorTimeStepSize     = anotherSolve._correctorTimeStepSize;
  _predictorMinTimeStamp     = anotherSolve._predictorMinTimeStamp;
  _predictorTimeStepSize     = anotherSolve._predictorTimeStepSize;
  _nextPredictorTimeStepSize = anotherSolve._nextPredictorTimeStepSize;
}

const int exahype::solvers::Solve::getSolverNumber () const {
  return _solverNumber;
}

void exahype::solvers::Solve::setActive(bool active) {
  _active = active;
}

bool exahype::solvers::Solve::isActive() const {
  return _active;
}

int exahype::solvers::Solve::getParentSolve () const {
  return _parentSolve;
}

void exahype::solvers::Solve::setParentSolve (int parentSolve) {
  _parentSolve = parentSolve;
}

void exahype::solvers::Solve::setType (Type type) {
  _type = type;
}

exahype::solvers::Solve::Type exahype::solvers::Solve::getType () const {
  return _type;
}

exahype::solvers::Solve::TimeStepping exahype::solvers::Solve::getTimeStepping() const {
  return _timeStepping;
}

bool exahype::solvers::Solve::isSameTimeStepSize() const {
  return _sameTimeStepSize;
}

double exahype::solvers::Solve::getPredictorMinTimeStamp () const {
  return _predictorMinTimeStamp;
}

void exahype::solvers::Solve::setPredictorMinTimeStamp (
    double predictorMinTimeStamp) {
  _predictorMinTimeStamp = predictorMinTimeStamp;
}

double exahype::solvers::Solve::getCorrectorMinTimeStamp () const {
  if (_sameTimeStepSize)
    return _predictorMinTimeStamp;

  return _correctorMinTimeStamp;
}

void exahype::solvers::Solve::setCorrectorMinTimeStamp (double correctorMinTimeStamp) {
  _correctorMinTimeStamp = correctorMinTimeStamp;
}

double exahype::solvers::Solve::getCorrectorTimeStepSize() const {
  if (_sameTimeStepSize)
    return _predictorTimeStepSize;

  return _correctorTimeStepSize;
}

double exahype::solvers::Solve::getPredictorTimeStepSize() const {
  return _predictorTimeStepSize;
}

void exahype::solvers::Solve::updateNextPredictorTimeStepSize(const double& nextPredictorTimeStepSize) {
  _nextPredictorTimeStepSize = std::min( _nextPredictorTimeStepSize, nextPredictorTimeStepSize );
}

void exahype::solvers::Solve::startNewTimeStep() {
  _correctorMinTimeStamp     = _predictorMinTimeStamp;
  _correctorTimeStepSize     = _predictorTimeStepSize;

  _predictorTimeStepSize     = _nextPredictorTimeStepSize;
  _predictorMinTimeStamp     = _predictorMinTimeStamp+_nextPredictorTimeStepSize;

  _nextPredictorTimeStepSize = std::numeric_limits<double>::max();
}
