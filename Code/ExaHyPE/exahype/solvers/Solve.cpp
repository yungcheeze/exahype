/*
//$<EXAHYPE_SOURCE_FILE_COPYRIGHT_NOTE>$
#include "exahype/solvers/Solve.h"

#include "tarch/Assertions.h"

#include <limits>

exahype::solvers::Solve::Solve (
    std::string name,
    int solverNumber,
    exahype::solvers::Solve::Type type,
    exahype::solvers::Solve::TimeStepping timeStepping,
    bool correctorTimeLagging,
    bool active,
    double correctorTimeStamp,
    double predictorTimeStamp,
    double correctorTimeStepSize,
    double predictorTimeStepSize,
    double nextPredictorTimeStepSize
):
  _name                     (name),
  _solverNumber             (solverNumber),
  _parentSolve              (InvalidParentSolveIdentifier),
  _type                     (type),
  _timeStepping             (timeStepping),
  _correctorTimeLagging     (correctorTimeLagging),
  _active                   (active),
  _correctorTimeStamp       (correctorTimeStamp),
  _predictorTimeStamp       (predictorTimeStamp),
  _correctorTimeStepSize    (correctorTimeStepSize),
  _predictorTimeStepSize    (predictorTimeStepSize),
  _nextPredictorTimeStepSize(nextPredictorTimeStepSize)
{}

exahype::solvers::Solve::Solve (
    int solverNumber,
    exahype::solvers::Solve::Type type,
    exahype::solvers::Solve::TimeStepping timeStepping,
    bool correctorTimeLagging,
    bool active
):
  _solverNumber    (solverNumber),
  _timeStepping    (timeStepping),
  _correctorTimeLagging(correctorTimeLagging)
{
  _parentSolve               = exahype::solvers::Solve::InvalidParentSolveIdentifier;
  _active                    = active;
  _type                      = type;
  _predictorTimeStamp        = std::numeric_limits<double>::max();
  _correctorTimeStamp        = std::numeric_limits<double>::max();
  _correctorTimeStepSize     = std::numeric_limits<double>::max();
  _predictorTimeStepSize     = std::numeric_limits<double>::max();
  _nextPredictorTimeStepSize = std::numeric_limits<double>::max();
}

exahype::solvers::Solve::Solve (const exahype::solvers::Solve& otherSolve):
  _solverNumber    (otherSolve._solverNumber),
  _timeStepping    (otherSolve._timeStepping),
  _correctorTimeLagging(otherSolve._correctorTimeLagging)
{
  _parentSolve               = otherSolve._parentSolve; // todo might be replaced in future versions.
  _active                    = otherSolve._active;
  _type                      = otherSolve._type;
  _timeStepping              = otherSolve._timeStepping;
  _correctorTimeStamp        = otherSolve._correctorTimeStamp;
  _correctorTimeStepSize     = otherSolve._correctorTimeStepSize;
  _predictorTimeStamp        = otherSolve._predictorTimeStamp;
  _predictorTimeStepSize     = otherSolve._predictorTimeStepSize;
  _nextPredictorTimeStepSize = otherSolve._nextPredictorTimeStepSize;
}

std::string exahype::solvers::Solve::getName () const {
  return _name;
}

int exahype::solvers::Solve::getSolverNumber () const {
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

bool exahype::solvers::Solve::isCorrectorTimeLagging() const {
  return _correctorTimeLagging;
}

double exahype::solvers::Solve::getPredictorTimeStamp () const {
  return _predictorTimeStamp;
}

void exahype::solvers::Solve::setPredictorTimeStamp (
    double predictorTimeStamp) {
  _predictorTimeStamp = predictorTimeStamp;
}

double exahype::solvers::Solve::getCorrectorTimeStamp () const {
    return _correctorTimeStamp;
}

void exahype::solvers::Solve::setCorrectorTimeStamp (double correctorTimeStamp) {
  _correctorTimeStamp = correctorTimeStamp;
}

double exahype::solvers::Solve::getCorrectorTimeStepSize() const {
    return _correctorTimeStepSize;
}

void exahype::solvers::Solve::setCorrectorTimeStepSize (double correctorTimeStepSize) {
  _correctorTimeStepSize = correctorTimeStepSize;
}

double exahype::solvers::Solve::getPredictorTimeStepSize() const {
  return _predictorTimeStepSize;
}

void exahype::solvers::Solve::setPredictorTimeStepSize (double predictorTimeStepSize) {
  _predictorTimeStepSize = predictorTimeStepSize;
}

double exahype::solvers::Solve::getNextPredictorTimeStepSize () const {
  return _nextPredictorTimeStepSize;
}

void exahype::solvers::Solve::updateNextPredictorTimeStepSize (const double& nextPredictorTimeStepSize) {
  _nextPredictorTimeStepSize = std::min( _nextPredictorTimeStepSize, nextPredictorTimeStepSize );
}

void exahype::solvers::Solve::startNewTimeStep () {
  _correctorTimeStamp     = _predictorTimeStamp;
  _correctorTimeStepSize  = _predictorTimeStepSize;

  _predictorTimeStepSize  = _nextPredictorTimeStepSize;
  _predictorTimeStamp     = _predictorTimeStamp+_nextPredictorTimeStepSize;

  _nextPredictorTimeStepSize = std::numeric_limits<double>::max();
}

void exahype::solvers::Solve::merge (const exahype::solvers::Solve& otherSolve) {
#if defined(Debug) || defined(Asserts)
  assertionMsg(otherSolve._solverNumber         == _solverNumber        ,"Solver numbers must be the same!");
  assertionMsg(otherSolve._timeStepping         == _timeStepping        ,"Time stepping mode must be the same!");
  assertionMsg(otherSolve._correctorTimeLagging == _correctorTimeLagging,"Time step size selection must be the same!");
#endif
//  _parentSolve               = otherSolve._parentSolve;
//  _type                      = otherSolve._type;
  _active                    = _active & otherSolve._active;
  _correctorTimeStamp        = std::min( _correctorTimeStamp       , otherSolve._correctorTimeStamp        );
  _correctorTimeStepSize     = std::min( _correctorTimeStepSize    , otherSolve._correctorTimeStepSize     );
  _predictorTimeStamp        = std::min( _predictorTimeStamp       , otherSolve._predictorTimeStamp        );
  _predictorTimeStepSize     = std::min( _predictorTimeStepSize    , otherSolve._predictorTimeStepSize     );
  _nextPredictorTimeStepSize = std::min( _nextPredictorTimeStepSize, otherSolve._nextPredictorTimeStepSize );
}
*/
