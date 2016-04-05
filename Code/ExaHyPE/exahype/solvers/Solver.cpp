#include "exahype/solvers/Solver.h"

std::vector<exahype::solvers::Solver*> exahype::solvers::RegisteredSolvers;

exahype::solvers::Solver::Solver(const std::string& identifier, Type type,
                                 int kernelNumber, int numberOfVariables,
                                 int nodesPerCoordinateAxis,
                                 TimeStepping timeStepping)
: _identifier(identifier),
  _type(type),
  _kernelNumber(kernelNumber),
  _numberOfVariables(numberOfVariables),
  _nodesPerCoordinateAxis(nodesPerCoordinateAxis),
  _unknownsPerFace(numberOfVariables *
                   power(nodesPerCoordinateAxis, DIMENSIONS - 1)),
  _unknownsPerCellBoundary(DIMENSIONS_TIMES_TWO * _unknownsPerFace),
  _unknownsPerCell(numberOfVariables *
                         power(nodesPerCoordinateAxis, DIMENSIONS + 0)),
  _fluxUnknownsPerCell(_unknownsPerCell * DIMENSIONS),
  _spaceTimeUnknownsPerCell(numberOfVariables *
                            power(nodesPerCoordinateAxis, DIMENSIONS + 1)),
  _spaceTimeFluxUnknownsPerCell(_spaceTimeUnknownsPerCell * DIMENSIONS),
  _timeStepping(timeStepping),
  _minCorrectorTimeStamp(std::numeric_limits<double>::max()),
  _minPredictorTimeStamp(std::numeric_limits<double>::max()),
  _minCorrectorTimeStepSize(std::numeric_limits<double>::max()),
  _minPredictorTimeStepSize(std::numeric_limits<double>::max()),
  _minNextPredictorTimeStepSize(std::numeric_limits<double>::max()) {
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

/*
 * \todo 16/04/02:Dominic Etienne Charrier:
 * non-virtual method is only for now; The refinement criterion must be
 * specified by the user. So replace method by virtual one later.
 */
bool exahype::solvers::Solver::refinementCriterion(
    const double* luh,
    const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    double t,
    const int level) {
  assertion(level>=getMinimumTreeDepth()+1);

  if (level<getMinimumTreeDepth()+4) {
    if (center[0]<0.5 && center[0]>0.25) {
      if (center[1]<0.5 && center[1]>0.25) {
        return true;
      }
    }
  }
  return false;
}

void exahype::solvers::Solver::synchroniseTimeStepping(
    exahype::records::ADERDGCellDescription& p) const {
  // todo 16/02/27:Dominic Etienne Charrier
  // in case we use optimistic time stepping:
  // if last predictor time step size is larger
  // as admissibleTimeStepSize + tolerance:
  // make sure that corrector time step size
  // will equal predictor time step size in next
  // sweep.
  // Extra attention must be paid to time stamps.
  // All this should be done by the solver.

  //  if (p.getNextPredictorTimeStepSize() < )

  if (_timeStepping == GlobalTimeStepping) {
    p.setCorrectorTimeStamp(_minCorrectorTimeStamp);
    p.setCorrectorTimeStepSize(_minCorrectorTimeStepSize);
    p.setPredictorTimeStamp(_minPredictorTimeStamp);
    p.setPredictorTimeStepSize(_minPredictorTimeStepSize);

    /*
        assertionNumericalEquals1(p.getCorrectorTimeStamp()
       ,solve.getCorrectorTimeStamp(),   1e-12); // todo precision
        assertionNumericalEquals1(p.getCorrectorTimeStepSize(),solve.getCorrectorTimeStepSize(),1e-12);
        assertionNumericalEquals1(p.getPredictorTimeStamp()
       ,solve.getPredictorTimeStamp(),   1e-12);
        assertionNumericalEquals1(p.getPredictorTimeStepSize(),solve.getPredictorTimeStepSize(),1e-12);
     */
  }
  /*  if (!solve.isCorrectorTimeLagging()) {*/
  //    p.setCorrectorTimeStamp   (p.getPredictorTimeStamp   ());
  //    p.setCorrectorTimeStepSize(p.getPredictorTimeStepSize());
  //  }

#if defined(Debug) || defined(Asserts)
  // @ŧodo Wieder reinnehmen
  /*
  if (solver.getTimeStepping()==exahype::solvers::Solver::GLOBAL &&
  !solver.isCorrectorTimeLagging()) {
    // Note that the solve time stamps and time step sizes are not modified if
  corrector time lagging
    // is deactivated. Thus, solve.getPredictor... and solve.getCorrector... are
  not the same in general
    // for any value of solve.isCorrectorTimeLagging().
    assertionNumericalEquals1(p.getPredictorTimeStamp()
  ,solver.getPredictorTimeStamp(),   1e-12); // todo precision
    assertionNumericalEquals1(p.getPredictorTimeStepSize(),solver.getPredictorTimeStepSize(),1e-12);
    assertionNumericalEquals1(p.getCorrectorTimeStamp()
  ,solver.getPredictorTimeStamp(),   1e-12);
    assertionNumericalEquals1(p.getCorrectorTimeStepSize(),solver.getPredictorTimeStepSize(),1e-12);
  }
   */
#endif
}

void exahype::solvers::Solver::startNewTimeStep() {
  _minCorrectorTimeStamp = _minPredictorTimeStamp;
  _minCorrectorTimeStepSize = _minPredictorTimeStepSize;

  _minPredictorTimeStepSize = _minNextPredictorTimeStepSize;
  _minPredictorTimeStamp =
      _minPredictorTimeStamp + _minNextPredictorTimeStepSize;

  _minNextPredictorTimeStepSize = std::numeric_limits<double>::max();
}

void exahype::solvers::Solver::updateMinNextPredictorTimeStepSize(
    const double& minNextPredictorTimeStepSize) {
  _minNextPredictorTimeStepSize =
      std::min(_minNextPredictorTimeStepSize, minNextPredictorTimeStepSize);
}

double exahype::solvers::Solver::getMinNextPredictorTimeStepSize() const {
  return _minNextPredictorTimeStepSize;
}

void exahype::solvers::Solver::setMinCorrectorTimeStamp(
    double minCorrectorTimeStamp) {
  _minCorrectorTimeStamp = minCorrectorTimeStamp;
}

double exahype::solvers::Solver::getMinCorrectorTimeStamp() const {
  return _minCorrectorTimeStamp;
}

void exahype::solvers::Solver::setMinPredictorTimeStamp(
    double minPredictorTimeStamp) {
  _minPredictorTimeStamp = minPredictorTimeStamp;
}

double exahype::solvers::Solver::getMinPredictorTimeStamp() const {
  return _minPredictorTimeStamp;
}

double exahype::solvers::Solver::getMinCorrectorTimeStepSize() const {
  return _minCorrectorTimeStepSize;
}

void exahype::solvers::Solver::setMinPredictorTimeStepSize(
    double minPredictorTimeStepSize) {
  _minPredictorTimeStepSize = minPredictorTimeStepSize;
}

double exahype::solvers::Solver::getMinPredictorTimeStepSize() const {
  return _minPredictorTimeStepSize;
}
