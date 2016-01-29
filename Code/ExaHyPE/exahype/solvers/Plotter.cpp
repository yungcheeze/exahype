#include "exahype/solvers/Plotter.h"


std::vector<exahype::solvers::Plotter*>  exahype::solvers::RegisteredPlotters;


exahype::solvers::Plotter::Plotter( int solver, const std::string& identifier, double time, double repeat, const std::string& filename ):
  _solver(solver),
  _identifier(identifier),
  _time(time),
  _repeat(repeat),
  _filename(filename),
  _isActive(false) {
  assertion( _time>=0.0 );
}


bool exahype::solvers::Plotter::checkWetherSolverBecomesActive( double currentTimeStamp ) {
  _isActive = (_time>=0.0) && tarch::la::greaterEquals( currentTimeStamp, _time );
  return _isActive;
}


bool exahype::solvers::Plotter::isActive() const {
  return _isActive;
}


bool exahype::solvers::Plotter::plotDataFromSolver( int solver ) const {
  return _isActive && _solver==solver;
}


void exahype::solvers::Plotter::plotPatch(
  const tarch::la::Vector<DIMENSIONS,double>&  offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS,double>&  sizeOfPatch,
  double* u
) {

}


void exahype::solvers::Plotter::finishedPlotting() {
  if (_repeat>0.0) {
    _time += _repeat;
  }
  else {
    _time = -1.0;
  }
  _isActive = false;
}


void exahype::solvers::Plotter::open() {

}


void exahype::solvers::Plotter::close() {

}


bool exahype::solvers::isAPlotterActive(double currentTimeStep) {
  bool result = false;
  for (
    std::vector<Plotter*>::const_iterator p = RegisteredPlotters.begin();
    p!=RegisteredPlotters.end();
    p++
  ) {
    result |= (*p)->checkWetherSolverBecomesActive(currentTimeStep);
  }
  return result;
}


void exahype::solvers::finishedPlotting() {
  for (
    std::vector<Plotter*>::iterator p = RegisteredPlotters.begin();
    p!=RegisteredPlotters.end();
    p++
  ) {
    if ((*p)->isActive()) {
      (*p)->finishedPlotting();
    }
  }
}
