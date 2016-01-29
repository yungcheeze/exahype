#include "exahype/solvers/Plotter.h"


exahype::solvers::Plotter::Plotter( int solver, const std::string& identifier, double time, double repeat, const std::string& filename ):
  _solver(solver),
  _identifier(identifier),
  _time(time),
  _repeat(repeat),
  _filename(filename) {
  assertion( _time>=0.0 );
}


bool exahype::solvers::Plotter::isActive( double currentTimeStamp ) const {
  return (_time>=0.0) && tarch::la::greaterEquals( currentTimeStamp, _time );
}


void exahype::solvers::Plotter::finishedPlotting() {
  if (_repeat>0.0) {

  }
  else {
    _time = -1.0;
  }
}


bool exahype::solvers::isAPlotterActive(double currentTimeStep) {
  bool result = false;
  for (
    std::vector<Plotter*>::const_iterator p = RegisteredPlotters.begin();
    p!=RegisteredPlotters.end();
    p++
  ) {
    result |= (*p)->isActive(currentTimeStep);
  }
  return result;
}


void exahype::solvers::finishedPlotting(double currentTimeStep) {
  for (
    std::vector<Plotter*>::iterator p = RegisteredPlotters.begin();
    p!=RegisteredPlotters.end();
    p++
  ) {
    if ((*p)->isActive(currentTimeStep)) {
      (*p)->finishedPlotting();
    }
  }
}
