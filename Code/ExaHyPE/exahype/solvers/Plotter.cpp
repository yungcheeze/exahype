#include "exahype/solvers/Plotter.h"


exahype::solvers::Plotter::Plotter( int solver, const std::string& identifier, double time, double repeat, const std::string& filename ):
  _solver(solver),
  _identifier(identifier),
  _time(time),
  _repeat(repeat),
  _filename(filename) {
}


bool exahype::solvers::Plotter::isActive( double currentTimeStamp ) const {
  return tarch::la::greaterEquals( currentTimeStamp, _time );
}


bool exahype::solvers::isAPlotterActive(double currentTimeStep) {
  result = false;
  for (
    std::vector<Plotter*>::const_iterator p = RegisteredPlotters.begin();
    p!=RegisteredPlotters.end();
    p++
  ) {
    result |= (*p)->isActive(currentTimeStep);
  }
  return result;
}


void exahype::solvers::finishedPlotting() {

}
