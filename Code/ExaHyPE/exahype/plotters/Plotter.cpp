#include "exahype/plotters/Plotter.h"
#include "exahype/solvers/Solver.h"


std::vector<exahype::plotters::Plotter*>  exahype::plotters::RegisteredPlotters;


tarch::logging::Log  exahype::plotters::Plotter::_log( "exahype::solvers::Plotter" );


exahype::plotters::Plotter::Plotter( int solver, const std::string& identifier, double time, double repeat, const std::string& filename ):
  _solver(solver),
  _identifier(identifier),
  _time(time),
  _repeat(repeat),
  _filename(filename),
  _device(nullptr) {
  assertion( _time>=0.0 );
}


bool exahype::plotters::Plotter::checkWetherSolverBecomesActive( double currentTimeStamp ) {
  if ( (_time>=0.0) && tarch::la::greaterEquals( currentTimeStamp, _time ) ) {
    assertion( _solver<solvers::RegisteredSolvers.size() );
    switch ( solvers::RegisteredSolvers[ _solver ]->getType() ) {
      case solvers::Solver::ADER_DG:
        if ( _identifier.compare( "vtk::binary" )==0 ) {
          logInfo( "open()", "create vtk::binary plotter for " << solvers::RegisteredSolvers[ _solver ]->getIdentifier() );
        }
        else {
          logError( "open()", "unknown plotter type " << _identifier << " for " << solvers::RegisteredSolvers[ _solver ]->getIdentifier() );
        }
        break;
    }
  }
  return isActive();
}


bool exahype::plotters::Plotter::isActive() const {
  return _device != nullptr;
}


bool exahype::plotters::Plotter::plotDataFromSolver( int solver ) const {
  return isActive() && _solver==solver;
}


void exahype::plotters::Plotter::plotPatch(
  const tarch::la::Vector<DIMENSIONS,double>&  offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS,double>&  sizeOfPatch,
  double* u
) {

}


void exahype::plotters::Plotter::finishedPlotting() {
  assertion(isActive());
  if (_repeat>0.0) {
    _time += _repeat;
  }
  else {
    _time = -1.0;
  }
  delete _device;
  _device = nullptr;
}


bool exahype::plotters::isAPlotterActive(double currentTimeStep) {
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


void exahype::plotters::finishedPlotting() {
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
