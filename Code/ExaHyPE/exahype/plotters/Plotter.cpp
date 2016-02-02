#include "exahype/plotters/Plotter.h"
#include "exahype/solvers/Solver.h"

#include "exahype/plotters/ADERDG2BinaryVTK.h"
#include "exahype/plotters/ADERDG2AsciiVTK.h"


std::vector<exahype::plotters::Plotter*>  exahype::plotters::RegisteredPlotters;


tarch::logging::Log  exahype::plotters::Plotter::_log( "exahype::solvers::Plotter" );


/*
exahype::plotters::Plotter::Plotter( int solver, const std::string& identifier, double time, double repeat, const std::string& filename ):
  _solver(solver),
  _identifier(identifier),
  _time(time),
  _repeat(repeat),
  _filename(filename),
  _device(nullptr) {
  assertion( _time>=0.0 );
}
*/

exahype::plotters::Plotter::Plotter( int solver, int plotterCount, const exahype::Parser& parser ):
  _solver(solver),
  _identifier(parser.getIdentifierForPlotter(solver,plotterCount)),
  _time(parser.getFirstSnapshotTimeForPlotter(solver,plotterCount)),
  _repeat(parser.getRepeatTimeForPlotter(solver,plotterCount)),
  _filename(parser.getFilenameForPlotter(solver,plotterCount)),
  _device(nullptr) {
 if (_time<0.0) {
   logError( "Plotter(...)", "plotter's first snapshot time is set to negative value " << _time );
 }
 if (_repeat<0.0) {
   logError( "Plotter(...)", "plotter's repeat time is set to negative value " << _repeat );
 }
 logInfo( "Plotter(...)", "write snapshot to file " << _filename << " every " << _repeat << " time units with first snapshot at " << _time << ". plotter type is " << _identifier );
}


bool exahype::plotters::Plotter::checkWetherSolverBecomesActive( double currentTimeStamp ) {
  if ( (_time>=0.0) && tarch::la::greaterEquals( currentTimeStamp, _time ) ) {
    assertion( _solver < static_cast<int>(solvers::RegisteredSolvers.size()) );
    switch ( solvers::RegisteredSolvers[ _solver ]->getType() ) {
      case solvers::Solver::ADER_DG:
        if ( _identifier.compare( "vtk::binary" )==0 ) {
          logInfo( "open()", "create vtk::binary plotter for " << solvers::RegisteredSolvers[ _solver ]->getIdentifier() );
          _device = new ADERDG2BinaryVTK( _filename, solvers::RegisteredSolvers[ _solver ]->getNodesPerCoordinateAxis()-1, solvers::RegisteredSolvers[ _solver ]->getNumberOfVariables() );
        }
        else if( _identifier.compare( "vtk::ascii" )==0 ) {
          logInfo( "open()", "create vtk::ascii plotter for " << solvers::RegisteredSolvers[ _solver ]->getIdentifier() );
          _device = new ADERDG2AsciiVTK( _filename, solvers::RegisteredSolvers[ _solver ]->getNodesPerCoordinateAxis()-1, solvers::RegisteredSolvers[ _solver ]->getNumberOfVariables() );
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
  double* u, double timeStamp
) {
  assertion( _device!=nullptr );
  _device->plotPatch( offsetOfPatch, sizeOfPatch, u, timeStamp );
}


void exahype::plotters::Plotter::finishedPlotting() {
  assertion(isActive());
  if (_repeat>0.0) {
    std::cerr << _repeat << std::endl;
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
