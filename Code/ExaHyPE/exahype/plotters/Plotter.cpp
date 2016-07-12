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
 
#include "exahype/plotters/Plotter.h"

#include "ADERDG2ProbeAscii.h"
#include "ADERDG2VTKAscii.h"
#include "ADERDG2VTKBinary.h"
#include "exahype/solvers/Solver.h"


std::vector<exahype::plotters::Plotter*> exahype::plotters::RegisteredPlotters;

tarch::logging::Log exahype::plotters::Plotter::_log(
    "exahype::solvers::Plotter");


exahype::plotters::Plotter::Plotter(int solver, int plotterCount,
                                    const exahype::Parser& parser)
    : _solver(solver),
      _identifier(parser.getIdentifierForPlotter(solver, plotterCount)),
      _time(parser.getFirstSnapshotTimeForPlotter(solver, plotterCount)),
      _repeat(parser.getRepeatTimeForPlotter(solver, plotterCount)),
      _filename(parser.getFilenameForPlotter(solver, plotterCount)),
      _select(parser.getSelectorForPlotter(solver, plotterCount)),
      _isActive(false),
      _device(nullptr) {
  if (_time < 0.0) {
    logError("Plotter(...)",
             "plotter's first snapshot time is set to negative value "
                 << _time);
  }
  if (_repeat < 0.0) {
    logError("Plotter(...)", "plotter's repeat time is set to negative value "
                                 << _repeat);
  }
  logInfo("Plotter(...)", "write snapshot to file "
                              << _filename << " every " << _repeat
                              << " time units with first snapshot at " << _time
                              << ". plotter type is " << _identifier);

  assertion(_solver < static_cast<int>(solvers::RegisteredSolvers.size()));
  switch (solvers::RegisteredSolvers[_solver]->getType()) {
    case exahype::solvers::Solver::Type::ADER_DG:
      /**
       * This is actually some kind of switch expression though switches do
       * not work for strings, so we map it onto an if-then-else cascade.
       */
      if (_identifier.compare( ADERDG2VTKAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2VTKAscii();
      }
      else if (_identifier.compare( ADERDG2VTKBinary::getIdentifier() ) == 0) {
        _device = new ADERDG2VTKBinary();
      }
      else if (_identifier.compare( ADERDG2ProbeAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2ProbeAscii();
      }


      if (_device==nullptr) {
        logError(
          "checkWetherSolverBecomesActive(double)",
          "unknown plotter type "
              << _identifier << " for "
              << solvers::RegisteredSolvers[_solver]->getIdentifier()
        );
      }
      else {
        _device->init(
            _filename,
            // Internally, we always use the nodes per coordinate axis, i.e.,
            // "order+1"
            // Consider to pass the nodes per coordinate axis instead of the
            // order
            solvers::RegisteredSolvers[_solver]->getNodesPerCoordinateAxis() -
            1,
            solvers::RegisteredSolvers[_solver]->getNumberOfVariables(),
            _select
        );
      }
    break;
  }
}


exahype::plotters::Plotter::~Plotter() {
  if (_device!=nullptr) {
    delete _device;
    _device = nullptr;
  }
}


bool exahype::plotters::Plotter::checkWetherSolverBecomesActive(double currentTimeStamp) {
  if ((_time >= 0.0) && tarch::la::greaterEquals(currentTimeStamp, _time)) {
    assertion(_device!=nullptr);
    _isActive = true;
    _device->startPlotting(currentTimeStamp);
  }
  return isActive();
}


bool exahype::plotters::Plotter::isActive() const {
  return _isActive;
}


bool exahype::plotters::Plotter::plotDataFromSolver(int solver) const {
  return isActive() && _solver == solver;
}


void exahype::plotters::Plotter::plotPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
  double timeStamp
) {
  assertion(_device != nullptr);
  if (_device!=nullptr) {
    _device->plotPatch(offsetOfPatch, sizeOfPatch, u, timeStamp);
  }
}


void exahype::plotters::Plotter::finishedPlotting() {
  assertion(isActive());
  if (_repeat > 0.0) {
    _time += _repeat;
  } else {
    _time = -1.0;
  }
  if (_device!=nullptr) {
    _device->finishPlotting();
  }
  _isActive = false;
}


bool exahype::plotters::isAPlotterActive(double currentTimeStep) {
  bool result = false;
  for (const auto& p : RegisteredPlotters) {
    result |= p->checkWetherSolverBecomesActive(currentTimeStep);
  }
  return result;
}


void exahype::plotters::finishedPlotting() {
  for (auto& p : RegisteredPlotters) {
    if (p->isActive()) {
      p->finishedPlotting();
    }
  }
}


std::string exahype::plotters::Plotter::getFileName() const {
  return _filename;
}
