/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released unter the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#ifndef _EXAHYPE_PLOTTERS_PLOTTER_H_
#define _EXAHYPE_PLOTTERS_PLOTTER_H_

#include <string>
#include <vector>

#include "exahype/Parser.h"
#include "peano/utils/Globals.h"

#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

namespace exahype {
namespace plotters {
class Plotter;

extern std::vector<Plotter*> RegisteredPlotters;

bool isAPlotterActive(double currentTimeStep);
void finishedPlotting();
}
}

class exahype::plotters::Plotter {
 public:
  class Device {
   public:
    virtual ~Device() {}

    virtual void plotPatch(
        const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
        const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
        double timeStamp) = 0;
  };

 private:
  static tarch::logging::Log _log;

  const int _solver;
  const std::string _identifier;
  double _time;
  const double _repeat;
  const std::string _filename;

  Device* _device;

 public:
  Plotter(int solver, int plotterCount, const exahype::Parser& parser);

  // Disallow copy and assignment
  Plotter(const Plotter& other) = delete;
  Plotter& operator=(const Plotter& other) = delete;

  /**
   * Checks whether there should be a plotter according to this class.
   * If it should become open, it is opened
   */
  bool checkWetherSolverBecomesActive(double currentTimeStamp);
  bool isActive() const;
  bool plotDataFromSolver(int solver) const;
  void finishedPlotting();

  void plotPatch(const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
                 const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
                 double* u, double timeStamp);

  std::string getFileName() const;
};

#endif
