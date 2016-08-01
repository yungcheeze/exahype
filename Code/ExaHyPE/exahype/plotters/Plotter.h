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


/**
 * Central plotter class
 *
 * The ExaHyPE kernel holds one plotter instance per plotter specified in the
 * config file.
 *
 * @author Tobias Weinzierl
 */
class exahype::plotters::Plotter {
 public:

  class UserOnTheFlyPostProcessing {
    public:
      virtual ~UserOnTheFlyPostProcessing() {}

      virtual void startPlotting( double time ) = 0;
      virtual void finishPlotting() = 0;

      virtual void mapQuantities(
        const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
        const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        double* Q,
        double* outputQuantities,
        double timeStamp) = 0;
  };

  /**
   * Actual device chosen by a plotter in the config file. If you implement
   * your own device, please also add a
   *
   * static std::string getIdentifier();
   *
   * routine that returns the device's identifier.
   */
  class Device {
   protected:
    UserOnTheFlyPostProcessing*  _postProcessing;
   public:
    Device(UserOnTheFlyPostProcessing* postProcessing):
      _postProcessing(postProcessing) {}
    virtual ~Device() {}

    /**
     * Configure the plotter. Is invoked directly after the constructor is
     * called.
     */
    virtual void init(const std::string& filename, int order, int unknowns, int writtenUnknowns, const std::string& select) = 0;

    /**
     * Hand a patch over to the plotter. Feel free to ignore the passed data if
     * you don't want to plot it.
     */
    virtual void plotPatch(
        const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
        const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
        double timeStamp) = 0;

    virtual void startPlotting( double time ) = 0;
    virtual void finishPlotting() = 0;
  };

 private:
  static tarch::logging::Log _log;

  const int              _solver;
  const std::string      _identifier;
  int                    _writtenUnknowns;
  double                 _time;
  const double           _repeat;
  const std::string      _filename;
  const std::string      _select;
  bool                   _isActive;

  Device*                      _device;

 public:
  /**
   * @param solver Number of the underlying solver. This number is important to
   *               parse the file: the constructor asks parser for the solverth
   *               solver section.
   * @param plotterCount Same story: Required to tell the parser which tag in
   *               the file is to be read.
   */
  Plotter(int solver, int plotterCount, const exahype::Parser& parser, UserOnTheFlyPostProcessing* postProcessing);
  ~Plotter();

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
