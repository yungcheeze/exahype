#ifndef _EXAHYPE_SOLVERS_PLOTTER_H_
#define _EXAHYPE_SOLVERS_PLOTTER_H_

#include <string>
#include <vector>


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
  protected:
    class Device {
      public:
        virtual ~Device() {}

        virtual void plotPatch(
          const tarch::la::Vector<DIMENSIONS,double>&  offsetOfPatch,
          const tarch::la::Vector<DIMENSIONS,double>&  sizeOfPatch,
          double* u
        ) = 0;
    };
  private:
    static tarch::logging::Log  _log;

    const int          _solver;
    const std::string  _identifier;
    double             _time;
    const double       _repeat;
    const std::string  _filename;

    Device*            _device;
  public:
    Plotter( int solver, const std::string& identifier, double time, double repeat, const std::string& filename );

    /**
     * Checks whether there should be a plotter according to this class.
     * If it should become open, it is opened
     */
    bool checkWetherSolverBecomesActive( double currentTimeStamp );
    bool isActive() const;
    bool plotDataFromSolver( int solver ) const;
    void finishedPlotting();

    void plotPatch(
      const tarch::la::Vector<DIMENSIONS,double>&  offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS,double>&  sizeOfPatch,
      double* u
    );
};

#endif
