#ifndef _EXAHYPE_SOLVERS_PLOTTER_H_
#define _EXAHYPE_SOLVERS_PLOTTER_H_

#include <string>
#include <vector>


#include "peano/utils/Globals.h"
#include "tarch/la/Vector.h"


namespace exahype {
  namespace solvers {
    class Plotter;

    extern std::vector<Plotter*> RegisteredPlotters;

    bool isAPlotterActive(double currentTimeStep);
    void finishedPlotting();
  }
}


class exahype::solvers::Plotter {
  private:
    const int          _solver;
    const std::string  _identifier;
    double             _time;
    const double       _repeat;
    const std::string  _filename;
    bool               _isActive;
  public:
    Plotter( int solver, const std::string& identifier, double time, double repeat, const std::string& filename );

    /**
     * Checks whether there should be a plotter according to this class.
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

    void open();
    void close();
};

#endif
