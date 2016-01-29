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
    void finishedPlotting(double currentTimeStep);
  }
}


class exahype::solvers::Plotter {
  private:
    const int          _solver;
    const std::string  _identifier;
    double             _time;
    const double       _repeat;
    const std::string  _filename;
  public:
    Plotter( int solver, const std::string& identifier, double time, double repeat, const std::string& filename );

    /**
     * Checks whether there should be a plotter according to this class.
     */
    bool isActive( double currentTimeStamp ) const;
    void finishedPlotting();
};

#endif
