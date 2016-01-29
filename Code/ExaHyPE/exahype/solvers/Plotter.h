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
  }
}


class exahype::solvers::Plotter {
  public:
    Plotter( int solver, const std::string& identifier, double time, double repeat, const std::string& filename );
};

#endif
