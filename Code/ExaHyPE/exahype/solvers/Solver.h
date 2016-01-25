#ifndef _EXAHYPE_SOLVERS_SOLVER_H_
#define _EXAHYPE_SOLVERS_SOLVER_H_


#include <string>
#include <vector>


#include "peano/utils/Globals.h"
#include "tarch/la/Vector.h"


namespace exahype {
  namespace solvers {
    class Solver;

    extern std::vector<Solver*> RegisteredSolvers;
  }
}



/**
 * Describes one solver.
 */
class exahype::solvers::Solver {
  protected:
    /**
     * Each solver has an identifier/name. It is used for debug purposes only.
     */
    const std::string _identifier;

    /**
     * Each solver has a kernel number that says which kernel is to be
     * executed. Typically this is an ascending index starting from 0.
     */
    const int         _kernelNumber;
  public:
    Solver(const std::string& identifier, int kernelNumber);

    /**
     * Identify minimal mesh width at a certain point in the domain. This
     * minimal mesh width is used both as a constraint on the AMR as well
     * as to set up the initial grid. If you return 0, you indicate that
     * this PDE might not exist in the domain.
     */
    virtual int getMinimumTreeDepth() const = 0;
};

#endif

