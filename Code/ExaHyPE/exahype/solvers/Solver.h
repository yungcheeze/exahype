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
  public:
    struct Plot {
      int          variable;
      double       nextSnapshot;
      bool         repeat;
      std::string  filename;

      Plot( int variable_, double nextSnapshot_, bool repeat_, const std::string& filename_);
    };


    // @todo TW Raus brauchen wir nicht
    enum Type {
      ADER_DG
    };
  protected:
    /**
     * Each solver has an identifier/name. It is used for debug purposes only.
     */
    const std::string _identifier;

    const Type        _type;

    /**
     * Each solver has a kernel number that says which kernel is to be
     * executed. Typically this is an ascending index starting from 0.
     */
    const int         _kernelNumber;

    const int         _numberOfVariables;

    const int         _nodesPerCoordinateAxis;
  public:
    Solver(const std::string& identifier, Type type, int kernelNumber, int numberOfVariables, int nodesPerCoordinateAxis);

    virtual ~Solver() {}

    /**
     * Identify minimal mesh width at a certain point in the domain. This
     * minimal mesh width is used both as a constraint on the AMR as well
     * as to set up the initial grid. If you return 0, you indicate that
     * this PDE might not exist in the domain.
     */
    virtual int getMinimumTreeDepth() const = 0;

    Type getType() const;

    int getNumberOfVariables() const;

    /**
     * If you use a higher order method, then this operation returns the
     * polynomial degree plus one. If you use a Finite Volume method, it
     * returns the number of cells within a patch per coordinate axis.
     */
    int getNodesPerCoordinateAxis() const;



    // Nur der ADER-DG Loeser
    virtual void spaceTimePredictor(
        double * lQi,
        double * lFi,
        const double * const luh, // const
        double * lQhi,
        double * lFhi,
        double * lQhbnd,
        double * lFhbnd,
        const tarch::la::Vector<DIMENSIONS,double>& dx,
        const double dt
    ) = 0;

};

#endif

