#ifndef EXAHYPE_TIMESTEPPING_SYNCHRONISATION_H_
#define EXAHYPE_TIMESTEPPING_SYNCHRONISATION_H_

#include "exahype/solvers/Solve.h"
#include "exahype/records/ADERDGCellDescription.h"

namespace exahype {
  // forward declarations
  namespace solvers {
    class Solve;
  }
  namespace records {
    class ADERDGCellDescription;
  }

  namespace timestepping {
    void synchroniseTimeStepping(const exahype::solvers::Solve& solve,exahype::records::ADERDGCellDescription& p);
  }
}

#endif /* EXAHYPE_TIMESTEPPING_SYNCHRONISATION_H_ */
