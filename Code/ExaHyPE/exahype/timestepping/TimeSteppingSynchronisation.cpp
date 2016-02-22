#include "exahype/timestepping/TimeSteppingSynchronization.h"

#include "exahype/solvers/Solve.h"
#include "exahype/records/ADERDGCellDescription.h"

void exahype::timestepping::synchroniseTimeStepping(const exahype::solvers::Solve& solve,exahype::records::ADERDGCellDescription& p) {
  if (solve.getTimeStepping()==exahype::solvers::Solve::GLOBAL) {
    p.setCorrectorTimeStamp   (solve.getCorrectorTimeStamp   ());
    p.setCorrectorTimeStepSize(solve.getCorrectorTimeStepSize());
    p.setPredictorTimeStamp   (solve.getPredictorTimeStamp   ());
    p.setPredictorTimeStepSize(solve.getPredictorTimeStepSize());

    assertionNumericalEquals1(p.getCorrectorTimeStamp()   ,solve.getCorrectorTimeStamp(),   1e-12); // todo precision
    assertionNumericalEquals1(p.getCorrectorTimeStepSize(),solve.getCorrectorTimeStepSize(),1e-12);
    assertionNumericalEquals1(p.getPredictorTimeStamp()   ,solve.getPredictorTimeStamp(),   1e-12);
    assertionNumericalEquals1(p.getPredictorTimeStepSize(),solve.getPredictorTimeStepSize(),1e-12);
  }
  if (!solve.isCorrectorTimeLagging()) {
    p.setCorrectorTimeStamp   (p.getPredictorTimeStamp   ());
    p.setCorrectorTimeStepSize(p.getPredictorTimeStepSize());
  }

#if defined(Debug) || defined(Asserts)
  if (solve.getTimeStepping()==exahype::solvers::Solve::GLOBAL && !solve.isCorrectorTimeLagging()) {
    // Note that the solve time stamps and time step sizes are not modified if corrector time lagging
    // is deactivated. Thus, solve.getPredictor... and solve.getCorrector... are not the same in general
    // for any value of solve.isCorrectorTimeLagging().
    assertionNumericalEquals1(p.getPredictorTimeStamp()   ,solve.getPredictorTimeStamp(),   1e-12); // todo precision
    assertionNumericalEquals1(p.getPredictorTimeStepSize(),solve.getPredictorTimeStepSize(),1e-12);
    assertionNumericalEquals1(p.getCorrectorTimeStamp()   ,solve.getPredictorTimeStamp(),   1e-12);
    assertionNumericalEquals1(p.getCorrectorTimeStepSize(),solve.getPredictorTimeStepSize(),1e-12);
  }
#endif
}
