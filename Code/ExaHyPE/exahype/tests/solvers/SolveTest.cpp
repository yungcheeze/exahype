#include "exahype/tests/solvers/SolveTest.h"


#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"

#include "exahype/solvers/Solve.h"

#include <limits>

registerTest(exahype::tests::solvers::SolveTest)


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",off)
#endif


exahype::tests::solvers::SolveTest::SolveTest():
tarch::tests::TestCase ( "exahype::tests::solvers::SolveTest" ) {
}


exahype::tests::solvers::SolveTest::~SolveTest() {
}


void exahype::tests::solvers::SolveTest::run() {
  testMethod( testSolve );
}


void exahype::tests::solvers::SolveTest::testSolve() {
  exahype::solvers::Solve solve(
      0, // solverNumber
      exahype::solvers::Solve::InvalidParentSolveIdentifier,
      exahype::solvers::Solve::SOLVE,
      exahype::solvers::Solve::GLOBAL,
      true, // corrector time lagging
      false  // active
  );

  validateEquals(solve.getSolverNumber()       ,0);
  validateEquals(solve.getParentSolve()        ,exahype::solvers::Solve::InvalidParentSolveIdentifier);
  validateEquals(solve.getType()               ,exahype::solvers::Solve::SOLVE);
  validateEquals(solve.getTimeStepping()       ,exahype::solvers::Solve::GLOBAL);
  validateEquals(solve.isCorrectorTimeLagging(),true);
  validateEquals(solve.isActive()              ,false);

  validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getPredictorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getCorrectorTimeStamp(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getPredictorTimeStamp(),std::numeric_limits<double>::max(),1e-10);

  /*
   * check the bucket chain shifting of the time step sizes and stamps
   */
  // initialise predictor time stamp
  solve.setPredictorTimeStamp(0.);
  validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getPredictorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getCorrectorTimeStamp()   ,std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getPredictorTimeStamp()   ,0.                                ,1e-10);

  // step 1
  solve.updateNextPredictorTimeStepSize(5.);
  solve.startNewTimeStep();

  validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getCorrectorTimeStamp()   ,0.                                ,1e-10);

  validateNumericalEqualsWithEps(solve.getPredictorTimeStepSize(),5.                                ,1e-10);
  validateNumericalEqualsWithEps(solve.getPredictorTimeStamp()   ,5.                                ,1e-10);

  // step 2
  solve.updateNextPredictorTimeStepSize(20.);
  solve.startNewTimeStep();

  validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize(),5. ,1e-10);
  validateNumericalEqualsWithEps(solve.getCorrectorTimeStamp()   ,5. ,1e-10);

  validateNumericalEqualsWithEps(solve.getPredictorTimeStepSize(),20.,1e-10);
  validateNumericalEqualsWithEps(solve.getPredictorTimeStamp()   ,25.,1e-10);

  /*
   *  copy the solver
   */
  exahype::solvers::Solve solveCopy(solve);
  solveCopy.setParentSolve(0);
  solveCopy.setType(exahype::solvers::Solve::SUBSOLVE);

  // check everything again for new solver
  validateEquals(solveCopy.getSolverNumber()       ,0);
  validateEquals(solveCopy.getParentSolve()        ,0);
  validateEquals(solveCopy.getType()               ,exahype::solvers::Solve::SUBSOLVE);
  validateEquals(solveCopy.getTimeStepping()       ,exahype::solvers::Solve::GLOBAL);
  validateEquals(solveCopy.isCorrectorTimeLagging(),true);
  validateEquals(solveCopy.isActive()              ,false);

  validateNumericalEqualsWithEps(solveCopy.getCorrectorTimeStepSize(),5. ,1e-10);
  validateNumericalEqualsWithEps(solveCopy.getCorrectorTimeStamp()   ,5. ,1e-10);

  validateNumericalEqualsWithEps(solveCopy.getPredictorTimeStepSize(),20.,1e-10);
  validateNumericalEqualsWithEps(solveCopy.getPredictorTimeStamp()   ,25.,1e-10);

  /*
   * merge two solvers
   */
  // create an additional solver
  exahype::solvers::Solve otherSolve(
      0, // solverNumber
      exahype::solvers::Solve::InvalidParentSolveIdentifier,
      exahype::solvers::Solve::SOLVE,
      exahype::solvers::Solve::GLOBAL,
      true, // corrector time lagging
      true  // active
  );

  otherSolve.setCorrectorTimeStepSize(0.1);
  otherSolve.setCorrectorTimeStamp(0.2);

  otherSolve.setPredictorTimeStamp(0.3);
  otherSolve.setPredictorTimeStepSize(0.4);
  otherSolve.updateNextPredictorTimeStepSize(0.5);

  // merge two solvers
  solve.merge(otherSolve);

  // check properties of solve again
  validate(solve.isActive()==false);
  validateNumericalEqualsWithEps(solveCopy.getCorrectorTimeStepSize()    ,0.1,1e-10);
  validateNumericalEqualsWithEps(solveCopy.getCorrectorTimeStamp()       ,0.2,1e-10);

  validateNumericalEqualsWithEps(solveCopy.getPredictorTimeStepSize()    ,0.3,1e-10);
  validateNumericalEqualsWithEps(solveCopy.getPredictorTimeStamp()       ,0.4,1e-10);

  validateNumericalEqualsWithEps(solveCopy.getNextPredictorTimeStepSize(),0.5,1e-10);
}

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",on)
#endif
