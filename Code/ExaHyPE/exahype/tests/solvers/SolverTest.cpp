#include "exahype/tests/solvers/SolverTest.h"


#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"

#include "exahype/solvers/Solver.h"

#include <limits>

registerTest(exahype::tests::solvers::SolverTest)


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",off)
#endif


exahype::tests::solvers::SolverTest::SolverTest():
tarch::tests::TestCase ( "exahype::tests::solvers::SolverTest" ) {
}


exahype::tests::solvers::SolverTest::~SolverTest() {
}


void exahype::tests::solvers::SolverTest::run() {
  testMethod( testSolve );
}


void exahype::tests::solvers::SolverTest::testSolve() {
/*  exahype::solvers::Solve solve(
      0, // solverNumber
      true, // corrector time lagging
      false  // active
  );

  validateEquals(solve.getSolverNumber()       ,0);
  validateEquals(solve.isCorrectorTimeLagging(),true);
  validateEquals(solve.isActive()              ,false);

  validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getCorrectorTimeStamp(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getTimeStamp(),std::numeric_limits<double>::max(),1e-10);

  //
  // check the bucket chain shifting of the time step sizes and stamps
  //
  // initialise predictor time stamp
  solve.setTimeStamp(0.);
  validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getCorrectorTimeStamp()   ,std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getTimeStamp()   ,0.                                ,1e-10);

  // step 1
  solve.updateNextTimeStepSize(5.);
  solve.startNewTimeStep();

  validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getCorrectorTimeStamp()   ,0.                                ,1e-10);

  validateNumericalEqualsWithEps(solve.getTimeStepSize(),5.                                ,1e-10);
  validateNumericalEqualsWithEps(solve.getTimeStamp()   ,5.                                ,1e-10);

  // step 2
  solve.updateNextTimeStepSize(20.);
  solve.startNewTimeStep();

  validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize(),5. ,1e-10);
  validateNumericalEqualsWithEps(solve.getCorrectorTimeStamp()   ,5. ,1e-10);

  validateNumericalEqualsWithEps(solve.getTimeStepSize(),20.,1e-10);
  validateNumericalEqualsWithEps(solve.getTimeStamp()   ,25.,1e-10);

  //
  //  copy the solver
  //
  exahype::solvers::Solve solveCopy(solve);
  solveCopy.setParentSolve(0);

  // check everything again for new solver
  validateEquals(solveCopy.getSolverNumber()       ,0);
  validateEquals(solveCopy.getParentSolve()        ,0);
  validateEquals(solveCopy.isCorrectorTimeLagging(),true);
  validateEquals(solveCopy.isActive()              ,false);

  validateNumericalEqualsWithEps(solveCopy.getCorrectorTimeStepSize(),5. ,1e-10);
  validateNumericalEqualsWithEps(solveCopy.getCorrectorTimeStamp()   ,5. ,1e-10);

  validateNumericalEqualsWithEps(solveCopy.getTimeStepSize(),20.,1e-10);
  validateNumericalEqualsWithEps(solveCopy.getTimeStamp()   ,25.,1e-10);

  //merge two solvers
  //
  // create an additional solver
  exahype::solvers::Solve otherSolve(
      0, // solverNumber
      true, // corrector time lagging
      true  // active
  );

  otherSolve.setCorrectorTimeStamp(0.1);
  otherSolve.setCorrectorTimeStepSize(0.2);

  otherSolve.setTimeStamp(0.3);
  otherSolve.setTimeStepSize(0.4);
  otherSolve.updateNextTimeStepSize(0.5);

  // merge two solvers
  solve.merge(otherSolve);

  // check properties of solve again
  validate(solve.isActive()==false);
  validateNumericalEqualsWithEps(solve.getCorrectorTimeStamp()       ,0.1,1e-10);
  validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize()    ,0.2,1e-10);

  validateNumericalEqualsWithEps(solve.getTimeStamp()       ,0.3,1e-10);
  validateNumericalEqualsWithEps(solve.getTimeStepSize()    ,0.4,1e-10);

  validateNumericalEqualsWithEps(solve.getNextTimeStepSize(),0.5,1e-10);*/
}

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",on)
#endif
