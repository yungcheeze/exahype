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
  // @todo If you have further tests, add them here
  testMethod( testSolve );
}


void exahype::tests::solvers::SolveTest::testSolve() {
  exahype::solvers::Solve solve(
                          0, // solverNumber
                          exahype::solvers::Solve::InvalidParentSolveIdentifier,
                          exahype::solvers::Solve::SOLVE,
                          exahype::solvers::Solve::GLOBAL,
                          false, // _useSameTimeStepSize
                          false  // active
  );

  validateEquals(solve.getSolverNumber()   ,0);
  validateEquals(solve.getParentSolve()    ,exahype::solvers::Solve::InvalidParentSolveIdentifier);
  validateEquals(solve.getType()           ,exahype::solvers::Solve::SOLVE);
  validateEquals(solve.getTimeStepping()   ,exahype::solvers::Solve::GLOBAL);
  validateEquals(solve.isSameTimeStepSize(),false);
  validateEquals(solve.isActive()          ,false);

  validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getPredictorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getCorrectorMinTimeStamp(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getPredictorMinTimeStamp(),std::numeric_limits<double>::max(),1e-10);

  // initialise predictor time stamp
 solve.setPredictorMinTimeStamp(0.);
 validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
 validateNumericalEqualsWithEps(solve.getPredictorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
 validateNumericalEqualsWithEps(solve.getCorrectorMinTimeStamp(),std::numeric_limits<double>::max(),1e-10);
 validateNumericalEqualsWithEps(solve.getPredictorMinTimeStamp(),0.                                ,1e-10);

  // check the bucket chain shifting of the time step sizes and stamps
  solve.updateNextPredictorTimeStepSize(5.);
  solve.startNewTimeStep();

  validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize(),std::numeric_limits<double>::max(),1e-10);
  validateNumericalEqualsWithEps(solve.getPredictorTimeStepSize(),5.                                ,1e-10);
  validateNumericalEqualsWithEps(solve.getCorrectorMinTimeStamp(),0.                                ,1e-10);
  validateNumericalEqualsWithEps(solve.getPredictorMinTimeStamp(),5.                                ,1e-10);

  solve.updateNextPredictorTimeStepSize(20.);
  solve.startNewTimeStep();

  validateNumericalEqualsWithEps(solve.getCorrectorTimeStepSize(),5. ,1e-10);
  validateNumericalEqualsWithEps(solve.getPredictorTimeStepSize(),20.,1e-10);
  validateNumericalEqualsWithEps(solve.getCorrectorMinTimeStamp(),5. ,1e-10);
  validateNumericalEqualsWithEps(solve.getPredictorMinTimeStamp(),25.,1e-10);

  // copy the solver
  exahype::solvers::Solve solve2(solve);
  solve2.setParentSolve(0);                          // todo create subsolve factory function
  solve2.setType(exahype::solvers::Solve::SUBSOLVE);

  validateEquals(solve2.getSolverNumber()   ,0);
  validateEquals(solve2.getParentSolve()    ,0);
  validateEquals(solve2.getType()           ,exahype::solvers::Solve::SUBSOLVE);
  validateEquals(solve2.getTimeStepping()   ,exahype::solvers::Solve::GLOBAL);
  validateEquals(solve2.isSameTimeStepSize(),false);
  validateEquals(solve2.isActive()          ,false);

  validateNumericalEqualsWithEps(solve2.getCorrectorTimeStepSize(),5. ,1e-10);
  validateNumericalEqualsWithEps(solve2.getPredictorTimeStepSize(),20.,1e-10);
  validateNumericalEqualsWithEps(solve2.getCorrectorMinTimeStamp(),5. ,1e-10);
  validateNumericalEqualsWithEps(solve2.getPredictorMinTimeStamp(),25.,1e-10);
}


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",on)
#endif
