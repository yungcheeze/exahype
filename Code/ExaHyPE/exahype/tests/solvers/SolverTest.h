#ifndef _EXAHYPE_TESTS_TEST_CASE_H_
#define _EXAHYPE_TESTS_TEST_CASE_H_

#include "tarch/tests/TestCase.h"

namespace exahype {
namespace tests {
namespace solvers {
class SolverTest;
}
}
}

/**
 * This is just a default test case that demonstrated how to write unit tests
 * in Peano. Feel free to rename, remove, or duplicate it.
 */
class exahype::tests::solvers::SolverTest : public tarch::tests::TestCase {
 private:
  /**
   * Tests the properties of a Solve
   * under various conditions.
   */
  void testSolve();

 public:
  SolverTest();
  virtual ~SolverTest();

  virtual void run();
};

#endif
