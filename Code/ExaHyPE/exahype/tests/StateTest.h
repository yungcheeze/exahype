#ifndef _EXAHYPE_TESTS_TEST_CASE_H_
#define _EXAHYPE_TESTS_TEST_CASE_H_

#include "tarch/tests/TestCase.h"

namespace exahype {
namespace tests {
class StateTest;
}
}

/**
 * This is just a default test case that demonstrated how to write unit tests
 * in Peano. Feel free to rename, remove, or duplicate it.
 */
class exahype::tests::StateTest : public tarch::tests::TestCase {
 private:
  /**
   * Tests the properties of a Solve
   * under various conditions.
   */
  void testState();

 public:
  StateTest();
  virtual ~StateTest();

  virtual void run();
};

#endif
