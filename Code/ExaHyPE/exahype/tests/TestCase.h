#define _EXAHYPE_TESTS_TEST_CASE_H_

#include "tarch/tests/TestCase.h"

namespace exahype {
namespace tests {
class TestCase;
}
}

/**
 * This is just a default test case that demonstrated how to write unit tests
 * in Peano. Feel free to rename, remove, or duplicate it.
 */
class exahype::tests::TestCase : public tarch::tests::TestCase {
 private:
  /**
   * These operation usually implement the real tests.
   */
  void test1();

  /**
   * These operation usually implement the real tests.
   */
  void test2();

  /**
   * These operation usually implement the real tests.
   */
  void test3();

 public:
  TestCase();
  virtual ~TestCase();

  virtual void run();
};

#endif
