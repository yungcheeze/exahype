// @todo ExaHyPE LIzenz
#ifndef _EXAHYPE_TESTS_GENERIC_EULER_KERNEL_TEST_H_
#define _EXAHYPE_TESTS_GENERIC_EULER_KERNEL_TEST_H_

#include "tarch/tests/TestCase.h"
#include "peano/utils/Globals.h"

namespace exahype {
namespace tests {
namespace fortran {
class GenericEulerKernelTest;
}
}
}

/**
 * This is just a default test case that demonstrated how to write unit tests
 * in Peano. Feel free to rename, remove, or duplicate it.
 */
class exahype::tests::fortran::GenericEulerKernelTest : public tarch::tests::TestCase {
 private:
 
  void testPDEFluxes();
  void testSpaceTimePredictor();
  void testVolumeIntegral();
  void testRiemannSolver();
  void testSurfaceIntegral();
  void testSolutionUpdate();

  const double eps = 1.0e-10;  // for quick adaption of the test cases (say,
                               // switch to single precision)

 public:
  GenericEulerKernelTest();
  virtual ~GenericEulerKernelTest();

  virtual void run();
};

#endif
