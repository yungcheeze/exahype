#include "exahype/tests/kernels/fortran/GenericEulerKernelTest.h"

#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"

#include "kernels/aderdg/generic/Kernels.h"

using std::cout;
using std::endl;

registerTest(exahype::tests::fortran::GenericEulerKernelTest)
#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", off)
#endif

    exahype::tests::fortran::GenericEulerKernelTest::GenericEulerKernelTest()
    : tarch::tests::TestCase("exahype::tests::fortran::GenericEulerKernelTest") {
}

exahype::tests::fortran::GenericEulerKernelTest::~GenericEulerKernelTest() {}

void exahype::tests::fortran::GenericEulerKernelTest::run() {

  testMethod(testPDEFluxes);

  testMethod(testSpaceTimePredictor);
  testMethod(testVolumeIntegral);
  testMethod(testRiemannSolver);
  testMethod(testSurfaceIntegral);

  testMethod(testSolutionUpdate);

}

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", on)
#endif
//#else
// todo VV TestCase
//#endif
