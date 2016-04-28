#include "exahype/tests/kernels/c/GenericEulerKernelTest.h"

#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"

#include "kernels/aderdg/generic/Kernels.h"

using std::cout;
using std::endl;

registerTest(exahype::tests::c::GenericEulerKernelTest)
#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", off)
#endif

    exahype::tests::c::GenericEulerKernelTest::GenericEulerKernelTest()
    : tarch::tests::TestCase("exahype::tests::c::GenericEulerKernelTest") {
}

exahype::tests::c::GenericEulerKernelTest::~GenericEulerKernelTest() {}

void exahype::tests::c::GenericEulerKernelTest::run() {

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
