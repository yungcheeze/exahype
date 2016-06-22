/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released unter the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#include "exahype/tests/kernels/c/GenericEulerKernelTest.h"

#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"

#include "kernels/aderdg/generic/Kernels.h"

using std::cout;
using std::endl;

// TODO: Do not conclude macro definitions with a semicolon?!
//       (https://goo.gl/22Ac4j)
// clang-format off
registerTest(exahype::tests::c::GenericEulerKernelTest)

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", off)
#endif

namespace exahype {
namespace tests {
namespace c {

GenericEulerKernelTest::GenericEulerKernelTest()
    : tarch::tests::TestCase("exahype::tests::c::GenericEulerKernelTest") {}

GenericEulerKernelTest::~GenericEulerKernelTest() {}

void GenericEulerKernelTest::run() {
  testMethod(testPDEFluxes);

  testMethod(testSpaceTimePredictorLinear);
  testMethod(testSpaceTimePredictorNonlinear);
  testMethod(testVolumeIntegralLinear);
  testMethod(testVolumeIntegralNonlinear);
  testMethod(testRiemannSolverLinear);
  testMethod(testRiemannSolverNonlinear);
  testMethod(testSurfaceIntegralLinear);
  testMethod(testSurfaceIntegralNonlinear);
  testMethod(testFaceUnknownsProjection);
  testMethod(testVolumeUnknownsProjection);

  testMethod(testSolutionUpdate);
}

}  // namespace c
}  // namespace tests
}  // namespace exahype

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", on)
#endif
//#else
// todo VV TestCase
//#endif

// clang-format on
