/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
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
