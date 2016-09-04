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

#include "exahype/tests/kernels/c/LimiterKernelTest.h"


registerTest(exahype::tests::c::LimiterKernelTest)

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", off)
#endif

tarch::logging::Log exahype::tests::c::LimiterKernelTest::_log( "exahype::tests::c::LimiterKernelTest" );

namespace exahype {
namespace tests {
namespace c {

const double LimiterKernelTest::eps = 1.0e-9;

LimiterKernelTest::LimiterKernelTest()
    : tarch::tests::TestCase("exahype::tests::c::LimiterKernelTest") {}

LimiterKernelTest::~LimiterKernelTest() {}

void LimiterKernelTest::run() {
  _log.info("LimiterKernelTest::run()", "LimiterKernelTest is active");

  
}

void LimiterKernelTest::testGetGaussLobattoData(const double* const luh, const int numberOfVariables, const int basisSize, const double* const expectedLob) {
  
}

void LimiterKernelTest::testGetFVMData(const double* const luh, const int numberOfVariables, const int basisSize, const int expectedbasisSizeLim, const double* const expectedLim) {
  
}
void LimiterKernelTest::testUpdateSubcellWithLimiterData(const double* const lim, const int numberOfVariables, const int basisSizeLim, const int basisSize, const double* const expectedLuh){
  
}

void LimiterKernelTest::testFindCellLocallocalMinlocalMax(const double* const luh, const int numberOfVariables, const int basisSize, const double* const expectedLocalMin, const double* const expectedLocalMax){
  
}

void LimiterKernelTest::testIsTroubledCell(const double* const luh, const int numberOfVariables, const int basisSize, const double* const troubledMin, const double* const troubledMax){
  
}

}  // namespace c
}  // namespace tests
}  // namespace exahype

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", on)
#endif


// clang-format on
