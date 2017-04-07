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

#include "exahype/tests/kernels/finitevolumes/riemannsolvers/c/RusanovTest.h"

#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"

#include "kernels/finitevolumes/riemannsolvers/c/riemannsolvers.h"

registerTest(exahype::tests::c::RusanovTest)

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", off)
#endif

tarch::logging::Log exahype::tests::c::RusanovTest::_log( "exahype::tests::c::RusanovTest" );

// ZeroFluxZeroNCP
exahype::tests::kernels::finitevolumes::riemannsolvers::c::ZeroFluxZeroNCP::ZeroFluxZeroNCP() {
  // do nothing
}

void
exahype::tests::kernels::finitevolumes::riemannsolvers::c::ZeroFluxZeroNCP::flux(const double* const Q, double** F) {
// do nothing
}

void
exahype::tests::kernels::finitevolumes::riemannsolvers::c::ZeroFluxZeroNCP::eigenvalues(
    const double* const Q,const int normalNonZeroIndex, double* lambda) {
  // always return 1.0
  std::fill_n(lambda, NumberOfVariables, 1.0);
}

bool
exahype::tests::kernels::finitevolumes::riemannsolvers::c::ZeroFluxZeroNCP::useNonConservativeProduct() {
  return false;
}

void exahype::tests::kernels::finitevolumes::riemannsolvers::c::ZeroFluxZeroNCP::nonConservativeProduct(
       const double* const Q, const double* const gradQ, double* BgradQ) {
  // do nothing
}


// RusanovTest
exahype::tests::kernels::finitevolumes::riemannsolvers::c::RusanovTest::RusanovTest()
    : tarch::tests::TestCase("exahype::tests::c::RusanovTest") {}

RusanovTest::~RusanovTest() {}

void RusanovTest::run() {
  _log.info("RusanovTest::run()", "RusanovTest is active");


}

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", on)
#endif
