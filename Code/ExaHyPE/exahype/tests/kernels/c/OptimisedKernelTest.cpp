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

#ifdef TEST_OPT_KERNEL

#include <sstream>

#include "exahype/tests/kernels/c/OptimisedKernelTest.h"


registerTest(exahype::tests::c::OptimisedKernelTest)

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", off)
#endif

tarch::logging::Log exahype::tests::c::OptimisedKernelTest::_log( "exahype::tests::c::OptimisedKernelTest" );

namespace exahype {
namespace tests {
namespace c {

const double OptimisedKernelTest::eps = 1.0e-13;
#ifdef Dim2
const std::string OptimisedKernelTest::dim = "2";
#endif
#ifdef Dim3
const std::string OptimisedKernelTest::dim = "3";
#endif

OptimisedKernelTest::OptimisedKernelTest()
    : tarch::tests::TestCase("exahype::tests::c::OptimisedKernelTest") {
      
  _numberOfVariables = kernels::aderdg::optimised::getNumberOfVariable();
  _basisSize         = kernels::aderdg::optimised::getBasisSize();
  _order             = kernels::aderdg::optimised::getBasisSize() - 1;
}

OptimisedKernelTest::~OptimisedKernelTest() {}


void OptimisedKernelTest::adjustedSolutionValues(const double* const x,
                                                  const double w,
                                                  const double t,
                                                  const double dt, double* Q) {

  double GAMMA = 1.4;
  Q[0] = 1.;
  Q[1] = 0.;
  Q[2] = 0.;
  Q[3] = 0.;
  Q[4] = 1. / (GAMMA -1) +
        std::exp(-((x[0] -0.5) *(x[0] -0.5) + (x[1] -0.5) *(x[1] -0.5) 
#if DIMENSIONS == 3     
          + (x[2] -0.5) *(x[2] -0.5)
#endif        
        ) /
        (0.5 *0.5
#if DIMENSIONS == 3
          *0.5
#endif        
        )) *
        1.0e-1;
}

int OptimisedKernelTest::getNumberOfVariables() {
  return _numberOfVariables;
}

int OptimisedKernelTest::getNodesPerCoordinateAxis() {
  return _basisSize; //basisSize
}



void OptimisedKernelTest::run() {
  _log.info("OptimisedKernelTest::run()", "OptimisedKernelTest is active");
  testMethod(testSolutionAdjustment);
  

}


void OptimisedKernelTest::testSolutionAdjustment() {
  std::ostringstream out;
  out << "Test solutionAdjustment, ORDER="<< _order <<", NVAR=" << _numberOfVariables;
  logInfo("OptimisedKernelTest::testSolutionAdjustment()", out.str());
  
#if DIMENSIONS == 2     
  const int luhSize = _numberOfVariables*_basisSize*_basisSize;
  const double dx[2] = {1.0, 1.0};
  const double center[2] = {0.5, 0.5};
#else
  const int luhSize = _numberOfVariables*_basisSize*_basisSize*_basisSize;
  const double dx[2] = {1.0, 1.0, 1.0};
  const double center[2] = {0.5, 0.5, 0.5};
#endif 
  
  double* luh_generic = new double[luhSize]();
  double* luh_optimised = new double[luhSize]();

  double t = 0.0;
  double dt = 0.0;
  
  kernels::aderdg::generic::c::solutionAdjustment<OptimisedKernelTest>( *this, luh_generic, center[0], dx[0], t, dt );
  kernels::aderdg::optimised::solutionAdjustment<OptimisedKernelTest::adjustedSolutionValues>( luh_optimised, center[0], dx[0], t, dt );
  
  
  for(int i=0; i<luhSize; i++) {
    validateNumericalEqualsWithEps(luh_generic[i], luh_optimised[i], eps);
  }
  
  delete[] luh_generic;
  delete[] luh_optimised;
}

}  // namespace c
}  // namespace tests
}  // namespace exahype

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", on)
#endif

#endif //TEST_OPT_KERNEL