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
#include <iomanip> 
#include <random>
#include <cstring>

#include "exahype/tests/kernels/c/OptimisedKernelTest.h"


registerTest(exahype::tests::c::OptimisedKernelTest)

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", off)
#endif

tarch::logging::Log exahype::tests::c::OptimisedKernelTest::_log( "exahype::tests::c::OptimisedKernelTest" );

namespace exahype {
namespace tests {
namespace c {

int OptimisedKernelTest::_numberOfVariables;
int OptimisedKernelTest::_basisSize;
int OptimisedKernelTest::_order; 
 
const double OptimisedKernelTest::eps  = 1.0e-13;
const double OptimisedKernelTest::eps2 = 1.0e-12; //for known reordered operations
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
  Q[0] = 1. / (GAMMA -1) +
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
  
  // fill the rest with 
  for(int i=1; i<_numberOfVariables; i++) {
    Q[i] = i*x[0]*x[1]
#if DIMENSIONS == 3
          *x[2]
#endif       
      ;
  }
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
  testMethod(testSolutionUpdate);

}


void OptimisedKernelTest::testSolutionAdjustment() {
  std::ostringstream out;
  out << "Test solutionAdjustment with gaussian pulse on Q[0], ORDER="<< _order <<", NVAR=" << _numberOfVariables;
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


void OptimisedKernelTest::testSolutionUpdate() {
  std::ostringstream out;
  out << "Test testSolutionUpdate with random values, ORDER="<< _order <<", NVAR=" << _numberOfVariables;
  logInfo("OptimisedKernelTest::testSolutionUpdate()", out.str());
  
#if DIMENSIONS == 2     
  const int luhSize = _numberOfVariables*_basisSize*_basisSize;
#else
  const int luhSize = _numberOfVariables*_basisSize*_basisSize*_basisSize;
#endif 

  const double dt = 0.05;
  double* luh_generic = new double[luhSize];
  double* luh_optimised = new double[luhSize];
  double* lduh_generic = new double[luhSize];
  double* lduh_optimised = new double[luhSize];
  
  std::random_device rd; //to generate a randome seed
  std::mt19937 mt(rd()); //mersenne twister random number generator with random seed
  std::uniform_real_distribution<double> dist(-1.0, 1.0); // [-1.0,1.0)
  
  for(int i=0; i<luhSize; i++) {
    luh_generic[i] = dist(mt);
    lduh_generic[i] = dist(mt);
  }

  std::memcpy(luh_optimised, luh_generic, luhSize*sizeof(double));
  std::memcpy(lduh_optimised, lduh_generic, luhSize*sizeof(double));

  kernels::aderdg::generic::c::solutionUpdate( luh_generic, lduh_generic, dt, _numberOfVariables, 0, _basisSize );
  kernels::aderdg::optimised::solutionUpdate( luh_optimised, lduh_optimised, dt );
  
  for(int i=0; i<luhSize; i++) {
    validateNumericalEqualsWithEps(luh_generic[i], luh_optimised[i], eps2);
  }
  
  delete[] luh_generic;
  delete[] luh_optimised;
  delete[] lduh_generic;
  delete[] lduh_optimised;

}





}  // namespace c
}  // namespace tests
}  // namespace exahype

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", on)
#endif

#endif //TEST_OPT_KERNEL