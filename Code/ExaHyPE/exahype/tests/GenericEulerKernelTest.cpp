#include "exahype/tests/GenericEulerKernelTest.h"


#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"


#include "kernels/aderdg/generic/Kernels.h"



registerTest(exahype::tests::GenericEulerKernelTest)


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",off)
#endif

 
exahype::tests::GenericEulerKernelTest::GenericEulerKernelTest():
  tarch::tests::TestCase( "exahype::tests::GenericEulerKernelTest" ) {
}


exahype::tests::GenericEulerKernelTest::~GenericEulerKernelTest() {
}


void exahype::tests::GenericEulerKernelTest::testFlux(const double* const Q, double* f, double* g) {
  const double GAMMA = 1.4;

  const double irho = 1.0/Q[0];
  const double p = (GAMMA-1)*( Q[4] - 0.5* (Q[1]*Q[1] + Q[2]*Q[2]) * irho );

  f[0] = Q[1];
  f[1] = irho*Q[1]*Q[1] + p;
  f[2] = irho*Q[1]*Q[2];
  f[3] = irho*Q[1]*Q[3];
  f[4] = irho*Q[1]*(Q[4]+p);

  g[0] = Q[2];
  g[1] = irho*Q[2]*Q[1];
  g[2] = irho*Q[2]*Q[2] + p;
  g[3] = irho*Q[2]*Q[3];
  g[4] = irho*Q[2]*(Q[4]+p);
}


void exahype::tests::GenericEulerKernelTest::run() {
  // @todo If you have further tests, add them here
  testMethod( testSpaceTimePredictor2d );
}


void exahype::tests::GenericEulerKernelTest::testSpaceTimePredictor2d() {
  // @todo wieder rein
//  #if Dim2

  double *luh = new double[80]();  // space DOF
  for(int i=0; i<16; i++) {
    luh[5*i+0] = 1.00000000000000000000e+00;
    luh[5*i+4] = 2.50000000000000044409e+00;
  }
  // @todo Sollte kein Array sein, weil wir nur Wuerfel unterstuetzen
  const double dx[2]        = {3.70370370370370349811e-02, 3.70370370370370349811e-02};
  const double timeStepSize =  1.40831757919882352703e-03;

  // local:
  double *lQi = new double[320]; // space-time DOF
  double *lFi = new double[640];
  double *rhs0 = new double[320];
  double *rhs  = new double[320];
  double *tmp  = new double[20];

  // output:
  double *lQhi = new double[80];
  double *lFhi = new double[160];
  double *lQhbnd = new double[80];
  double *lFhbnd = new double[80];

  // Original aus Angelikas Code
  //
  //dg::spaceTimePredictor<2>(lQi, lFi, luh, lQhi, lFhi, lQhbnd, lFhbnd, rhs0, rhs, tmp, dx, timeStepSize);

  kernels::aderdg::generic::spaceTimePredictor<testFlux>( lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd, luh, dx[0], timeStepSize,
    5, // getNumberOfVariables(),
    4  // getNodesPerCoordinateAxis()
  );


  //validateNumericalEquals(<return vector>, <referencesolution>);


  // @todo Das +3 ist falsch, das hab ich nur reingebaut, um den Test abschmieren zu lassen
  validateNumericalEqualsWithEps(lQhi[0], 1+3, eps);

  // @todo wieder rein
//  #endif
}


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",on)
#endif
