#include "exahype/tests/kernels/GenericEulerKernelTest.h"


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

void exahype::tests::GenericEulerKernelTest::testEigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  const double GAMMA = 1.4;

  double irho = 1.0/Q[0];
  double p    = (GAMMA-1)*( Q[4] - 0.5* (Q[1]*Q[1] + Q[2]*Q[2]) * irho );

  double u_n = Q[normalNonZeroIndex+1] * irho;
  double c  = std::sqrt(GAMMA * p * irho);

  lambda[0] = u_n-c;
  lambda[1] = u_n;
  lambda[2] = u_n;
  lambda[3] = u_n;
  lambda[4] = u_n+c;
}

void exahype::tests::GenericEulerKernelTest::run() {
  // @todo If you have further tests, add them here
  testMethod( testPDEFluxes2d );

  testMethod( testVolumeIntegral2d );
  testMethod( testRiemannSolver2d );
  testMethod( testSurfaceIntegral2d );
  testMethod( testSolutionUpdate2d     );

  testMethod( testSpaceTimePredictor2d );
}

void exahype::tests::GenericEulerKernelTest::testPDEFluxes2d() {
  std::cout << "Test PDE-related functions" << std::endl;

  // 2D
  double Q[5] = {1., 0.1, 0.2, 0.3, 3.5}; // pressure = 1.39
  double f[5], g[5];

  // Angelika's old code
  //  problem::PDEFlux(Q,f,g);
  testFlux(Q,f,g);

  validateNumericalEquals(0.1, f[0]);
  validateNumericalEquals(1.4, f[1]);
  validateNumericalEqualsWithParams1(0.02, f[2], g[1]);
  validateNumericalEquals(0.03, f[3]);
  validateNumericalEquals(0.489, f[4]);

  validateNumericalEquals(0.2, g[0]);
  validateNumericalEquals(1.43, g[2]);
  validateNumericalEquals(0.06, g[3]);
  validateNumericalEquals(0.978,g[4]);

  //  // 3D
  //  double h[5];
  //  problem::(Q,f,g,h);
  //
  //  validateNumericalEquals(0.1, f[0]);
  //  validateNumericalEquals(1.4, f[1]);
  //  validateNumericalEqualsWithParams1(0.02, f[2], g[1]);
  //  validateNumericalEquals(0.03, f[3]);
  //  validateNumericalEquals(0.489, f[4]);
  //
  //  validateNumericalEquals(0.2, g[0]);
  //  validateNumericalEquals(1.43, g[2]);
  //  validateNumericalEquals(0.06, g[3]);
  //  validateNumericalEquals(0.978,g[4]);
  //
  //  validateNumericalEquals(0.3, h[0]);
  //  validateNumericalEquals(0.03,h[1]);
  //  validateNumericalEquals(0.06,h[2]);
  //  validateNumericalEquals(1.48,h[3]);
  //  validateNumericalEquals(1.4670,h[4]);
} // testPDEFluxes2d

void exahype::tests::GenericEulerKernelTest::testSolutionUpdate2d() {
  std::cout << "Test solution update, ORDER=3, DIM=2" << std::endl;

  // in/out:
  double * luh  = new double[80]();
  for(int i=0;i<80;i+=5) {
    luh[i]   = 1.0;
    luh[i+4] = 2.5;
  }

  // inputs:
  const double dt = 1.40831757919882352703e-03;
  double lduh[80] = {
      6.93136751845e-310,
      6.93136751845e-310,
      -3.85240924537,
      -1.22497305535e-16,
      5.34573771861,
      5.34573771861,
      2.85142157753,
      -7.22235167637,
      -7.74641761194e-17,
      -5.34573771861,
      5.34573771861,
      -2.85142157753,
      -7.22235167637,
      2.85795830992e-17,
      7.22235167637,
      2.85142157753,
      3.85240924537,
      -3.85240924537,
      2.93091460302e-17,
      -7.22235167637,
      8.69555536681e-322,
      8.69555536681e-322,
      2.85142157753,
      -3.85240924537,
      -3.85240924537,
      -1.22497305535e-16,
      5.34573771861,
      5.34573771861,
      2.85142157753,
      -7.22235167637,
      -7.74641761194e-17,
      -5.34573771861,
      5.34573771861,
      -2.85142157753,
      -7.22235167637,
      2.85795830992e-17,
      7.22235167637,
      2.85142157753,
      3.85240924537,
      -3.85240924537,
      2.93091460302e-17,
      -7.22235167637,
      -5.70284315505,
      3.33519014225e-319,
      2.05164022212e-16,
      -1.54928352239e-16,
      10.6914754372,
      -10.6914754372,
      0,
      -5.42249232835e-16,
      -6.4862093408e-17,
      -10.6914754372,
      -10.6914754372,
      0,
      -2.27017326925e-16,
      3.71888723946e-17,
      14.4447033527,
      -5.70284315505,
      0,
      1.3016105338e-16,
      3.915089398e-17,
      -7.70481849073,
      7.70481849073,
      0,
      1.3702812893e-16,
      5.71591661983e-17,
      5.70284315505,
      14.4447033527,
      0,
      2.00057081693e-16,
      3.71888723946e-17,
      -5.70284315505,
      14.4447033527,
      0,
      1.3016105338e-16,
      5.14845862726e-17,
      7.70481849073,
      7.70481849073,
      0,
      1.80196051954e-16};

  // Angelika's old code.
  // dg::updateSolution<2>(luh,lduh,dt);
  kernels::aderdg::generic::solutionUpdate( luh, lduh, dt,
                                            5,   //getNumberOfVariables(),
                                            4   // getNodesPerCoordinateAxis()
  );

  validateNumericalEqualsWithEps(luh[0], 1, eps);
  validateNumericalEqualsWithEps(luh[1], 3.22688437998e-311, eps);
  validateNumericalEqualsWithEps(luh[2], -0.179348147189, eps);
  validateNumericalEqualsWithEps(luh[3], -5.70283772674e-18, eps);
  validateNumericalEqualsWithEps(luh[4], 2.74886975763, eps);
  validateNumericalEqualsWithEps(luh[5], 1.13274736515, eps);
  validateNumericalEqualsWithEps(luh[6], 0.0708075706777, eps);
  validateNumericalEqualsWithEps(luh[7], -0.179348147189, eps);
  validateNumericalEqualsWithEps(luh[8], -1.92361949169e-18, eps);
  validateNumericalEqualsWithEps(luh[9], 2.36725263485, eps);
  validateNumericalEqualsWithEps(luh[10], 1.13274736515, eps);
  validateNumericalEqualsWithEps(luh[11], -0.0708075706777, eps);
  validateNumericalEqualsWithEps(luh[12], -0.179348147189, eps);
  validateNumericalEqualsWithEps(luh[13], 7.09698932695e-19, eps);
  validateNumericalEqualsWithEps(luh[14], 2.67934814719, eps);
  validateNumericalEqualsWithEps(luh[15], 1.13274736515, eps);
  validateNumericalEqualsWithEps(luh[16], 0.179348147189, eps);
  validateNumericalEqualsWithEps(luh[17], -0.179348147189, eps);
  validateNumericalEqualsWithEps(luh[18], 1.36448147157e-18, eps);
  validateNumericalEqualsWithEps(luh[19], 2.16376485233, eps);
  validateNumericalEqualsWithEps(luh[20], 1, eps);
  validateNumericalEqualsWithEps(luh[21], 0, eps);
  validateNumericalEqualsWithEps(luh[22], 0.0708075706777, eps);
  validateNumericalEqualsWithEps(luh[23], -0.0956644720902, eps);
  validateNumericalEqualsWithEps(luh[24], 2.40433552791, eps);
  validateNumericalEqualsWithEps(luh[25], 1, eps);
  validateNumericalEqualsWithEps(luh[26], 0.0708075706776, eps);
  validateNumericalEqualsWithEps(luh[27], 0.0708075706776, eps);
  validateNumericalEqualsWithEps(luh[28], 0.0377688254662, eps);
  validateNumericalEqualsWithEps(luh[29], 2.40433552791, eps);
  validateNumericalEqualsWithEps(luh[30], 1, eps);
  validateNumericalEqualsWithEps(luh[31], -0.0708075706776, eps);
  validateNumericalEqualsWithEps(luh[32], 0.0708075706776, eps);
  validateNumericalEqualsWithEps(luh[33], -0.0377688254662, eps);
  validateNumericalEqualsWithEps(luh[34], 2.40433552791, eps);
  validateNumericalEqualsWithEps(luh[35], 1, eps);
  validateNumericalEqualsWithEps(luh[36], 0.179348147189, eps);
  validateNumericalEqualsWithEps(luh[37], 0.0708075706777, eps);
  validateNumericalEqualsWithEps(luh[38], 0.0956644720902, eps);
  validateNumericalEqualsWithEps(luh[39], 2.40433552791, eps);
  validateNumericalEqualsWithEps(luh[40], 1, eps);
  validateNumericalEqualsWithEps(luh[41], -0.179348147189, eps);
  validateNumericalEqualsWithEps(luh[42], -0.141615141355, eps);
  validateNumericalEqualsWithEps(luh[43], 8.27559956784e-321, eps);
  validateNumericalEqualsWithEps(luh[44], 2.5, eps);
  validateNumericalEqualsWithEps(luh[45], 1, eps);
  validateNumericalEqualsWithEps(luh[46], 0.141615141355, eps);
  validateNumericalEqualsWithEps(luh[47], -0.141615141355, eps);
  validateNumericalEqualsWithEps(luh[48], 0, eps);
  validateNumericalEqualsWithEps(luh[49], 2.5, eps);
  validateNumericalEqualsWithEps(luh[50], 1, eps);
  validateNumericalEqualsWithEps(luh[51], -0.141615141355, eps);
  validateNumericalEqualsWithEps(luh[52], -0.141615141355, eps);
  validateNumericalEqualsWithEps(luh[53], 0, eps);
  validateNumericalEqualsWithEps(luh[54], 2.5, eps);
  validateNumericalEqualsWithEps(luh[55], 1, eps);
  validateNumericalEqualsWithEps(luh[56], 0.358696294376, eps);
  validateNumericalEqualsWithEps(luh[57], -0.141615141355, eps);
  validateNumericalEqualsWithEps(luh[58], 0, eps);
  validateNumericalEqualsWithEps(luh[59], 2.5, eps);
  validateNumericalEqualsWithEps(luh[60], 1, eps);
  validateNumericalEqualsWithEps(luh[61], -0.358696294377, eps);
  validateNumericalEqualsWithEps(luh[62], 0.358696294377, eps);
  validateNumericalEqualsWithEps(luh[63], 0, eps);
  validateNumericalEqualsWithEps(luh[64], 2.5, eps);
  validateNumericalEqualsWithEps(luh[65], 1, eps);
  validateNumericalEqualsWithEps(luh[66], 0.141615141355, eps);
  validateNumericalEqualsWithEps(luh[67], 0.358696294376, eps);
  validateNumericalEqualsWithEps(luh[68], 0, eps);
  validateNumericalEqualsWithEps(luh[69], 2.5, eps);
  validateNumericalEqualsWithEps(luh[70], 1, eps);
  validateNumericalEqualsWithEps(luh[71], -0.141615141355, eps);
  validateNumericalEqualsWithEps(luh[72], 0.358696294376, eps);
  validateNumericalEqualsWithEps(luh[73], 0, eps);
  validateNumericalEqualsWithEps(luh[74], 2.5, eps);
  validateNumericalEqualsWithEps(luh[75], 1, eps);
  validateNumericalEqualsWithEps(luh[76], 0.358696294377, eps);
  validateNumericalEqualsWithEps(luh[77], 0.358696294377, eps);
  validateNumericalEqualsWithEps(luh[78], 0, eps);
  validateNumericalEqualsWithEps(luh[79], 2.5, eps);

  delete[] luh;
} // testSolutionUpdate2d

void exahype::tests::GenericEulerKernelTest::testSurfaceIntegral2d() {
  std::cout << "Test surface integral, ORDER=3, DIM=2" << std::endl;

  // inputs:
  const double dx[2] = {0.1, 0.1};  // mesh spacing
  double lFhbnd[4*20] = {0.};

  double* FLeft  = &lFhbnd[ 0];
  double* FRight = &lFhbnd[20];
  double* FFront = &lFhbnd[40];
  double* FBack  = &lFhbnd[60];

  //  Angelika's original code
  //  double FLeft[20] = {0.};
  //  double FRight[20] = {0.};
  //  double FFront[20] = {0.};
  //  double FBack[20] = {0.};
  for(int i=0;i<20; i+=5) {
    // in x orientation 1
    FLeft[ i+1] = 1.;
    FRight[i+1] = 1.;
    // in y orientation 1
    FFront[i+2] = 1.;
    FBack[ i+2] = 1.;
  }

  // in/out:
  double lduh[80] {
    2.68172016875e-17,
    -7.70481849073,
    -7.70481849073,
    0,
    9.38602059063e-17,
    7.85885858643e-17,
    5.70284315505,
    -14.4447033527,
    0,
    2.75060050525e-16,
    5.86182920605e-17,
    -5.70284315505,
    -14.4447033527,
    0,
    2.05164022212e-16,
    3.915089398e-17,
    7.70481849073,
    -7.70481849073,
    0,
    1.3702812893e-16,
    7.85885858643e-17,
    -14.4447033527,
    5.70284315505,
    0,
    2.75060050525e-16,
    -2.4499461107e-16,
    10.6914754372,
    10.6914754372,
    0,
    -8.57481138746e-16,
    -1.54928352239e-16,
    -10.6914754372,
    10.6914754372,
    0,
    -5.42249232835e-16,
    5.71591661981e-17,
    14.4447033527,
    5.70284315505,
    0,
    2.00057081693e-16,
    5.86182920605e-17,
    -14.4447033527,
    -5.70284315505,
    0,
    2.05164022212e-16,
    -1.54928352239e-16,
    10.6914754372,
    -10.6914754372,
    0,
    -5.42249232835e-16,
    -6.48620934071e-17,
    -10.6914754372,
    -10.6914754372,
    0,
    -2.27017326925e-16,
    3.71888723943e-17,
    14.4447033527,
    -5.70284315505,
    0,
    1.3016105338e-16,
    3.915089398e-17,
    -7.70481849073,
    7.70481849073,
    0,
    1.3702812893e-16,
    5.71591661981e-17,
    5.70284315505,
    14.4447033527,
    0,
    2.00057081693e-16,
    3.71888723943e-17,
    -5.70284315505,
    14.4447033527,
    0,
    1.3016105338e-16,
    5.14845862726e-17,
    7.70481849073,
    7.70481849073,
    0,
    1.80196051954e-16
  };

  // Angelika's original code
  //    dg::surfaceIntegral(lduh, dx, FLeft, FRight, FFront, FBack);

//    kernels::aderdg::generic::surfaceIntegral( lduh, lFhbnd, dx[0],
//                                               5,   //getNumberOfVariables(),
//                                               4   // getNodesPerCoordinateAxis()
//    );

  kernels::aderdg::generic::surfaceIntegral2( lduh, lFhbnd, dx[0],
                                              5,   //getNumberOfVariables(),
                                              4   // getNodesPerCoordinateAxis()
  );

  validateNumericalEqualsWithEps(lduh[0],3.40976703141e-17, eps);
  validateNumericalEqualsWithEps(lduh[1],-4.85118201268, eps);
  validateNumericalEqualsWithEps(lduh[2],-4.85118201268, eps);
  validateNumericalEqualsWithEps(lduh[3],0, eps);
  validateNumericalEqualsWithEps(lduh[4],1.19341846099e-16, eps);
  validateNumericalEqualsWithEps(lduh[5],-2.88746795166e-16, eps);
  validateNumericalEqualsWithEps(lduh[6],3.59067902355, eps);
  validateNumericalEqualsWithEps(lduh[7],-9.0948132221, eps);
  validateNumericalEqualsWithEps(lduh[8],0, eps);
  validateNumericalEqualsWithEps(lduh[9],2.19250415768e-16, eps);
  validateNumericalEqualsWithEps(lduh[10],-3.03328335072e-16, eps);
  validateNumericalEqualsWithEps(lduh[11],-3.59067902355, eps);
  validateNumericalEqualsWithEps(lduh[12],-9.0948132221, eps);
  validateNumericalEqualsWithEps(lduh[13],0, eps);
  validateNumericalEqualsWithEps(lduh[14],1.68215026097e-16, eps);
  validateNumericalEqualsWithEps(lduh[15],3.915089398e-17, eps);
  validateNumericalEqualsWithEps(lduh[16],4.85118201268, eps);
  validateNumericalEqualsWithEps(lduh[17],-4.85118201268, eps);
  validateNumericalEqualsWithEps(lduh[18],0, eps);
  validateNumericalEqualsWithEps(lduh[19],1.3702812893e-16, eps);
  validateNumericalEqualsWithEps(lduh[20],-2.88746795166e-16, eps);
  validateNumericalEqualsWithEps(lduh[21],-9.0948132221, eps);
  validateNumericalEqualsWithEps(lduh[22],3.59067902355, eps);
  validateNumericalEqualsWithEps(lduh[23],0, eps);
  validateNumericalEqualsWithEps(lduh[24],2.19250415768e-16, eps);
  validateNumericalEqualsWithEps(lduh[25],2.94794991158e-16, eps);
  validateNumericalEqualsWithEps(lduh[26],6.73166971973, eps);
  validateNumericalEqualsWithEps(lduh[27],6.73166971973, eps);
  validateNumericalEqualsWithEps(lduh[28],0, eps);
  validateNumericalEqualsWithEps(lduh[29],-7.88824372364e-16, eps);
  validateNumericalEqualsWithEps(lduh[30],-1.54928352239e-16, eps);
  validateNumericalEqualsWithEps(lduh[31],-6.73166971973, eps);
  validateNumericalEqualsWithEps(lduh[32],6.73166971973, eps);
  validateNumericalEqualsWithEps(lduh[33],0, eps);
  validateNumericalEqualsWithEps(lduh[34],-5.42249232835e-16, eps);
  validateNumericalEqualsWithEps(lduh[35],4.19105793331e-16, eps);
  validateNumericalEqualsWithEps(lduh[36],9.0948132221, eps);
  validateNumericalEqualsWithEps(lduh[37],3.59067902355, eps);
  validateNumericalEqualsWithEps(lduh[38],0, eps);
  validateNumericalEqualsWithEps(lduh[39],2.37006077808e-16, eps);
  validateNumericalEqualsWithEps(lduh[40],-3.03328335072e-16, eps);
  validateNumericalEqualsWithEps(lduh[41],-9.0948132221, eps);
  validateNumericalEqualsWithEps(lduh[42],-3.59067902355, eps);
  validateNumericalEqualsWithEps(lduh[43],0, eps);
  validateNumericalEqualsWithEps(lduh[44],1.68215026097e-16, eps);
  validateNumericalEqualsWithEps(lduh[45],-1.54928352239e-16, eps);
  validateNumericalEqualsWithEps(lduh[46],6.73166971973, eps);
  validateNumericalEqualsWithEps(lduh[47],-6.73166971973, eps);
  validateNumericalEqualsWithEps(lduh[48],0, eps);
  validateNumericalEqualsWithEps(lduh[49],-5.42249232835e-16, eps);
  validateNumericalEqualsWithEps(lduh[50],-6.04651695635e-16, eps);
  validateNumericalEqualsWithEps(lduh[51],-6.73166971973, eps);
  validateNumericalEqualsWithEps(lduh[52],-6.73166971973, eps);
  validateNumericalEqualsWithEps(lduh[53],0, eps);
  validateNumericalEqualsWithEps(lduh[54],-2.95674093307e-16, eps);
  validateNumericalEqualsWithEps(lduh[55],4.04524253424e-16, eps);
  validateNumericalEqualsWithEps(lduh[56],9.0948132221, eps);
  validateNumericalEqualsWithEps(lduh[57],-3.59067902355, eps);
  validateNumericalEqualsWithEps(lduh[58],0, eps);
  validateNumericalEqualsWithEps(lduh[59],1.85970688137e-16, eps);
  validateNumericalEqualsWithEps(lduh[60],3.915089398e-17, eps);
  validateNumericalEqualsWithEps(lduh[61],-4.85118201268, eps);
  validateNumericalEqualsWithEps(lduh[62],4.85118201268, eps);
  validateNumericalEqualsWithEps(lduh[63],0, eps);
  validateNumericalEqualsWithEps(lduh[64],1.3702812893e-16, eps);
  validateNumericalEqualsWithEps(lduh[65],4.19105793331e-16, eps);
  validateNumericalEqualsWithEps(lduh[66],3.59067902355, eps);
  validateNumericalEqualsWithEps(lduh[67],9.0948132221, eps);
  validateNumericalEqualsWithEps(lduh[68],0, eps);
  validateNumericalEqualsWithEps(lduh[69],2.37006077808e-16, eps);
  validateNumericalEqualsWithEps(lduh[70],4.04524253424e-16, eps);
  validateNumericalEqualsWithEps(lduh[71],-3.59067902355, eps);
  validateNumericalEqualsWithEps(lduh[72],9.0948132221, eps);
  validateNumericalEqualsWithEps(lduh[73],0, eps);
  validateNumericalEqualsWithEps(lduh[74],1.85970688137e-16, eps);
  validateNumericalEqualsWithEps(lduh[75],4.4204117646e-17, eps);
  validateNumericalEqualsWithEps(lduh[76],4.85118201268, eps);
  validateNumericalEqualsWithEps(lduh[77],4.85118201268, eps);
  validateNumericalEqualsWithEps(lduh[78],0, eps);
  validateNumericalEqualsWithEps(lduh[79],1.54714411761e-16, eps);
} // testSurfaceIntegral2d

void exahype::tests::GenericEulerKernelTest::testRiemannSolver2d() {
  // Rusanov
  std::cout << "Test Riemann Solver (Rusanov), ORDER=3, DIM=2" << std::endl;
  // input:
  double QL[20] = {1.,0.,0.,0.,2.5,1.,0.,0.,0.,2.5,1.,0.,0.,0.,2.5,1.,0.,0.,0.,2.5};
  double QR[20] = {1.,0.,0.,0.,2.5,1.,0.,0.,0.,2.5,1.,0.,0.,0.,2.5,1.,0.,0.,0.,2.5};

  //    const double nx[3]= { 1., 0., 0. }; // normal vector in x direction

  // 04/02/16:Dominic Etienne Charrier
  // Angelika's original code
  // Temporary arrays are now allocated by the kernel function.
  //  // local:
  //  double QavL[5]    __attribute__((aligned(16/*ALIGNMENT*/)));
  //  double QavR[5]    __attribute__((aligned(16)));
  //  double lambdaL[5] __attribute__((aligned(16)));
  //  double lambdaR[5] __attribute__((aligned(16)));

  // output:
  double *FL = new double[20];
  double *FR = new double[20];

  // Angelika;s original code
  // dg::solveRiemannProblem<2>(FL, FR, QL, QR, QavL, QavR, lambdaL, lambdaR, 0.0 /*unused*/, 0.0 /*unused*/, nx);

  kernels::aderdg::generic::riemannSolver<testEigenvalues>( FL, FR, QL, QR,
                                                            0.0, // dt
                                                            0,   // normalNonZero
                                                            5,   //getNumberOfVariables(),
                                                            4   // getNodesPerCoordinateAxis()
  );

  // (a) FL == FR, element by element
  for(int i=0;i<20;i++) {
    validateEquals(FL[i],FR[i]);
  }

  // 04/02/16:Dominic Etienne Charrier
  // Angelika's original code
  // Moved this check as a comment inside of the kernel.
  // (b) check max speed
  //  double sMax = 0;
  //  for(int ivar=0; ivar < 5; ivar++) {
  //    sMax = std::max(sMax, std::max(fabs(lambdaL[ivar]),fabs(lambdaR[ivar])));
  //  }
  //  validateNumericalEqualsWithEps(sMax, 1.18321595661992, eps);

  delete[] FL;
  delete[] FR;
} // testRiemannSolver

void exahype::tests::GenericEulerKernelTest::testVolumeIntegral2d() {
  std::cout << "Test volume integral, ORDER=3, DIM=2" << std::endl;

  // output:
  double *lduh = new double[80];

  // input:
  const double dx[2] = {3.70370370370370349811e-02, 3.70370370370370349811e-02}; // mesh spacing
  const double lFhi[160] = {
      -5.78725778411e-18,
      1,
      3.90655760696e-35,
      0,
      -2.02554022444e-17,
      3.32798659798e-18,
      1,
      -6.28505465062e-35,
      0,
      1.16479530929e-17,
      -3.07890876627e-18,
      1,
      6.02226203114e-35,
      0,
      -1.07761806819e-17,
      7.66579600935e-18,
      1,
      -5.27718213855e-35,
      0,
      2.68302860327e-17,
      -1.48132747725e-17,
      1,
      -6.28505465062e-35,
      0,
      -5.18464617038e-17,
      7.63123721751e-18,
      1,
      1.09880295341e-34,
      0,
      2.67093302613e-17,
      -4.17981679914e-18,
      1,
      -7.10324058069e-35,
      0,
      -1.4629358797e-17,
      1.23156350651e-17,
      1,
      5.00229224796e-35,
      0,
      4.31047227278e-17,
      -1.48132747725e-17,
      1,
      6.02226203114e-35,
      0,
      -5.18464617038e-17,
      7.63123721751e-18,
      1,
      -7.10324058069e-35,
      0,
      2.67093302613e-17,
      -4.17981679914e-18,
      1,
      4.8794232613e-35,
      0,
      -1.4629358797e-17,
      1.23156350651e-17,
      1,
      -4.73179353928e-35,
      0,
      4.31047227278e-17,
      -5.78725778411e-18,
      1,
      -5.27718213855e-35,
      0,
      -2.02554022444e-17,
      3.32798659798e-18,
      1,
      5.00229224796e-35,
      0,
      1.16479530929e-17,
      -3.07890876627e-18,
      1,
      -4.73179353928e-35,
      0,
      -1.07761806819e-17,
      7.66579600935e-18,
      1,
      7.37687416339e-35,
      0,
      2.68302860327e-17,
      -5.78725778411e-18,
      3.90655760696e-35,
      1,
      0,
      -2.02554022444e-17,
      -1.48132747725e-17,
      -6.28505465062e-35,
      1,
      0,
      -5.18464617038e-17,
      -1.48132747725e-17,
      6.02226203114e-35,
      1,
      0,
      -5.18464617038e-17,
      -5.78725778411e-18,
      -5.27718213855e-35,
      1,
      0,
      -2.02554022444e-17,
      3.32798659798e-18,
      -6.28505465062e-35,
      1,
      0,
      1.16479530929e-17,
      7.63123721751e-18,
      1.09880295341e-34,
      1,
      0,
      2.67093302613e-17,
      7.63123721751e-18,
      -7.10324058069e-35,
      1,
      0,
      2.67093302613e-17,
      3.32798659798e-18,
      5.00229224796e-35,
      1,
      0,
      1.16479530929e-17,
      -3.07890876627e-18,
      6.02226203114e-35,
      1,
      0,
      -1.07761806819e-17,
      -4.17981679914e-18,
      -7.10324058069e-35,
      1,
      0,
      -1.4629358797e-17,
      -4.17981679914e-18,
      4.8794232613e-35,
      1,
      0,
      -1.4629358797e-17,
      -3.07890876627e-18,
      -4.73179353928e-35,
      1,
      0,
      -1.07761806819e-17,
      7.66579600935e-18,
      -5.27718213855e-35,
      1,
      0,
      2.68302860327e-17,
      1.23156350651e-17,
      5.00229224796e-35,
      1,
      0,
      4.31047227278e-17,
      1.23156350651e-17,
      -4.73179353928e-35,
      1,
      0,
      4.31047227278e-17,
      7.66579600935e-18,
      7.37687416339e-35,
      1,
      0,
      2.68302860327e-17
  };

  kernels::aderdg::generic::volumeIntegral( lduh, lFhi,dx[0],
                                            5, // getNumberOfVariables(),
                                            4 //getNodesPerCoordinateAxis()
  );

  validateNumericalEqualsWithEps(lduh[0], 2.68172016875e-17, eps);
  validateNumericalEqualsWithEps(lduh[1], -7.70481849073, eps);
  validateNumericalEqualsWithEps(lduh[2], -7.70481849073, eps);
  validateNumericalEqualsWithEps(lduh[3], 0, eps);
  validateNumericalEqualsWithEps(lduh[4], 9.38602059063e-17, eps);
  validateNumericalEqualsWithEps(lduh[5], 7.85885858643e-17, eps);
  validateNumericalEqualsWithEps(lduh[6], 5.70284315505, eps);
  validateNumericalEqualsWithEps(lduh[7], -14.4447033527, eps);
  validateNumericalEqualsWithEps(lduh[8], 0, eps);
  validateNumericalEqualsWithEps(lduh[9], 2.75060050525e-16, eps);
  validateNumericalEqualsWithEps(lduh[10], 5.86182920605e-17, eps);
  validateNumericalEqualsWithEps(lduh[11], -5.70284315505, eps);
  validateNumericalEqualsWithEps(lduh[12], -14.4447033527, eps);
  validateNumericalEqualsWithEps(lduh[13], 0, eps);
  validateNumericalEqualsWithEps(lduh[14], 2.05164022212e-16, eps);
  validateNumericalEqualsWithEps(lduh[15], 3.915089398e-17, eps);
  validateNumericalEqualsWithEps(lduh[16], 7.70481849073, eps);
  validateNumericalEqualsWithEps(lduh[17], -7.70481849073, eps);
  validateNumericalEqualsWithEps(lduh[18], 0, eps);
  validateNumericalEqualsWithEps(lduh[19], 1.3702812893e-16, eps);
  validateNumericalEqualsWithEps(lduh[20], 7.85885858643e-17, eps);
  validateNumericalEqualsWithEps(lduh[21], -14.4447033527, eps);
  validateNumericalEqualsWithEps(lduh[22], 5.70284315505, eps);
  validateNumericalEqualsWithEps(lduh[23], 0, eps);
  validateNumericalEqualsWithEps(lduh[24], 2.75060050525e-16, eps);
  validateNumericalEqualsWithEps(lduh[25], -2.4499461107e-16, eps);
  validateNumericalEqualsWithEps(lduh[26], 10.6914754372, eps);
  validateNumericalEqualsWithEps(lduh[27], 10.6914754372, eps);
  validateNumericalEqualsWithEps(lduh[28], 0, eps);
  validateNumericalEqualsWithEps(lduh[29], -8.57481138746e-16, eps);
  validateNumericalEqualsWithEps(lduh[30], -1.54928352239e-16, eps);
  validateNumericalEqualsWithEps(lduh[31], -10.6914754372, eps);
  validateNumericalEqualsWithEps(lduh[32], 10.6914754372, eps);
  validateNumericalEqualsWithEps(lduh[33], 0, eps);
  validateNumericalEqualsWithEps(lduh[34], -5.42249232835e-16, eps);
  validateNumericalEqualsWithEps(lduh[35], 5.71591661981e-17, eps);
  validateNumericalEqualsWithEps(lduh[36], 14.4447033527, eps);
  validateNumericalEqualsWithEps(lduh[37], 5.70284315505, eps);
  validateNumericalEqualsWithEps(lduh[38], 0, eps);
  validateNumericalEqualsWithEps(lduh[39], 2.00057081693e-16, eps);
  validateNumericalEqualsWithEps(lduh[40], 5.86182920605e-17, eps);
  validateNumericalEqualsWithEps(lduh[41], -14.4447033527, eps);
  validateNumericalEqualsWithEps(lduh[42], -5.70284315505, eps);
  validateNumericalEqualsWithEps(lduh[43], 0, eps);
  validateNumericalEqualsWithEps(lduh[44], 2.05164022212e-16, eps);
  validateNumericalEqualsWithEps(lduh[45], -1.54928352239e-16, eps);
  validateNumericalEqualsWithEps(lduh[46], 10.6914754372, eps);
  validateNumericalEqualsWithEps(lduh[47], -10.6914754372, eps);
  validateNumericalEqualsWithEps(lduh[48], 0, eps);
  validateNumericalEqualsWithEps(lduh[49], -5.42249232835e-16, eps);
  validateNumericalEqualsWithEps(lduh[50], -6.48620934071e-17, eps);
  validateNumericalEqualsWithEps(lduh[51], -10.6914754372, eps);
  validateNumericalEqualsWithEps(lduh[52], -10.6914754372, eps);
  validateNumericalEqualsWithEps(lduh[53], 0, eps);
  validateNumericalEqualsWithEps(lduh[54], -2.27017326925e-16, eps);
  validateNumericalEqualsWithEps(lduh[55], 3.71888723943e-17, eps);
  validateNumericalEqualsWithEps(lduh[56], 14.4447033527, eps);
  validateNumericalEqualsWithEps(lduh[57], -5.70284315505, eps);
  validateNumericalEqualsWithEps(lduh[58], 0, eps);
  validateNumericalEqualsWithEps(lduh[59], 1.3016105338e-16, eps);
  validateNumericalEqualsWithEps(lduh[60], 3.915089398e-17, eps);
  validateNumericalEqualsWithEps(lduh[61], -7.70481849073, eps);
  validateNumericalEqualsWithEps(lduh[62], 7.70481849073, eps);
  validateNumericalEqualsWithEps(lduh[63], 0, eps);
  validateNumericalEqualsWithEps(lduh[64], 1.3702812893e-16, eps);
  validateNumericalEqualsWithEps(lduh[65], 5.71591661981e-17, eps);
  validateNumericalEqualsWithEps(lduh[66], 5.70284315505, eps);
  validateNumericalEqualsWithEps(lduh[67], 14.4447033527, eps);
  validateNumericalEqualsWithEps(lduh[68], 0, eps);
  validateNumericalEqualsWithEps(lduh[69], 2.00057081693e-16, eps);
  validateNumericalEqualsWithEps(lduh[70], 3.71888723943e-17, eps);
  validateNumericalEqualsWithEps(lduh[71], -5.70284315505, eps);
  validateNumericalEqualsWithEps(lduh[72], 14.4447033527, eps);
  validateNumericalEqualsWithEps(lduh[73], 0, eps);
  validateNumericalEqualsWithEps(lduh[74], 1.3016105338e-16, eps);
  validateNumericalEqualsWithEps(lduh[75], 5.14845862726e-17, eps);
  validateNumericalEqualsWithEps(lduh[76], 7.70481849073, eps);
  validateNumericalEqualsWithEps(lduh[77], 7.70481849073, eps);
  validateNumericalEqualsWithEps(lduh[78], 0, eps);
  validateNumericalEqualsWithEps(lduh[79], 1.80196051954e-16, eps);

  delete[] lduh;
} // testVolumeIntegral2d

void exahype::tests::GenericEulerKernelTest::testSpaceTimePredictor2d() {
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
  //dg::spaceTimePredictor<2>(lQi, lFi, luh, lQhi, lFhi, lQhbnd, lFhbnd, rhs0, rhs, tmp, dx, timeStepSize);

  // todo 10/02/16:Dominic Etienne Charrier
  // REMOVED
  /*
  kernels::aderdg::generic::spaceTimePredictor<testFlux>( lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd, luh, dx[0], timeStepSize,
                                                          5, // getNumberOfVariables(),
                                                          4  // getNodesPerCoordinateAxis()
  );
  */
  // ADDED
  kernels::aderdg::generic::spaceTimePredictor<testFlux>( lQi, lFi, luh, dx[0], timeStepSize,
                                                          5, // getNumberOfVariables(),
                                                          4  // getNodesPerCoordinateAxis()
  );
  kernels::aderdg::generic::predictor( lQhi, lFhi, lQi, lFi, timeStepSize,
                                                          5, // getNumberOfVariables(),
                                                          4  // getNodesPerCoordinateAxis()
  );
  kernels::aderdg::generic::extrapolatedPredictor( lQhbnd, lFhbnd, lQhi, lFhi, timeStepSize,
                                                          5, // getNumberOfVariables(),
                                                          4  // getNodesPerCoordinateAxis()
  );


  //validateNumericalEquals(<return vector>, <referencesolution>);
  validateNumericalEqualsWithEps(lQhi[0], 1, eps);
  validateNumericalEqualsWithEps(lQhi[1], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lQhi[2], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lQhi[3], 0, eps);
  validateNumericalEqualsWithEps(lQhi[4], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[5], 1, eps);
  validateNumericalEqualsWithEps(lQhi[6], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lQhi[7], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lQhi[8], 0, eps);
  validateNumericalEqualsWithEps(lQhi[9], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[10], 1, eps);
  validateNumericalEqualsWithEps(lQhi[11], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lQhi[12], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lQhi[13], 0, eps);
  validateNumericalEqualsWithEps(lQhi[14], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[15], 1, eps);
  validateNumericalEqualsWithEps(lQhi[16], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lQhi[17], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lQhi[18], 0, eps);
  validateNumericalEqualsWithEps(lQhi[19], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[20], 1, eps);
  validateNumericalEqualsWithEps(lQhi[21], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lQhi[22], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lQhi[23], 0, eps);
  validateNumericalEqualsWithEps(lQhi[24], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[25], 1, eps);
  validateNumericalEqualsWithEps(lQhi[26], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lQhi[27], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lQhi[28], 0, eps);
  validateNumericalEqualsWithEps(lQhi[29], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[30], 1, eps);
  validateNumericalEqualsWithEps(lQhi[31], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lQhi[32], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lQhi[33], 0, eps);
  validateNumericalEqualsWithEps(lQhi[34], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[35], 1, eps);
  validateNumericalEqualsWithEps(lQhi[36], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lQhi[37], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lQhi[38], 0, eps);
  validateNumericalEqualsWithEps(lQhi[39], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[40], 1, eps);
  validateNumericalEqualsWithEps(lQhi[41], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lQhi[42], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lQhi[43], 0, eps);
  validateNumericalEqualsWithEps(lQhi[44], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[45], 1, eps);
  validateNumericalEqualsWithEps(lQhi[46], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lQhi[47], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lQhi[48], 0, eps);
  validateNumericalEqualsWithEps(lQhi[49], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[50], 1, eps);
  validateNumericalEqualsWithEps(lQhi[51], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lQhi[52], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lQhi[53], 0, eps);
  validateNumericalEqualsWithEps(lQhi[54], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[55], 1, eps);
  validateNumericalEqualsWithEps(lQhi[56], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lQhi[57], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lQhi[58], 0, eps);
  validateNumericalEqualsWithEps(lQhi[59], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[60], 1, eps);
  validateNumericalEqualsWithEps(lQhi[61], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lQhi[62], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lQhi[63], 0, eps);
  validateNumericalEqualsWithEps(lQhi[64], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[65], 1, eps);
  validateNumericalEqualsWithEps(lQhi[66], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lQhi[67], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lQhi[68], 0, eps);
  validateNumericalEqualsWithEps(lQhi[69], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[70], 1, eps);
  validateNumericalEqualsWithEps(lQhi[71], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lQhi[72], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lQhi[73], 0, eps);
  validateNumericalEqualsWithEps(lQhi[74], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[75], 1, eps);
  validateNumericalEqualsWithEps(lQhi[76], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lQhi[77], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lQhi[78], 0, eps);
  validateNumericalEqualsWithEps(lQhi[79], 2.5, eps);

  // lFhi_x
  validateNumericalEqualsWithEps(lFhi[0], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lFhi[1], 1, eps);
  validateNumericalEqualsWithEps(lFhi[2], 3.90655760696e-35, eps);
  validateNumericalEqualsWithEps(lFhi[3], 0, eps);
  validateNumericalEqualsWithEps(lFhi[4], -2.02554022444e-17, eps);
  validateNumericalEqualsWithEps(lFhi[5], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lFhi[6], 1, eps);
  validateNumericalEqualsWithEps(lFhi[7], -6.28505465062e-35, eps);
  validateNumericalEqualsWithEps(lFhi[8], 0, eps);
  validateNumericalEqualsWithEps(lFhi[9], 1.16479530929e-17, eps);
  validateNumericalEqualsWithEps(lFhi[10], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lFhi[11], 1, eps);
  validateNumericalEqualsWithEps(lFhi[12], 6.02226203114e-35, eps);
  validateNumericalEqualsWithEps(lFhi[13], 0, eps);
  validateNumericalEqualsWithEps(lFhi[14], -1.07761806819e-17, eps);
  validateNumericalEqualsWithEps(lFhi[15], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lFhi[16], 1, eps);
  validateNumericalEqualsWithEps(lFhi[17], -5.27718213855e-35, eps);
  validateNumericalEqualsWithEps(lFhi[18], 0, eps);
  validateNumericalEqualsWithEps(lFhi[19], 2.68302860327e-17, eps);
  validateNumericalEqualsWithEps(lFhi[20], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lFhi[21], 1, eps);
  validateNumericalEqualsWithEps(lFhi[22], -6.28505465062e-35, eps);
  validateNumericalEqualsWithEps(lFhi[23], 0, eps);
  validateNumericalEqualsWithEps(lFhi[24], -5.18464617038e-17, eps);
  validateNumericalEqualsWithEps(lFhi[25], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lFhi[26], 1, eps);
  validateNumericalEqualsWithEps(lFhi[27], 1.09880295341e-34, eps);
  validateNumericalEqualsWithEps(lFhi[28], 0, eps);
  validateNumericalEqualsWithEps(lFhi[29], 2.67093302613e-17, eps);
  validateNumericalEqualsWithEps(lFhi[30], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lFhi[31], 1, eps);
  validateNumericalEqualsWithEps(lFhi[32], -7.10324058069e-35, eps);
  validateNumericalEqualsWithEps(lFhi[33], 0, eps);
  validateNumericalEqualsWithEps(lFhi[34], -1.4629358797e-17, eps);
  validateNumericalEqualsWithEps(lFhi[35], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lFhi[36], 1, eps);
  validateNumericalEqualsWithEps(lFhi[37], 5.00229224796e-35, eps);
  validateNumericalEqualsWithEps(lFhi[38], 0, eps);
  validateNumericalEqualsWithEps(lFhi[39], 4.31047227278e-17, eps);
  validateNumericalEqualsWithEps(lFhi[40], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lFhi[41], 1, eps);
  validateNumericalEqualsWithEps(lFhi[42], 6.02226203114e-35, eps);
  validateNumericalEqualsWithEps(lFhi[43], 0, eps);
  validateNumericalEqualsWithEps(lFhi[44], -5.18464617038e-17, eps);
  validateNumericalEqualsWithEps(lFhi[45], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lFhi[46], 1, eps);
  validateNumericalEqualsWithEps(lFhi[47], -7.10324058069e-35, eps);
  validateNumericalEqualsWithEps(lFhi[48], 0, eps);
  validateNumericalEqualsWithEps(lFhi[49], 2.67093302613e-17, eps);
  validateNumericalEqualsWithEps(lFhi[50], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lFhi[51], 1, eps);
  validateNumericalEqualsWithEps(lFhi[52], 4.8794232613e-35, eps);
  validateNumericalEqualsWithEps(lFhi[53], 0, eps);
  validateNumericalEqualsWithEps(lFhi[54], -1.4629358797e-17, eps);
  validateNumericalEqualsWithEps(lFhi[55], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lFhi[56], 1, eps);
  validateNumericalEqualsWithEps(lFhi[57], -4.73179353928e-35, eps);
  validateNumericalEqualsWithEps(lFhi[58], 0, eps);
  validateNumericalEqualsWithEps(lFhi[59], 4.31047227278e-17, eps);
  validateNumericalEqualsWithEps(lFhi[60], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lFhi[61], 1, eps);
  validateNumericalEqualsWithEps(lFhi[62], -5.27718213855e-35, eps);
  validateNumericalEqualsWithEps(lFhi[63], 0, eps);
  validateNumericalEqualsWithEps(lFhi[64], -2.02554022444e-17, eps);
  validateNumericalEqualsWithEps(lFhi[65], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lFhi[66], 1, eps);
  validateNumericalEqualsWithEps(lFhi[67], 5.00229224796e-35, eps);
  validateNumericalEqualsWithEps(lFhi[68], 0, eps);
  validateNumericalEqualsWithEps(lFhi[69], 1.16479530929e-17, eps);
  validateNumericalEqualsWithEps(lFhi[70], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lFhi[71], 1, eps);
  validateNumericalEqualsWithEps(lFhi[72], -4.73179353928e-35, eps);
  validateNumericalEqualsWithEps(lFhi[73], 0, eps);
  validateNumericalEqualsWithEps(lFhi[74], -1.07761806819e-17, eps);
  validateNumericalEqualsWithEps(lFhi[75], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lFhi[76], 1, eps);
  validateNumericalEqualsWithEps(lFhi[77], 7.37687416339e-35, eps);
  validateNumericalEqualsWithEps(lFhi[78], 0, eps);
  validateNumericalEqualsWithEps(lFhi[79], 2.68302860327e-17, eps);
  // lFhi_y
  validateNumericalEqualsWithEps(lFhi[80], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lFhi[81], 3.90655760696e-35, eps);
  validateNumericalEqualsWithEps(lFhi[82], 1, eps);
  validateNumericalEqualsWithEps(lFhi[83], 0, eps);
  validateNumericalEqualsWithEps(lFhi[84], -2.02554022444e-17, eps);
  validateNumericalEqualsWithEps(lFhi[85], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lFhi[86], -6.28505465062e-35, eps);
  validateNumericalEqualsWithEps(lFhi[87], 1, eps);
  validateNumericalEqualsWithEps(lFhi[88], 0, eps);
  validateNumericalEqualsWithEps(lFhi[89], -5.18464617038e-17, eps);
  validateNumericalEqualsWithEps(lFhi[90], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lFhi[91], 6.02226203114e-35, eps);
  validateNumericalEqualsWithEps(lFhi[92], 1, eps);
  validateNumericalEqualsWithEps(lFhi[93], 0, eps);
  validateNumericalEqualsWithEps(lFhi[94], -5.18464617038e-17, eps);
  validateNumericalEqualsWithEps(lFhi[95], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lFhi[96], -5.27718213855e-35, eps);
  validateNumericalEqualsWithEps(lFhi[97], 1, eps);
  validateNumericalEqualsWithEps(lFhi[98], 0, eps);
  validateNumericalEqualsWithEps(lFhi[99], -2.02554022444e-17, eps);
  validateNumericalEqualsWithEps(lFhi[100], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lFhi[101], -6.28505465062e-35, eps);
  validateNumericalEqualsWithEps(lFhi[102], 1, eps);
  validateNumericalEqualsWithEps(lFhi[103], 0, eps);
  validateNumericalEqualsWithEps(lFhi[104], 1.16479530929e-17, eps);
  validateNumericalEqualsWithEps(lFhi[105], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lFhi[106], 1.09880295341e-34, eps);
  validateNumericalEqualsWithEps(lFhi[107], 1, eps);
  validateNumericalEqualsWithEps(lFhi[108], 0, eps);
  validateNumericalEqualsWithEps(lFhi[109], 2.67093302613e-17, eps);
  validateNumericalEqualsWithEps(lFhi[110], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lFhi[111], -7.10324058069e-35, eps);
  validateNumericalEqualsWithEps(lFhi[112], 1, eps);
  validateNumericalEqualsWithEps(lFhi[113], 0, eps);
  validateNumericalEqualsWithEps(lFhi[114], 2.67093302613e-17, eps);
  validateNumericalEqualsWithEps(lFhi[115], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lFhi[116], 5.00229224796e-35, eps);
  validateNumericalEqualsWithEps(lFhi[117], 1, eps);
  validateNumericalEqualsWithEps(lFhi[118], 0, eps);
  validateNumericalEqualsWithEps(lFhi[119], 1.16479530929e-17, eps);
  validateNumericalEqualsWithEps(lFhi[120], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lFhi[121], 6.02226203114e-35, eps);
  validateNumericalEqualsWithEps(lFhi[122], 1, eps);
  validateNumericalEqualsWithEps(lFhi[123], 0, eps);
  validateNumericalEqualsWithEps(lFhi[124], -1.07761806819e-17, eps);
  validateNumericalEqualsWithEps(lFhi[125], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lFhi[126], -7.10324058069e-35, eps);
  validateNumericalEqualsWithEps(lFhi[127], 1, eps);
  validateNumericalEqualsWithEps(lFhi[128], 0, eps);
  validateNumericalEqualsWithEps(lFhi[129], -1.4629358797e-17, eps);
  validateNumericalEqualsWithEps(lFhi[130], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lFhi[131], 4.8794232613e-35, eps);
  validateNumericalEqualsWithEps(lFhi[132], 1, eps);
  validateNumericalEqualsWithEps(lFhi[133], 0, eps);
  validateNumericalEqualsWithEps(lFhi[134], -1.4629358797e-17, eps);
  validateNumericalEqualsWithEps(lFhi[135], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lFhi[136], -4.73179353928e-35, eps);
  validateNumericalEqualsWithEps(lFhi[137], 1, eps);
  validateNumericalEqualsWithEps(lFhi[138], 0, eps);
  validateNumericalEqualsWithEps(lFhi[139], -1.07761806819e-17, eps);
  validateNumericalEqualsWithEps(lFhi[140], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lFhi[141], -5.27718213855e-35, eps);
  validateNumericalEqualsWithEps(lFhi[142], 1, eps);
  validateNumericalEqualsWithEps(lFhi[143], 0, eps);
  validateNumericalEqualsWithEps(lFhi[144], 2.68302860327e-17, eps);
  validateNumericalEqualsWithEps(lFhi[145], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lFhi[146], 5.00229224796e-35, eps);
  validateNumericalEqualsWithEps(lFhi[147], 1, eps);
  validateNumericalEqualsWithEps(lFhi[148], 0, eps);
  validateNumericalEqualsWithEps(lFhi[149], 4.31047227278e-17, eps);
  validateNumericalEqualsWithEps(lFhi[150], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lFhi[151], -4.73179353928e-35, eps);
  validateNumericalEqualsWithEps(lFhi[152], 1, eps);
  validateNumericalEqualsWithEps(lFhi[153], 0, eps);
  validateNumericalEqualsWithEps(lFhi[154], 4.31047227278e-17, eps);
  validateNumericalEqualsWithEps(lFhi[155], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lFhi[156], 7.37687416339e-35, eps);
  validateNumericalEqualsWithEps(lFhi[157], 1, eps);
  validateNumericalEqualsWithEps(lFhi[158], 0, eps);
  validateNumericalEqualsWithEps(lFhi[159], 2.68302860327e-17, eps);

  // lQhbnd
  validateNumericalEqualsWithEps(lQhbnd[0], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[1], -1.3650848498e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[2], -2.06067776336e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[3], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[4], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[5], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[6], -3.19038266038e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[7], 1.55129951621e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[8], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[9], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[10], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[11], -3.19038266038e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[12], -2.6243758438e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[13], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[14], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[15], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[16], -1.3650848498e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[17], 5.74601263789e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[18], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[19], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[20], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[21], 1.620214355e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[22], -2.06067776336e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[23], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[24], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[25], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[26], 2.69499929145e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[27], 1.55129951621e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[28], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[29], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[30], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[31], 2.69499929145e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[32], -2.6243758438e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[33], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[34], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[35], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[36], 1.620214355e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[37], 5.74601263789e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[38], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[39], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[40], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[41], -2.06067776336e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[42], -1.3650848498e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[43], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[44], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[45], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[46], 1.55129951621e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[47], -3.19038266038e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[48], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[49], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[50], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[51], -2.6243758438e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[52], -3.19038266038e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[53], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[54], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[55], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[56], 5.74601263789e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[57], -1.3650848498e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[58], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[59], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[60], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[61], -2.06067776336e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[62], 1.620214355e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[63], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[64], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[65], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[66], 1.55129951621e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[67], 2.69499929145e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[68], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[69], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[70], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[71], -2.6243758438e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[72], 2.69499929145e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[73], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[74], 2.5, eps);
  validateNumericalEqualsWithEps(lQhbnd[75], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lQhbnd[76], 5.74601263789e-18, eps);
  validateNumericalEqualsWithEps(lQhbnd[77], 1.620214355e-17, eps);
  validateNumericalEqualsWithEps(lQhbnd[78], 0, eps);
  validateNumericalEqualsWithEps(lQhbnd[79], 2.5, eps);

  // lFhbnd
  validateNumericalEqualsWithEps(lFhbnd[0], -1.3650848498e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[1], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[2], 1.40928628571e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[3], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[4], -4.77779697431e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[5], -3.19038266038e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[6], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[7], -2.19527167955e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[8], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[9], -1.11663393113e-16, eps);
  validateNumericalEqualsWithEps(lFhbnd[10], -3.19038266038e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[11], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[12], 1.74686629283e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[13], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[14], -1.11663393113e-16, eps);
  validateNumericalEqualsWithEps(lFhbnd[15], -1.3650848498e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[16], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[17], -1.48638399144e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[18], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[19], -4.77779697431e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[20], 1.620214355e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[21], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[22], -1.59208789796e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[23], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[24], 5.6707502425e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[25], 2.69499929145e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[26], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[27], 1.85364226628e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[28], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[29], 9.43249752007e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[30], 2.69499929145e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[31], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[32], -1.47272479882e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[33], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[34], 9.43249752007e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[35], 1.620214355e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[36], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[37], 1.77187526833e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[38], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[39], 5.6707502425e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[40], -1.3650848498e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[41], 1.40928628571e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[42], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[43], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[44], -4.77779697431e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[45], -3.19038266038e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[46], -2.19527167955e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[47], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[48], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[49], -1.11663393113e-16, eps);
  validateNumericalEqualsWithEps(lFhbnd[50], -3.19038266038e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[51], 1.74686629283e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[52], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[53], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[54], -1.11663393113e-16, eps);
  validateNumericalEqualsWithEps(lFhbnd[55], -1.3650848498e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[56], -1.48638399144e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[57], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[58], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[59], -4.77779697431e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[60], 1.620214355e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[61], -1.59208789796e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[62], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[63], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[64], 5.6707502425e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[65], 2.69499929145e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[66], 1.85364226628e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[67], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[68], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[69], 9.43249752007e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[70], 2.69499929145e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[71], -1.47272479882e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[72], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[73], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[74], 9.43249752007e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[75], 1.620214355e-17, eps);
  validateNumericalEqualsWithEps(lFhbnd[76], 1.77187526833e-34, eps);
  validateNumericalEqualsWithEps(lFhbnd[77], 0.999999999999, eps);
  validateNumericalEqualsWithEps(lFhbnd[78], 0, eps);
  validateNumericalEqualsWithEps(lFhbnd[79], 5.6707502425e-17, eps);

  delete[] luh;
  delete[] lQi;
  delete[] lFi;
  delete[] rhs0;
  delete[] rhs;
  delete[] tmp;
  delete[] lQhi;
  delete[] lFhi;
  delete[] lQhbnd;
  delete[] lFhbnd;
  //  #endif
} // testSpaceTimePredictor2d


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",on)
#endif
