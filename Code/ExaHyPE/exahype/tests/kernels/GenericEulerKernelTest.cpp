#include "exahype/tests/kernels/GenericEulerKernelTest.h"


#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"


#include "kernels/aderdg/generic/Kernels.h"

using std::cout;
using std::endl;

registerTest(exahype::tests::GenericEulerKernelTest)


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",off)
#endif


exahype::tests::GenericEulerKernelTest::GenericEulerKernelTest():
tarch::tests::TestCase( "exahype::tests::GenericEulerKernelTest" ) {
}


exahype::tests::GenericEulerKernelTest::~GenericEulerKernelTest() {
}

void exahype::tests::GenericEulerKernelTest::run() {
#if DIMENSIONS == 2
  testMethod( testPDEFluxes2d );

  testMethod( testVolumeIntegral2d );
  testMethod( testRiemannSolver2d );
  testMethod( testSurfaceIntegral2d );
  testMethod( testSolutionUpdate2d     );

  testMethod( testSpaceTimePredictor2d );

#elif DIMENSIONS == 3
  testMethod( testPDEFluxes3d );

  testMethod( testVolumeIntegral3d );
  testMethod( testSurfaceIntegral3d );
  testMethod( testSolutionUpdate3d );

#endif

}

#if DIMENSIONS == 2
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



void exahype::tests::GenericEulerKernelTest::testPDEFluxes2d() {
  cout << "Test PDE-related functions, DIM=2" << endl;

  // 2D
  double Q[5] = {1., 0.1, 0.2, 0.3, 3.5}; // pressure = 1.39
  double f[5], g[5];

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

} // testPDEFluxes2d

void exahype::tests::GenericEulerKernelTest::testSolutionUpdate2d() {
  cout << "Test solution update, ORDER=3, DIM=2" << endl;

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
  cout << "Test surface integral, ORDER=3, DIM=2" << endl;

  // inputs:
  const double dx[2] = {0.1, 0.1};  // mesh spacing
  double lFhbnd[4*20] = {0.};

  double* FLeft  = &lFhbnd[ 0];
  double* FRight = &lFhbnd[20];
  double* FFront = &lFhbnd[40];
  double* FBack  = &lFhbnd[60];

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

  // lFhbnd = [ FLeft | FRight | FFront | FBack ]
  kernels::aderdg::generic::surfaceIntegral( lduh, lFhbnd, dx[0],
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
  cout << "Test Riemann Solver (Rusanov), ORDER=3, DIM=2" << endl;
  // input:
  double QL[20] = {1.,0.,0.,0.,2.5,1.,0.,0.,0.,2.5,1.,0.,0.,0.,2.5,1.,0.,0.,0.,2.5};
  double QR[20] = {1.,0.,0.,0.,2.5,1.,0.,0.,0.,2.5,1.,0.,0.,0.,2.5,1.,0.,0.,0.,2.5};

  // output:
  double *FL = new double[20];
  double *FR = new double[20];


  kernels::aderdg::generic::riemannSolver<testEigenvalues>( FL, FR, QL, QR,
                                                            0.0, // dt
                                                            0,   // normalNonZero
                                                            5,   //getNumberOfVariables(),
                                                            4   // getNodesPerCoordinateAxis()
  );

  // FL == FR, element by element
  for(int i=0;i<20;i++) {
    validateEquals(FL[i],FR[i]);
  }

  delete[] FL;
  delete[] FR;
} // testRiemannSolver

void exahype::tests::GenericEulerKernelTest::testVolumeIntegral2d() {
  cout << "Test volume integral, ORDER=3, DIM=2" << endl;

  { // first test

    // output:
    double *lduh = new double[80];

    // input:
    double dx[2] = {3.70370370370370349811e-02, 3.70370370370370349811e-02}; // mesh spacing
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
  } // scope limiter first test


  { // second test, analogous to 3d seed

    // input:
    double dx[2] = {0.05, 0.05};          // mesh spacing
    double * lFhi = new double[160]();  // nVar * dim * nDOFx * nDOFy
    // lFhi = [ lFhi_x | lFhi_y ]
    double * lFhi_x = &lFhi[0];
    double * lFhi_y = &lFhi[80];

    // seed direction
    for(int i=0; i<80; i+=5) {
      lFhi_x[i+1] = 1.;
      lFhi_y[i+2] = 1.;
    }

    // output:
    double *lduh = new double[80]; // intentionally left uninitialised

    kernels::aderdg::generic::volumeIntegral(lduh, lFhi, dx[0], 5, 4);

    validateNumericalEqualsWithEps(lduh[ 0],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[ 1],  -5.70727295609806     , eps);
    validateNumericalEqualsWithEps(lduh[ 2],  -5.70727295609806     , eps);
    validateNumericalEqualsWithEps(lduh[ 3],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[ 4],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[ 5],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[ 6],   4.22432826300142     , eps);
    validateNumericalEqualsWithEps(lduh[ 7],  -10.6997802612945     , eps);
    validateNumericalEqualsWithEps(lduh[ 8],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[ 9],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[10],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[11],  -4.22432826300142     , eps);
    validateNumericalEqualsWithEps(lduh[12],  -10.6997802612945     , eps);
    validateNumericalEqualsWithEps(lduh[13],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[14],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[15],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[16],   5.70727295609806     , eps);
    validateNumericalEqualsWithEps(lduh[17],  -5.70727295609806     , eps);
    validateNumericalEqualsWithEps(lduh[18],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[19],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[20],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[21],  -10.6997802612945     , eps);
    validateNumericalEqualsWithEps(lduh[22],   4.22432826300142     , eps);
    validateNumericalEqualsWithEps(lduh[23],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[24],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[25],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[26],   7.91961143498436     , eps);
    validateNumericalEqualsWithEps(lduh[27],   7.91961143498436     , eps);
    validateNumericalEqualsWithEps(lduh[28],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[29],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[30],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[31],  -7.91961143498436     , eps);
    validateNumericalEqualsWithEps(lduh[32],   7.91961143498436     , eps);
    validateNumericalEqualsWithEps(lduh[33],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[34],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[35],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[36],   10.6997802612945     , eps);
    validateNumericalEqualsWithEps(lduh[37],   4.22432826300142     , eps);
    validateNumericalEqualsWithEps(lduh[38],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[39],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[40],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[41],  -10.6997802612945     , eps);
    validateNumericalEqualsWithEps(lduh[42],  -4.22432826300142     , eps);
    validateNumericalEqualsWithEps(lduh[43],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[44],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[45],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[46],   7.91961143498436     , eps);
    validateNumericalEqualsWithEps(lduh[47],  -7.91961143498436     , eps);
    validateNumericalEqualsWithEps(lduh[48],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[49],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[50],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[51],  -7.91961143498436     , eps);
    validateNumericalEqualsWithEps(lduh[52],  -7.91961143498436     , eps);
    validateNumericalEqualsWithEps(lduh[53],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[54],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[55],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[56],   10.6997802612945     , eps);
    validateNumericalEqualsWithEps(lduh[57],  -4.22432826300142     , eps);
    validateNumericalEqualsWithEps(lduh[58],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[59],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[60],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[61],  -5.70727295609806     , eps);
    validateNumericalEqualsWithEps(lduh[62],   5.70727295609806     , eps);
    validateNumericalEqualsWithEps(lduh[63],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[64],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[65],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[66],   4.22432826300142     , eps);
    validateNumericalEqualsWithEps(lduh[67],   10.6997802612945     , eps);
    validateNumericalEqualsWithEps(lduh[68],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[69],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[70],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[71],  -4.22432826300142     , eps);
    validateNumericalEqualsWithEps(lduh[72],   10.6997802612945     , eps);
    validateNumericalEqualsWithEps(lduh[73],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[74],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[75],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[76],   5.70727295609806     , eps);
    validateNumericalEqualsWithEps(lduh[77],   5.70727295609806     , eps);
    validateNumericalEqualsWithEps(lduh[78],  0.000000000000000E+000, eps);
    validateNumericalEqualsWithEps(lduh[79],  0.000000000000000E+000, eps);



    delete[] lFhi;
    delete[] lduh;
  } // scope limiter second test

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

#elif DIMENSIONS == 3
// copied from MyEulerSolver.cpp
void exahype::tests::GenericEulerKernelTest::testFlux(const double* const Q, double* f, double* g, double *h) {
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

  h[0] = Q[3];
  h[1] = irho*Q[3]*Q[1];
  h[2] = irho*Q[3]*Q[2];
  h[3] = irho*Q[3]*Q[3] + p;
  h[4] = irho*Q[3]*(Q[4]+p);
}

void exahype::tests::GenericEulerKernelTest::testPDEFluxes3d() {
  cout << "Test PDE-related functions, DIM=3" << endl;

  double Q[5] = {1., 0.1, 0.2, 0.3, 3.5}; // pressure = 1.39
  double f[5], g[5], h[5];

  testFlux(Q,f,g,h);

  validateNumericalEquals(0.1, f[0]);
  validateNumericalEquals(1.4, f[1]);
  validateNumericalEqualsWithParams1(0.02, f[2], g[1]);
  validateNumericalEquals(0.03, f[3]);
  validateNumericalEquals(0.489, f[4]);

  validateNumericalEquals(0.2, g[0]);
  validateNumericalEquals(1.43, g[2]);
  validateNumericalEquals(0.06, g[3]);
  validateNumericalEquals(0.978,g[4]);

  validateNumericalEquals(0.3, h[0]);
  validateNumericalEquals(0.03,h[1]);
  validateNumericalEquals(0.06,h[2]);
  validateNumericalEquals(1.48,h[3]);
  validateNumericalEquals(1.4670,h[4]);
} // testPDEFluxes3d


void exahype::tests::GenericEulerKernelTest::testVolumeIntegral3d() {
  cout << "Test volume integral, ORDER=3, DIM=3" << endl;

  // output:
  double *lduh = new double[320]; // intentionally left uninitialised

  // input:
  const double dx[3] = {0.05, 0.05, 0.05}; // mesh spacing
  double * lFhi = new double[960](); // nVar * dim * nDOFx * nDOFy * nDOFz
  // lFhi = [ lFhi_x  | lFhi_y | lFhi_z ], 320 entries each
  double * lFhi_x = &lFhi[0];
  double * lFhi_y = &lFhi[320];
  double * lFhi_z = &lFhi[640];

  // seed direction
  for(int i=0;i<320;i+=5) {
    lFhi_x[i+1] = 1.;
    lFhi_y[i+2] = 1.;
    lFhi_z[i+3] = 1.;
  }

  kernels::aderdg::generic::volumeIntegral(lduh, lFhi, dx[0],
                                           5,  // getNumberOfVariables(),
                                           4); //getNodesPerCoordinateAxis()

  validateNumericalEqualsWithEps(lduh[  0],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[  1], -0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[  2], -0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[  3], -0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[  4],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[  5],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[  6],  0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[  7],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[  8],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[  9],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 10],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 11], -0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[ 12],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 13],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 14],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 15],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 16],  0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[ 17], -0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[ 18], -0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[ 19],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 20],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 21],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 22],  0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[ 23],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 24],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 25],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 26],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[ 27],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[ 28],  -3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[ 29],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 30],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 31],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[ 32],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[ 33],  -3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[ 34],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 35],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 36],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 37],  0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[ 38],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 39],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 40],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 41],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 42], -0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[ 43],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 44],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 45],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 46],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[ 47],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[ 48],  -3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[ 49],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 50],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 51],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[ 52],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[ 53],  -3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[ 54],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 55],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 56],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 57], -0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[ 58],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 59],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 60],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 61], -0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[ 62],  0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[ 63], -0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[ 64],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 65],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 66],  0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[ 67],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 68],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 69],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 70],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 71], -0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[ 72],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 73],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 74],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 75],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 76],  0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[ 77],  0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[ 78], -0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[ 79],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 80],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 81],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 82],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 83],  0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[ 84],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 85],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 86],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[ 87],  -3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[ 88],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[ 89],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 90],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 91],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[ 92],  -3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[ 93],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[ 94],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 95],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[ 96],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 97],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[ 98],  0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[ 99],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[100],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[101],  -3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[102],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[103],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[104],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[105],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[106],   2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[107],   2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[108],   2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[109],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[110],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[111],  -2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[112],   2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[113],   2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[114],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[115],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[116],   3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[117],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[118],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[119],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[120],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[121],  -3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[122],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[123],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[124],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[125],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[126],   2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[127],  -2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[128],   2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[129],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[130],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[131],  -2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[132],  -2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[133],   2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[134],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[135],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[136],   3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[137],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[138],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[139],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[140],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[141],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[142],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[143],  0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[144],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[145],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[146],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[147],   3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[148],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[149],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[150],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[151],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[152],   3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[153],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[154],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[155],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[156],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[157],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[158],  0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[159],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[160],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[161],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[162],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[163], -0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[164],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[165],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[166],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[167],  -3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[168],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[169],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[170],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[171],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[172],  -3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[173],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[174],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[175],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[176],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[177],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[178], -0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[179],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[180],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[181],  -3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[182],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[183],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[184],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[185],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[186],   2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[187],   2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[188],  -2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[189],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[190],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[191],  -2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[192],   2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[193],  -2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[194],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[195],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[196],   3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[197],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[198],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[199],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[200],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[201],  -3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[202],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[203],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[204],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[205],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[206],   2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[207],  -2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[208],  -2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[209],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[210],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[211],  -2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[212],  -2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[213],  -2.58236811285953     , eps);
  validateNumericalEqualsWithEps(lduh[214],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[215],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[216],   3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[217],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[218],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[219],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[220],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[221],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[222],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[223], -0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[224],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[225],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[226],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[227],   3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[228],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[229],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[230],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[231],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[232],   3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[233],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[234],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[235],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[236],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[237],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[238], -0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[239],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[240],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[241], -0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[242], -0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[243],  0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[244],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[245],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[246],  0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[247],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[248],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[249],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[250],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[251], -0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[252],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[253],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[254],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[255],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[256],  0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[257], -0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[258],  0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[259],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[260],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[261],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[262],  0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[263],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[264],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[265],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[266],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[267],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[268],   3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[269],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[270],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[271],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[272],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[273],   3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[274],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[275],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[276],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[277],  0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[278],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[279],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[280],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[281],  -1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[282], -0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[283],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[284],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[285],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[286],   1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[287],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[288],   3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[289],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[290],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[291],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[292],  -1.37743760463265     , eps);
  validateNumericalEqualsWithEps(lduh[293],   3.48890492774857     , eps);
  validateNumericalEqualsWithEps(lduh[294],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[295],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[296],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[297], -0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[298],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[299],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[300],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[301], -0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[302],  0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[303],  0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[304],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[305],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[306],  0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[307],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[308],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[309],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[310],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[311], -0.734726526868064     , eps);
  validateNumericalEqualsWithEps(lduh[312],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[313],   1.86098520289870     , eps);
  validateNumericalEqualsWithEps(lduh[314],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[315],  0.000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lduh[316],  0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[317],  0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[318],  0.992651275150334     , eps);
  validateNumericalEqualsWithEps(lduh[319],  0.000000000000000E+000, eps);

  delete[] lduh;
  delete[] lFhi;
} // testVolumeIntegral3d

void exahype::tests::GenericEulerKernelTest::testSurfaceIntegral3d() {
  cout << "Test surface integral, ORDER=3, DIM=3" << endl;

  // inputs:
  const double dx[3] = {0.1, 0.1, 0.1};  // mesh spacing
  double * lFhbnd  = new double[6*80]();   // 480
  double * FLeft   = &lFhbnd[0];
  double * FRight  = &lFhbnd[80];
  double * FFront  = &lFhbnd[160];
  double * FBack   = &lFhbnd[240];
  double * FBottom = &lFhbnd[320];
  double * FTop    = &lFhbnd[400];

  for(int i=0; i<80; i+=5) {
    // in x orientation 1
    FLeft[ i+1] = 1.;
    FRight[i+1] = 1.;
    // in y orientation 1
    FFront[i+2] = 1.;
    FBack[ i+2] = 1.;
    // in z direction 1
    FBottom[i+3] = 1.;
    FTop[   i+3] = 1.;
  }

  // in/out:
  double * lduh = new double[320];
  for(int i=0; i<320; i++) {
    lduh[i] = i/10.;
  }

  // lFhbnd = [ FLeft | FRight | FFront | FBack | FBottom | FTop ]
  kernels::aderdg::generic::surfaceIntegral( lduh, lFhbnd, dx[0],
                                               5,  //getNumberOfVariables(),
                                               4); // getNodesPerCoordinateAxis()

  validateNumericalEqualsWithEps(lduh[  0],  6.845895347955754E-016, eps);
  validateNumericalEqualsWithEps(lduh[  1],  0.596325637575167     , eps);
  validateNumericalEqualsWithEps(lduh[  2],  0.696325637575167     , eps);
  validateNumericalEqualsWithEps(lduh[  3],  0.796325637575167     , eps);
  validateNumericalEqualsWithEps(lduh[  4],  0.400000000000002     , eps);
  validateNumericalEqualsWithEps(lduh[  5],  0.500000000000001     , eps);
  validateNumericalEqualsWithEps(lduh[  6],  0.232636736565968     , eps);
  validateNumericalEqualsWithEps(lduh[  7],   1.63049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[  8],   1.73049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[  9],  0.900000000000002     , eps);
  validateNumericalEqualsWithEps(lduh[ 10],   1.00000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 11],   1.46736326343403     , eps);
  validateNumericalEqualsWithEps(lduh[ 12],   2.13049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[ 13],   2.23049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[ 14],   1.40000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 15],   1.50000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 16],   1.10367436242483     , eps);
  validateNumericalEqualsWithEps(lduh[ 17],   2.19632563757517     , eps);
  validateNumericalEqualsWithEps(lduh[ 18],   2.29632563757517     , eps);
  validateNumericalEqualsWithEps(lduh[ 19],   1.90000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 20],   2.00000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 21],   3.03049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[ 22],   1.83263673656597     , eps);
  validateNumericalEqualsWithEps(lduh[ 23],   3.23049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[ 24],   2.40000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 25],   2.50000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 26],   1.91128119768368     , eps);
  validateNumericalEqualsWithEps(lduh[ 27],   2.01128119768368     , eps);
  validateNumericalEqualsWithEps(lduh[ 28],   4.54445246387428     , eps);
  validateNumericalEqualsWithEps(lduh[ 29],   2.90000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 30],   3.00000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 31],   3.78871880231632     , eps);
  validateNumericalEqualsWithEps(lduh[ 32],   2.51128119768368     , eps);
  validateNumericalEqualsWithEps(lduh[ 33],   5.04445246387428     , eps);
  validateNumericalEqualsWithEps(lduh[ 34],   3.40000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 35],   3.50000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 36],   2.66950739855065     , eps);
  validateNumericalEqualsWithEps(lduh[ 37],   3.33263673656597     , eps);
  validateNumericalEqualsWithEps(lduh[ 38],   4.73049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[ 39],   3.90000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 40],   4.00000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 41],   5.03049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[ 42],   4.56736326343403     , eps);
  validateNumericalEqualsWithEps(lduh[ 43],   5.23049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[ 44],   4.40000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 45],   4.50000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 46],   3.91128119768368     , eps);
  validateNumericalEqualsWithEps(lduh[ 47],   5.38871880231632     , eps);
  validateNumericalEqualsWithEps(lduh[ 48],   6.54445246387428     , eps);
  validateNumericalEqualsWithEps(lduh[ 49],   4.90000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 50],   5.00000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 51],   5.78871880231632     , eps);
  validateNumericalEqualsWithEps(lduh[ 52],   5.88871880231632     , eps);
  validateNumericalEqualsWithEps(lduh[ 53],   7.04445246387428     , eps);
  validateNumericalEqualsWithEps(lduh[ 54],   5.40000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 55],   5.50000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 56],   4.66950739855065     , eps);
  validateNumericalEqualsWithEps(lduh[ 57],   6.06736326343403     , eps);
  validateNumericalEqualsWithEps(lduh[ 58],   6.73049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[ 59],   5.90000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 60],   6.00000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 61],   6.59632563757517     , eps);
  validateNumericalEqualsWithEps(lduh[ 62],   5.70367436242483     , eps);
  validateNumericalEqualsWithEps(lduh[ 63],   6.79632563757517     , eps);
  validateNumericalEqualsWithEps(lduh[ 64],   6.40000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 65],   6.50000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 66],   6.23263673656597     , eps);
  validateNumericalEqualsWithEps(lduh[ 67],   5.76950739855065     , eps);
  validateNumericalEqualsWithEps(lduh[ 68],   7.73049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[ 69],   6.90000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 70],   7.00000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 71],   7.46736326343403     , eps);
  validateNumericalEqualsWithEps(lduh[ 72],   6.26950739855065     , eps);
  validateNumericalEqualsWithEps(lduh[ 73],   8.23049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[ 74],   7.40000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 75],   7.50000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 76],   7.10367436242483     , eps);
  validateNumericalEqualsWithEps(lduh[ 77],   7.20367436242483     , eps);
  validateNumericalEqualsWithEps(lduh[ 78],   8.29632563757517     , eps);
  validateNumericalEqualsWithEps(lduh[ 79],   7.90000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 80],   8.00000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 81],   9.03049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[ 82],   9.13049260144935     , eps);
  validateNumericalEqualsWithEps(lduh[ 83],   7.93263673656597     , eps);
  validateNumericalEqualsWithEps(lduh[ 84],   8.40000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 85],   8.50000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 86],   7.91128119768368     , eps);
  validateNumericalEqualsWithEps(lduh[ 87],   10.4444524638743     , eps);
  validateNumericalEqualsWithEps(lduh[ 88],   8.11128119768368     , eps);
  validateNumericalEqualsWithEps(lduh[ 89],   8.90000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 90],   9.00000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 91],   9.78871880231632     , eps);
  validateNumericalEqualsWithEps(lduh[ 92],   10.9444524638743     , eps);
  validateNumericalEqualsWithEps(lduh[ 93],   8.61128119768368     , eps);
  validateNumericalEqualsWithEps(lduh[ 94],   9.40000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 95],   9.50000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[ 96],   8.66950739855065     , eps);
  validateNumericalEqualsWithEps(lduh[ 97],   10.6304926014493     , eps);
  validateNumericalEqualsWithEps(lduh[ 98],   9.43263673656597     , eps);
  validateNumericalEqualsWithEps(lduh[ 99],   9.90000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[100],   10.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[101],   11.8444524638743     , eps);
  validateNumericalEqualsWithEps(lduh[102],   9.51128119768368     , eps);
  validateNumericalEqualsWithEps(lduh[103],   9.61128119768368     , eps);
  validateNumericalEqualsWithEps(lduh[104],   10.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[105],   10.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[106],   9.30881594357023     , eps);
  validateNumericalEqualsWithEps(lduh[107],   9.40881594357023     , eps);
  validateNumericalEqualsWithEps(lduh[108],   9.50881594357023     , eps);
  validateNumericalEqualsWithEps(lduh[109],   10.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[110],   11.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[111],   12.3911840564298     , eps);
  validateNumericalEqualsWithEps(lduh[112],   9.90881594357023     , eps);
  validateNumericalEqualsWithEps(lduh[113],   10.0088159435702     , eps);
  validateNumericalEqualsWithEps(lduh[114],   11.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[115],   11.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[116],   9.85554753612572     , eps);
  validateNumericalEqualsWithEps(lduh[117],   11.0112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[118],   11.1112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[119],   11.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[120],   12.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[121],   13.8444524638743     , eps);
  validateNumericalEqualsWithEps(lduh[122],   12.8887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[123],   11.6112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[124],   12.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[125],   12.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[126],   11.3088159435702     , eps);
  validateNumericalEqualsWithEps(lduh[127],   13.9911840564298     , eps);
  validateNumericalEqualsWithEps(lduh[128],   11.5088159435702     , eps);
  validateNumericalEqualsWithEps(lduh[129],   12.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[130],   13.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[131],   14.3911840564298     , eps);
  validateNumericalEqualsWithEps(lduh[132],   14.4911840564298     , eps);
  validateNumericalEqualsWithEps(lduh[133],   12.0088159435702     , eps);
  validateNumericalEqualsWithEps(lduh[134],   13.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[135],   13.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[136],   11.8555475361257     , eps);
  validateNumericalEqualsWithEps(lduh[137],   14.3887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[138],   13.1112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[139],   13.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[140],   14.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[141],   15.0304926014493     , eps);
  validateNumericalEqualsWithEps(lduh[142],   13.2695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[143],   13.9326367365660     , eps);
  validateNumericalEqualsWithEps(lduh[144],   14.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[145],   14.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[146],   13.9112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[147],   12.9555475361257     , eps);
  validateNumericalEqualsWithEps(lduh[148],   14.1112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[149],   14.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[150],   15.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[151],   15.7887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[152],   13.4555475361257     , eps);
  validateNumericalEqualsWithEps(lduh[153],   14.6112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[154],   15.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[155],   15.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[156],   14.6695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[157],   14.7695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[158],   15.4326367365660     , eps);
  validateNumericalEqualsWithEps(lduh[159],   15.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[160],   16.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[161],   17.0304926014493     , eps);
  validateNumericalEqualsWithEps(lduh[162],   17.1304926014493     , eps);
  validateNumericalEqualsWithEps(lduh[163],   16.6673632634340     , eps);
  validateNumericalEqualsWithEps(lduh[164],   16.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[165],   16.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[166],   15.9112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[167],   18.4444524638743     , eps);
  validateNumericalEqualsWithEps(lduh[168],   17.4887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[169],   16.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[170],   17.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[171],   17.7887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[172],   18.9444524638743     , eps);
  validateNumericalEqualsWithEps(lduh[173],   17.9887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[174],   17.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[175],   17.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[176],   16.6695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[177],   18.6304926014493     , eps);
  validateNumericalEqualsWithEps(lduh[178],   18.1673632634340     , eps);
  validateNumericalEqualsWithEps(lduh[179],   17.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[180],   18.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[181],   19.8444524638743     , eps);
  validateNumericalEqualsWithEps(lduh[182],   17.5112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[183],   18.9887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[184],   18.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[185],   18.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[186],   17.3088159435702     , eps);
  validateNumericalEqualsWithEps(lduh[187],   17.4088159435702     , eps);
  validateNumericalEqualsWithEps(lduh[188],   20.0911840564298     , eps);
  validateNumericalEqualsWithEps(lduh[189],   18.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[190],   19.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[191],   20.3911840564298     , eps);
  validateNumericalEqualsWithEps(lduh[192],   17.9088159435702     , eps);
  validateNumericalEqualsWithEps(lduh[193],   20.5911840564298     , eps);
  validateNumericalEqualsWithEps(lduh[194],   19.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[195],   19.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[196],   17.8555475361257     , eps);
  validateNumericalEqualsWithEps(lduh[197],   19.0112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[198],   20.4887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[199],   19.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[200],   20.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[201],   21.8444524638743     , eps);
  validateNumericalEqualsWithEps(lduh[202],   20.8887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[203],   20.9887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[204],   20.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[205],   20.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[206],   19.3088159435702     , eps);
  validateNumericalEqualsWithEps(lduh[207],   21.9911840564298     , eps);
  validateNumericalEqualsWithEps(lduh[208],   22.0911840564298     , eps);
  validateNumericalEqualsWithEps(lduh[209],   20.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[210],   21.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[211],   22.3911840564298     , eps);
  validateNumericalEqualsWithEps(lduh[212],   22.4911840564298     , eps);
  validateNumericalEqualsWithEps(lduh[213],   22.5911840564298     , eps);
  validateNumericalEqualsWithEps(lduh[214],   21.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[215],   21.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[216],   19.8555475361257     , eps);
  validateNumericalEqualsWithEps(lduh[217],   22.3887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[218],   22.4887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[219],   21.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[220],   22.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[221],   23.0304926014493     , eps);
  validateNumericalEqualsWithEps(lduh[222],   21.2695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[223],   22.6673632634340     , eps);
  validateNumericalEqualsWithEps(lduh[224],   22.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[225],   22.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[226],   21.9112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[227],   20.9555475361257     , eps);
  validateNumericalEqualsWithEps(lduh[228],   23.4887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[229],   22.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[230],   23.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[231],   23.7887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[232],   21.4555475361257     , eps);
  validateNumericalEqualsWithEps(lduh[233],   23.9887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[234],   23.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[235],   23.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[236],   22.6695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[237],   22.7695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[238],   24.1673632634340     , eps);
  validateNumericalEqualsWithEps(lduh[239],   23.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[240],   24.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[241],   24.5963256375752     , eps);
  validateNumericalEqualsWithEps(lduh[242],   24.6963256375752     , eps);
  validateNumericalEqualsWithEps(lduh[243],   23.8036743624248     , eps);
  validateNumericalEqualsWithEps(lduh[244],   24.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[245],   24.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[246],   24.2326367365660     , eps);
  validateNumericalEqualsWithEps(lduh[247],   25.6304926014493     , eps);
  validateNumericalEqualsWithEps(lduh[248],   23.8695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[249],   24.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[250],   25.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[251],   25.4673632634340     , eps);
  validateNumericalEqualsWithEps(lduh[252],   26.1304926014493     , eps);
  validateNumericalEqualsWithEps(lduh[253],   24.3695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[254],   25.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[255],   25.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[256],   25.1036743624248     , eps);
  validateNumericalEqualsWithEps(lduh[257],   26.1963256375752     , eps);
  validateNumericalEqualsWithEps(lduh[258],   25.3036743624248     , eps);
  validateNumericalEqualsWithEps(lduh[259],   25.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[260],   26.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[261],   27.0304926014493     , eps);
  validateNumericalEqualsWithEps(lduh[262],   25.8326367365660     , eps);
  validateNumericalEqualsWithEps(lduh[263],   25.3695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[264],   26.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[265],   26.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[266],   25.9112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[267],   26.0112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[268],   25.0555475361257     , eps);
  validateNumericalEqualsWithEps(lduh[269],   26.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[270],   27.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[271],   27.7887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[272],   26.5112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[273],   25.5555475361257     , eps);
  validateNumericalEqualsWithEps(lduh[274],   27.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[275],   27.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[276],   26.6695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[277],   27.3326367365660     , eps);
  validateNumericalEqualsWithEps(lduh[278],   26.8695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[279],   27.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[280],   28.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[281],   29.0304926014493     , eps);
  validateNumericalEqualsWithEps(lduh[282],   28.5673632634340     , eps);
  validateNumericalEqualsWithEps(lduh[283],   27.3695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[284],   28.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[285],   28.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[286],   27.9112811976837     , eps);
  validateNumericalEqualsWithEps(lduh[287],   29.3887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[288],   27.0555475361257     , eps);
  validateNumericalEqualsWithEps(lduh[289],   28.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[290],   29.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[291],   29.7887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[292],   29.8887188023163     , eps);
  validateNumericalEqualsWithEps(lduh[293],   27.5555475361257     , eps);
  validateNumericalEqualsWithEps(lduh[294],   29.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[295],   29.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[296],   28.6695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[297],   30.0673632634340     , eps);
  validateNumericalEqualsWithEps(lduh[298],   28.8695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[299],   29.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[300],   30.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[301],   30.5963256375752     , eps);
  validateNumericalEqualsWithEps(lduh[302],   29.7036743624248     , eps);
  validateNumericalEqualsWithEps(lduh[303],   29.8036743624248     , eps);
  validateNumericalEqualsWithEps(lduh[304],   30.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[305],   30.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[306],   30.2326367365660     , eps);
  validateNumericalEqualsWithEps(lduh[307],   29.7695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[308],   29.8695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[309],   30.9000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[310],   31.0000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[311],   31.4673632634340     , eps);
  validateNumericalEqualsWithEps(lduh[312],   30.2695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[313],   30.3695073985507     , eps);
  validateNumericalEqualsWithEps(lduh[314],   31.4000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[315],   31.5000000000000     , eps);
  validateNumericalEqualsWithEps(lduh[316],   31.1036743624248     , eps);
  validateNumericalEqualsWithEps(lduh[317],   31.2036743624248     , eps);
  validateNumericalEqualsWithEps(lduh[318],   31.3036743624248     , eps);
  validateNumericalEqualsWithEps(lduh[319],   31.9000000000000     , eps);

  delete[] lFhbnd;
  delete[] lduh;
} // testSurfaceIntegral3d


void exahype::tests::GenericEulerKernelTest::testSolutionUpdate3d() {
  cout << "Test solution update, ORDER=3, DIM=3" << endl;

  // in/out:
  double * luh  = new double[320]();
  for(int i=0;i<320;i+=5) {
    luh[i]   = 1.0;
    luh[i+4] = 2.5;
  }

  // inputs:
  const double dt = 1.40831757919882352703e-03;
  double * lduh = new double[320];
  for(int i=0;i<320;i++) {
    lduh[i] = i;
  }

  kernels::aderdg::generic::solutionUpdate( luh, lduh, dt,
                                            5,   //getNumberOfVariables(),
                                            4   // getNodesPerCoordinateAxis()
  );

  validateNumericalEqualsWithEps(luh[0],  1.00000000000000, eps);
  validateNumericalEqualsWithEps(luh[1],  0.267667977113063, eps);
  validateNumericalEqualsWithEps(luh[2],  0.535335954226126, eps);
  validateNumericalEqualsWithEps(luh[3],  0.803003931339190, eps);
  validateNumericalEqualsWithEps(luh[4],  3.57067190845225, eps);
  validateNumericalEqualsWithEps(luh[5],  1.71387176637496, eps);
  validateNumericalEqualsWithEps(luh[6],  0.856646119649957, eps);
  validateNumericalEqualsWithEps(luh[7],  0.999420472924950, eps);
  validateNumericalEqualsWithEps(luh[8],  1.14219482619994, eps);
  validateNumericalEqualsWithEps(luh[9],  3.78496917947494, eps);
  validateNumericalEqualsWithEps(luh[10],   2.42774353274993, eps);
  validateNumericalEqualsWithEps(luh[11],   1.57051788602492, eps);
  validateNumericalEqualsWithEps(luh[12],   1.71329223929991, eps);
  validateNumericalEqualsWithEps(luh[13],   1.85606659257491, eps);
  validateNumericalEqualsWithEps(luh[14],   4.49884094584990, eps);
  validateNumericalEqualsWithEps(luh[15],   5.01501965669595, eps);
  validateNumericalEqualsWithEps(luh[16],   4.28268763380901, eps);
  validateNumericalEqualsWithEps(luh[17],   4.55035561092207, eps);
  validateNumericalEqualsWithEps(luh[18],   4.81802358803514, eps);
  validateNumericalEqualsWithEps(luh[19],   7.58569156514820, eps);
  validateNumericalEqualsWithEps(luh[20],   3.85548706549986, eps);
  validateNumericalEqualsWithEps(luh[21],   2.99826141877485, eps);
  validateNumericalEqualsWithEps(luh[22],   3.14103577204984, eps);
  validateNumericalEqualsWithEps(luh[23],   3.28381012532484, eps);
  validateNumericalEqualsWithEps(luh[24],   5.92658447859983, eps);
  validateNumericalEqualsWithEps(luh[25],   2.90389939179034, eps);
  validateNumericalEqualsWithEps(luh[26],   1.98005536746195, eps);
  validateNumericalEqualsWithEps(luh[27],   2.05621134313357, eps);
  validateNumericalEqualsWithEps(luh[28],   2.13236731880518, eps);
  validateNumericalEqualsWithEps(luh[29],   4.70852329447679, eps);
  validateNumericalEqualsWithEps(luh[30],   3.28467927014841, eps);
  validateNumericalEqualsWithEps(luh[31],   2.36083524582002, eps);
  validateNumericalEqualsWithEps(luh[32],   2.43699122149164, eps);
  validateNumericalEqualsWithEps(luh[33],   2.51314719716325, eps);
  validateNumericalEqualsWithEps(luh[34],   5.08930317283486, eps);
  validateNumericalEqualsWithEps(luh[35],   5.99710236462475, eps);
  validateNumericalEqualsWithEps(luh[36],   5.13987671789974, eps);
  validateNumericalEqualsWithEps(luh[37],   5.28265107117474, eps);
  validateNumericalEqualsWithEps(luh[38],   5.42542542444973, eps);
  validateNumericalEqualsWithEps(luh[39],   8.06819977772472, eps);
  validateNumericalEqualsWithEps(luh[40],   6.71097413099972, eps);
  validateNumericalEqualsWithEps(luh[41],   5.85374848427471, eps);
  validateNumericalEqualsWithEps(luh[42],   5.99652283754970, eps);
  validateNumericalEqualsWithEps(luh[43],   6.13929719082469, eps);
  validateNumericalEqualsWithEps(luh[44],   8.78207154409969, eps);
  validateNumericalEqualsWithEps(luh[45],   4.42701890522261, eps);
  validateNumericalEqualsWithEps(luh[46],   3.50317488089423, eps);
  validateNumericalEqualsWithEps(luh[47],   3.57933085656584, eps);
  validateNumericalEqualsWithEps(luh[48],   3.65548683223745, eps);
  validateNumericalEqualsWithEps(luh[49],   6.23164280790907, eps);
  validateNumericalEqualsWithEps(luh[50],   4.80779878358068, eps);
  validateNumericalEqualsWithEps(luh[51],   3.88395475925229, eps);
  validateNumericalEqualsWithEps(luh[52],   3.96011073492391, eps);
  validateNumericalEqualsWithEps(luh[53],   4.03626671059552, eps);
  validateNumericalEqualsWithEps(luh[54],   6.61242268626713, eps);
  validateNumericalEqualsWithEps(luh[55],   8.85258943012461, eps);
  validateNumericalEqualsWithEps(luh[56],   7.99536378339960, eps);
  validateNumericalEqualsWithEps(luh[57],   8.13813813667459, eps);
  validateNumericalEqualsWithEps(luh[58],   8.28091248994959, eps);
  validateNumericalEqualsWithEps(luh[59],   10.9236868432246, eps);
  validateNumericalEqualsWithEps(luh[60],   17.0600786267838, eps);
  validateNumericalEqualsWithEps(luh[61],   16.3277466038969, eps);
  validateNumericalEqualsWithEps(luh[62],   16.5954145810099, eps);
  validateNumericalEqualsWithEps(luh[63],   16.8630825581230, eps);
  validateNumericalEqualsWithEps(luh[64],   19.6307505352360, eps);
  validateNumericalEqualsWithEps(luh[65],   10.2803329628745, eps);
  validateNumericalEqualsWithEps(luh[66],   9.42310731614953, eps);
  validateNumericalEqualsWithEps(luh[67],   9.56588166942452, eps);
  validateNumericalEqualsWithEps(luh[68],   9.70865602269952, eps);
  validateNumericalEqualsWithEps(luh[69],   12.3514303759745, eps);
  validateNumericalEqualsWithEps(luh[70],   10.9942047292495, eps);
  validateNumericalEqualsWithEps(luh[71],   10.1369790825245, eps);
  validateNumericalEqualsWithEps(luh[72],   10.2797534357995, eps);
  validateNumericalEqualsWithEps(luh[73],   10.4225277890745, eps);
  validateNumericalEqualsWithEps(luh[74],   13.0653021423495, eps);
  validateNumericalEqualsWithEps(luh[75],   21.0750982834797, eps);
  validateNumericalEqualsWithEps(luh[76],   20.3427662605928, eps);
  validateNumericalEqualsWithEps(luh[77],   20.6104342377059, eps);
  validateNumericalEqualsWithEps(luh[78],   20.8781022148189, eps);
  validateNumericalEqualsWithEps(luh[79],   23.6457701919320, eps);
  validateNumericalEqualsWithEps(luh[80],   12.4219482619994, eps);
  validateNumericalEqualsWithEps(luh[81],   11.5647226152744, eps);
  validateNumericalEqualsWithEps(luh[82],   11.7074969685494, eps);
  validateNumericalEqualsWithEps(luh[83],   11.8502713218244, eps);
  validateNumericalEqualsWithEps(luh[84],   14.4930456750994, eps);
  validateNumericalEqualsWithEps(luh[85],   7.47325793208716, eps);
  validateNumericalEqualsWithEps(luh[86],   6.54941390775877, eps);
  validateNumericalEqualsWithEps(luh[87],   6.62556988343038, eps);
  validateNumericalEqualsWithEps(luh[88],   6.70172585910200, eps);
  validateNumericalEqualsWithEps(luh[89],   9.27788183477361, eps);
  validateNumericalEqualsWithEps(luh[90],   7.85403781044523, eps);
  validateNumericalEqualsWithEps(luh[91],   6.93019378611684, eps);
  validateNumericalEqualsWithEps(luh[92],   7.00634976178845, eps);
  validateNumericalEqualsWithEps(luh[93],   7.08250573746007, eps);
  validateNumericalEqualsWithEps(luh[94],   9.65866171313168, eps);
  validateNumericalEqualsWithEps(luh[95],   14.5635635611243, eps);
  validateNumericalEqualsWithEps(luh[96],   13.7063379143993, eps);
  validateNumericalEqualsWithEps(luh[97],   13.8491122676743, eps);
  validateNumericalEqualsWithEps(luh[98],   13.9918866209493, eps);
  validateNumericalEqualsWithEps(luh[99],   16.6346609742243, eps);
  validateNumericalEqualsWithEps(luh[100],   8.61559756716136, eps);
  validateNumericalEqualsWithEps(luh[101],   7.69175354283298, eps);
  validateNumericalEqualsWithEps(luh[102],   7.76790951850459, eps);
  validateNumericalEqualsWithEps(luh[103],   7.84406549417620, eps);
  validateNumericalEqualsWithEps(luh[104],   10.4202214698478, eps);
  validateNumericalEqualsWithEps(luh[105],   5.26527532594805, eps);
  validateNumericalEqualsWithEps(luh[106],   4.30589699571898, eps);
  validateNumericalEqualsWithEps(luh[107],   4.34651866548992, eps);
  validateNumericalEqualsWithEps(luh[108],   4.38714033526085, eps);
  validateNumericalEqualsWithEps(luh[109],   6.92776200503179, eps);
  validateNumericalEqualsWithEps(luh[110],   5.46838367480272, eps);
  validateNumericalEqualsWithEps(luh[111],   4.50900534457365, eps);
  validateNumericalEqualsWithEps(luh[112],   4.54962701434459, eps);
  validateNumericalEqualsWithEps(luh[113],   4.59024868411552, eps);
  validateNumericalEqualsWithEps(luh[114],   7.13087035388646, eps);
  validateNumericalEqualsWithEps(luh[115],   9.75793720223557, eps);
  validateNumericalEqualsWithEps(luh[116],   8.83409317790718, eps);
  validateNumericalEqualsWithEps(luh[117],   8.91024915357879, eps);
  validateNumericalEqualsWithEps(luh[118],   8.98640512925041, eps);
  validateNumericalEqualsWithEps(luh[119],   11.5625611049220, eps);
  validateNumericalEqualsWithEps(luh[120],   10.1387170805936, eps);
  validateNumericalEqualsWithEps(luh[121],   9.21487305626525, eps);
  validateNumericalEqualsWithEps(luh[122],   9.29102903193686, eps);
  validateNumericalEqualsWithEps(luh[123],   9.36718500760847, eps);
  validateNumericalEqualsWithEps(luh[124],   11.9433409832801, eps);
  validateNumericalEqualsWithEps(luh[125],   6.07770872136673, eps);
  validateNumericalEqualsWithEps(luh[126],   5.11833039113766, eps);
  validateNumericalEqualsWithEps(luh[127],   5.15895206090859, eps);
  validateNumericalEqualsWithEps(luh[128],   5.19957373067953, eps);
  validateNumericalEqualsWithEps(luh[129],   7.74019540045046, eps);
  validateNumericalEqualsWithEps(luh[130],   6.28081707022140, eps);
  validateNumericalEqualsWithEps(luh[131],   5.32143873999233, eps);
  validateNumericalEqualsWithEps(luh[132],   5.36206040976326, eps);
  validateNumericalEqualsWithEps(luh[133],   5.40268207953420, eps);
  validateNumericalEqualsWithEps(luh[134],   7.94330374930513, eps);
  validateNumericalEqualsWithEps(luh[135],   11.2810567156678, eps);
  validateNumericalEqualsWithEps(luh[136],   10.3572126913395, eps);
  validateNumericalEqualsWithEps(luh[137],   10.4333686670111, eps);
  validateNumericalEqualsWithEps(luh[138],   10.5095246426827, eps);
  validateNumericalEqualsWithEps(luh[139],   13.0856806183543, eps);
  validateNumericalEqualsWithEps(luh[140],   20.9884094584990, eps);
  validateNumericalEqualsWithEps(luh[141],   20.1311838117740, eps);
  validateNumericalEqualsWithEps(luh[142],   20.2739581650490, eps);
  validateNumericalEqualsWithEps(luh[143],   20.4167325183240, eps);
  validateNumericalEqualsWithEps(luh[144],   23.0595068715990, eps);
  validateNumericalEqualsWithEps(luh[145],   12.0426164723840, eps);
  validateNumericalEqualsWithEps(luh[146],   11.1187724480556, eps);
  validateNumericalEqualsWithEps(luh[147],   11.1949284237272, eps);
  validateNumericalEqualsWithEps(luh[148],   11.2710843993988, eps);
  validateNumericalEqualsWithEps(luh[149],   13.8472403750704, eps);
  validateNumericalEqualsWithEps(luh[150],   12.4233963507420, eps);
  validateNumericalEqualsWithEps(luh[151],   11.4995523264137, eps);
  validateNumericalEqualsWithEps(luh[152],   11.5757083020853, eps);
  validateNumericalEqualsWithEps(luh[153],   11.6518642777569, eps);
  validateNumericalEqualsWithEps(luh[154],   14.2280202534285, eps);
  validateNumericalEqualsWithEps(luh[155],   23.1300247576239, eps);
  validateNumericalEqualsWithEps(luh[156],   22.2727991108989, eps);
  validateNumericalEqualsWithEps(luh[157],   22.4155734641739, eps);
  validateNumericalEqualsWithEps(luh[158],   22.5583478174489, eps);
  validateNumericalEqualsWithEps(luh[159],   25.2011221707239, eps);
  validateNumericalEqualsWithEps(luh[160],   23.8438965239989, eps);
  validateNumericalEqualsWithEps(luh[161],   22.9866708772739, eps);
  validateNumericalEqualsWithEps(luh[162],   23.1294452305488, eps);
  validateNumericalEqualsWithEps(luh[163],   23.2722195838238, eps);
  validateNumericalEqualsWithEps(luh[164],   25.9149939370988, eps);
  validateNumericalEqualsWithEps(luh[165],   13.5657359858162, eps);
  validateNumericalEqualsWithEps(luh[166],   12.6418919614879, eps);
  validateNumericalEqualsWithEps(luh[167],   12.7180479371595, eps);
  validateNumericalEqualsWithEps(luh[168],   12.7942039128311, eps);
  validateNumericalEqualsWithEps(luh[169],   15.3703598885027, eps);
  validateNumericalEqualsWithEps(luh[170],   13.9465158641743, eps);
  validateNumericalEqualsWithEps(luh[171],   13.0226718398459, eps);
  validateNumericalEqualsWithEps(luh[172],   13.0988278155175, eps);
  validateNumericalEqualsWithEps(luh[173],   13.1749837911892, eps);
  validateNumericalEqualsWithEps(luh[174],   15.7511397668608, eps);
  validateNumericalEqualsWithEps(luh[175],   25.9855118231238, eps);
  validateNumericalEqualsWithEps(luh[176],   25.1282861763987, eps);
  validateNumericalEqualsWithEps(luh[177],   25.2710605296737, eps);
  validateNumericalEqualsWithEps(luh[178],   25.4138348829487, eps);
  validateNumericalEqualsWithEps(luh[179],   28.0566092362237, eps);
  validateNumericalEqualsWithEps(luh[180],   14.7080756208905, eps);
  validateNumericalEqualsWithEps(luh[181],   13.7842315965621, eps);
  validateNumericalEqualsWithEps(luh[182],   13.8603875722337, eps);
  validateNumericalEqualsWithEps(luh[183],   13.9365435479053, eps);
  validateNumericalEqualsWithEps(luh[184],   16.5126995235769, eps);
  validateNumericalEqualsWithEps(luh[185],   8.51500890762276, eps);
  validateNumericalEqualsWithEps(luh[186],   7.55563057739369, eps);
  validateNumericalEqualsWithEps(luh[187],   7.59625224716462, eps);
  validateNumericalEqualsWithEps(luh[188],   7.63687391693556, eps);
  validateNumericalEqualsWithEps(luh[189],   10.1774955867065, eps);
  validateNumericalEqualsWithEps(luh[190],   8.71811725647743, eps);
  validateNumericalEqualsWithEps(luh[191],   7.75873892624836, eps);
  validateNumericalEqualsWithEps(luh[192],   7.79936059601929, eps);
  validateNumericalEqualsWithEps(luh[193],   7.83998226579023, eps);
  validateNumericalEqualsWithEps(luh[194],   10.3806039355612, eps);
  validateNumericalEqualsWithEps(luh[195],   15.8504152559647, eps);
  validateNumericalEqualsWithEps(luh[196],   14.9265712316363, eps);
  validateNumericalEqualsWithEps(luh[197],   15.0027272073079, eps);
  validateNumericalEqualsWithEps(luh[198],   15.0788831829795, eps);
  validateNumericalEqualsWithEps(luh[199],   17.6550391586511, eps);
  validateNumericalEqualsWithEps(luh[200],   16.2311951343227, eps);
  validateNumericalEqualsWithEps(luh[201],   15.3073511099943, eps);
  validateNumericalEqualsWithEps(luh[202],   15.3835070856660, eps);
  validateNumericalEqualsWithEps(luh[203],   15.4596630613376, eps);
  validateNumericalEqualsWithEps(luh[204],   18.0358190370092, eps);
  validateNumericalEqualsWithEps(luh[205],   9.32744230304143, eps);
  validateNumericalEqualsWithEps(luh[206],   8.36806397281237, eps);
  validateNumericalEqualsWithEps(luh[207],   8.40868564258330, eps);
  validateNumericalEqualsWithEps(luh[208],   8.44930731235423, eps);
  validateNumericalEqualsWithEps(luh[209],   10.9899289821252, eps);
  validateNumericalEqualsWithEps(luh[210],   9.53055065189610, eps);
  validateNumericalEqualsWithEps(luh[211],   8.57117232166704, eps);
  validateNumericalEqualsWithEps(luh[212],   8.61179399143797, eps);
  validateNumericalEqualsWithEps(luh[213],   8.65241566120890, eps);
  validateNumericalEqualsWithEps(luh[214],   11.1930373309798, eps);
  validateNumericalEqualsWithEps(luh[215],   17.3735347693969, eps);
  validateNumericalEqualsWithEps(luh[216],   16.4496907450685, eps);
  validateNumericalEqualsWithEps(luh[217],   16.5258467207402, eps);
  validateNumericalEqualsWithEps(luh[218],   16.6020026964118, eps);
  validateNumericalEqualsWithEps(luh[219],   19.1781586720834, eps);
  validateNumericalEqualsWithEps(luh[220],   32.4103577204984, eps);
  validateNumericalEqualsWithEps(luh[221],   31.5531320737734, eps);
  validateNumericalEqualsWithEps(luh[222],   31.6959064270484, eps);
  validateNumericalEqualsWithEps(luh[223],   31.8386807803234, eps);
  validateNumericalEqualsWithEps(luh[224],   34.4814551335984, eps);
  validateNumericalEqualsWithEps(luh[225],   18.1350945261131, eps);
  validateNumericalEqualsWithEps(luh[226],   17.2112505017847, eps);
  validateNumericalEqualsWithEps(luh[227],   17.2874064774563, eps);
  validateNumericalEqualsWithEps(luh[228],   17.3635624531279, eps);
  validateNumericalEqualsWithEps(luh[229],   19.9397184287995, eps);
  validateNumericalEqualsWithEps(luh[230],   18.5158744044711, eps);
  validateNumericalEqualsWithEps(luh[231],   17.5920303801427, eps);
  validateNumericalEqualsWithEps(luh[232],   17.6681863558144, eps);
  validateNumericalEqualsWithEps(luh[233],   17.7443423314860, eps);
  validateNumericalEqualsWithEps(luh[234],   20.3204983071576, eps);
  validateNumericalEqualsWithEps(luh[235],   34.5519730196233, eps);
  validateNumericalEqualsWithEps(luh[236],   33.6947473728983, eps);
  validateNumericalEqualsWithEps(luh[237],   33.8375217261733, eps);
  validateNumericalEqualsWithEps(luh[238],   33.9802960794483, eps);
  validateNumericalEqualsWithEps(luh[239],   36.6230704327233, eps);
  validateNumericalEqualsWithEps(luh[240],   65.2403145071352, eps);
  validateNumericalEqualsWithEps(luh[241],   64.5079824842482, eps);
  validateNumericalEqualsWithEps(luh[242],   64.7756504613613, eps);
  validateNumericalEqualsWithEps(luh[243],   65.0433184384744, eps);
  validateNumericalEqualsWithEps(luh[244],   67.8109864155874, eps);
  validateNumericalEqualsWithEps(luh[245],   35.9797165523733, eps);
  validateNumericalEqualsWithEps(luh[246],   35.1224909056482, eps);
  validateNumericalEqualsWithEps(luh[247],   35.2652652589232, eps);
  validateNumericalEqualsWithEps(luh[248],   35.4080396121982, eps);
  validateNumericalEqualsWithEps(luh[249],   38.0508139654732, eps);
  validateNumericalEqualsWithEps(luh[250],   36.6935883187482, eps);
  validateNumericalEqualsWithEps(luh[251],   35.8363626720232, eps);
  validateNumericalEqualsWithEps(luh[252],   35.9791370252982, eps);
  validateNumericalEqualsWithEps(luh[253],   36.1219113785732, eps);
  validateNumericalEqualsWithEps(luh[254],   38.7646857318482, eps);
  validateNumericalEqualsWithEps(luh[255],   69.2553341638311, eps);
  validateNumericalEqualsWithEps(luh[256],   68.5230021409442, eps);
  validateNumericalEqualsWithEps(luh[257],   68.7906701180572, eps);
  validateNumericalEqualsWithEps(luh[258],   69.0583380951703, eps);
  validateNumericalEqualsWithEps(luh[259],   71.8260060722834, eps);
  validateNumericalEqualsWithEps(luh[260],   38.1213318514982, eps);
  validateNumericalEqualsWithEps(luh[261],   37.2641062047731, eps);
  validateNumericalEqualsWithEps(luh[262],   37.4068805580481, eps);
  validateNumericalEqualsWithEps(luh[263],   37.5496549113231, eps);
  validateNumericalEqualsWithEps(luh[264],   40.1924292645981, eps);
  validateNumericalEqualsWithEps(luh[265],   21.1813335529776, eps);
  validateNumericalEqualsWithEps(luh[266],   20.2574895286492, eps);
  validateNumericalEqualsWithEps(luh[267],   20.3336455043208, eps);
  validateNumericalEqualsWithEps(luh[268],   20.4098014799924, eps);
  validateNumericalEqualsWithEps(luh[269],   22.9859574556641, eps);
  validateNumericalEqualsWithEps(luh[270],   21.5621134313357, eps);
  validateNumericalEqualsWithEps(luh[271],   20.6382694070073, eps);
  validateNumericalEqualsWithEps(luh[272],   20.7144253826789, eps);
  validateNumericalEqualsWithEps(luh[273],   20.7905813583505, eps);
  validateNumericalEqualsWithEps(luh[274],   23.3667373340221, eps);
  validateNumericalEqualsWithEps(luh[275],   40.2629471506230, eps);
  validateNumericalEqualsWithEps(luh[276],   39.4057215038980, eps);
  validateNumericalEqualsWithEps(luh[277],   39.5484958571730, eps);
  validateNumericalEqualsWithEps(luh[278],   39.6912702104480, eps);
  validateNumericalEqualsWithEps(luh[279],   42.3340445637230, eps);
  validateNumericalEqualsWithEps(luh[280],   40.9768189169980, eps);
  validateNumericalEqualsWithEps(luh[281],   40.1195932702730, eps);
  validateNumericalEqualsWithEps(luh[282],   40.2623676235480, eps);
  validateNumericalEqualsWithEps(luh[283],   40.4051419768230, eps);
  validateNumericalEqualsWithEps(luh[284],   43.0479163300980, eps);
  validateNumericalEqualsWithEps(luh[285],   22.7044530664099, eps);
  validateNumericalEqualsWithEps(luh[286],   21.7806090420815, eps);
  validateNumericalEqualsWithEps(luh[287],   21.8567650177531, eps);
  validateNumericalEqualsWithEps(luh[288],   21.9329209934247, eps);
  validateNumericalEqualsWithEps(luh[289],   24.5090769690963, eps);
  validateNumericalEqualsWithEps(luh[290],   23.0852329447679, eps);
  validateNumericalEqualsWithEps(luh[291],   22.1613889204396, eps);
  validateNumericalEqualsWithEps(luh[292],   22.2375448961112, eps);
  validateNumericalEqualsWithEps(luh[293],   22.3137008717828, eps);
  validateNumericalEqualsWithEps(luh[294],   24.8898568474544, eps);
  validateNumericalEqualsWithEps(luh[295],   43.1184342161229, eps);
  validateNumericalEqualsWithEps(luh[296],   42.2612085693979, eps);
  validateNumericalEqualsWithEps(luh[297],   42.4039829226729, eps);
  validateNumericalEqualsWithEps(luh[298],   42.5467572759479, eps);
  validateNumericalEqualsWithEps(luh[299],   45.1895316292229, eps);
  validateNumericalEqualsWithEps(luh[300],   81.3003931339190, eps);
  validateNumericalEqualsWithEps(luh[301],   80.5680611110320, eps);
  validateNumericalEqualsWithEps(luh[302],   80.8357290881451, eps);
  validateNumericalEqualsWithEps(luh[303],   81.1033970652582, eps);
  validateNumericalEqualsWithEps(luh[304],   83.8710650423712, eps);
  validateNumericalEqualsWithEps(luh[305],   44.5461777488728, eps);
  validateNumericalEqualsWithEps(luh[306],   43.6889521021478, eps);
  validateNumericalEqualsWithEps(luh[307],   43.8317264554228, eps);
  validateNumericalEqualsWithEps(luh[308],   43.9745008086978, eps);
  validateNumericalEqualsWithEps(luh[309],   46.6172751619728, eps);
  validateNumericalEqualsWithEps(luh[310],   45.2600495152478, eps);
  validateNumericalEqualsWithEps(luh[311],   44.4028238685228, eps);
  validateNumericalEqualsWithEps(luh[312],   44.5455982217978, eps);
  validateNumericalEqualsWithEps(luh[313],   44.6883725750728, eps);
  validateNumericalEqualsWithEps(luh[314],   47.3311469283478, eps);
  validateNumericalEqualsWithEps(luh[315],   85.3154127906149, eps);
  validateNumericalEqualsWithEps(luh[316],   84.5830807677280, eps);
  validateNumericalEqualsWithEps(luh[317],   84.8507487448410, eps);
  validateNumericalEqualsWithEps(luh[318],   85.1184167219541, eps);
  validateNumericalEqualsWithEps(luh[319],   87.8860846990672, eps);

  delete[] luh;
  delete[] lduh;
} // testSolutionUpdate3d

#endif // DIMENSIONS==3

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",on)
#endif
//#else
  //todo VV TestCase
//#endif
