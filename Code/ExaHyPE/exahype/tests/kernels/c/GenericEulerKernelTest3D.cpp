#include "exahype/tests/kernels/c/GenericEulerKernelTest.h"

#include "kernels/aderdg/generic/Kernels.h"
#include "../testdata/generic_euler_testdata.h"

using std::cout;
using std::endl;

#if DIMENSIONS == 3

namespace exahype {
namespace tests {
namespace c {

void GenericEulerKernelTest::testFlux(const double *const Q, double *f,
                                      double *g, double *h) {
  const double GAMMA = 1.4;

  const double irho = 1.0 / Q[0];
  const double p =
      (GAMMA - 1) *
      (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * irho);

  f[0] = Q[1];
  f[1] = irho * Q[1] * Q[1] + p;
  f[2] = irho * Q[1] * Q[2];
  f[3] = irho * Q[1] * Q[3];
  f[4] = irho * Q[1] * (Q[4] + p);

  g[0] = Q[2];
  g[1] = irho * Q[2] * Q[1];
  g[2] = irho * Q[2] * Q[2] + p;
  g[3] = irho * Q[2] * Q[3];
  g[4] = irho * Q[2] * (Q[4] + p);

  h[0] = Q[3];
  h[1] = irho * Q[3] * Q[1];
  h[2] = irho * Q[3] * Q[2];
  h[3] = irho * Q[3] * Q[3] + p;
  h[4] = irho * Q[3] * (Q[4] + p);
}

void GenericEulerKernelTest::testPDEFluxes() {
  cout << "Test PDE-related functions, DIM=3" << endl;

  double Q[5] = {1., 0.1, 0.2, 0.3, 3.5};  // pressure = 1.39
  double f[5], g[5], h[5];

  testFlux(Q, f, g, h);

  for (int i = 0; i < 5; i++) {
    validateNumericalEqualsWithParams1(
        f[i], ::exahype::tests::testdata::generic_euler::testPDEFluxes::f[i],
        i);
  }

  for (int i = 0; i < 5; i++) {
    validateNumericalEqualsWithParams1(
        g[i], ::exahype::tests::testdata::generic_euler::testPDEFluxes::g[i],
        i);
  }

  for (int i = 0; i < 5; i++) {
    validateNumericalEqualsWithParams1(
        h[i], ::exahype::tests::testdata::generic_euler::testPDEFluxes::h[i],
        i);
  }
}  // testPDEFluxes

void GenericEulerKernelTest::testVolumeIntegral() {
  cout << "Test volume integral, ORDER=3, DIM=3" << endl;

  // output:
  double *lduh = new double[320];  // intentionally left uninitialised

  // input:
  const double dx[3] = {0.05, 0.05, 0.05};  // mesh spacing
  double *lFhi = new double[960]();  // nVar * dim * nDOFx * nDOFy * nDOFz
  // lFhi = [ lFhi_x  | lFhi_y | lFhi_z ], 320 entries each
  double *lFhi_x = &lFhi[0];
  double *lFhi_y = &lFhi[320];
  double *lFhi_z = &lFhi[640];

  // seed direction
  for (int i = 0; i < 320; i += 5) {
    lFhi_x[i + 1] = 1.;
    lFhi_y[i + 2] = 1.;
    lFhi_z[i + 3] = 1.;
  }

  kernels::aderdg::generic::c::volumeIntegralNonlinear(
      lduh, lFhi, dx[0],
      5,   // getNumberOfVariables(),
      4);  // getNodesPerCoordinateAxis()

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lduh[i],
        ::exahype::tests::testdata::generic_euler::testVolumeIntegral::lduh[i],
        eps, i);
  }

  delete[] lduh;
  delete[] lFhi;
}  // testVolumeIntegral

void GenericEulerKernelTest::testSurfaceIntegral() {
  cout << "Test surface integral, ORDER=3, DIM=3" << endl;

  // inputs:
  const double dx[3] = {0.1, 0.1, 0.1};   // mesh spacing
  double *lFhbnd = new double[6 * 80]();  // 480
  double *FLeft = &lFhbnd[0];
  double *FRight = &lFhbnd[80];
  double *FFront = &lFhbnd[160];
  double *FBack = &lFhbnd[240];
  double *FBottom = &lFhbnd[320];
  double *FTop = &lFhbnd[400];

  for (int i = 0; i < 80; i += 5) {
    // in x orientation 1
    FLeft[i + 1] = 1.;
    FRight[i + 1] = 1.;
    // in y orientation 1
    FFront[i + 2] = 1.;
    FBack[i + 2] = 1.;
    // in z direction 1
    FBottom[i + 3] = 1.;
    FTop[i + 3] = 1.;
  }

  // in/out:
  double *lduh = new double[320];
  for (int i = 0; i < 320; i++) {
    lduh[i] = i / 10.;
  }

  // lFhbnd = [ FLeft | FRight | FFront | FBack | FBottom | FTop ]
  kernels::aderdg::generic::c::surfaceIntegralNonlinear(
      lduh, lFhbnd, dx[0],
      5,   // getNumberOfVariables(),
      4);  // getNodesPerCoordinateAxis()

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lduh[i],
        ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::lduh[i],
        eps, i);
  }

  delete[] lFhbnd;
  delete[] lduh;
}  // testSurfaceIntegral

void GenericEulerKernelTest::testRiemannSolver() {}

void GenericEulerKernelTest::testSolutionUpdate() {
  cout << "Test solution update, ORDER=3, DIM=3" << endl;

  // in/out:
  double *luh = new double[320]();
  for (int i = 0; i < 320; i += 5) {
    luh[i] = 1.0;
    luh[i + 4] = 2.5;
  }

  // inputs:
  const double dt = 1.40831757919882352703e-03;
  double *lduh = new double[320];
  for (int i = 0; i < 320; i++) {
    lduh[i] = i;
  }

  kernels::aderdg::generic::c::solutionUpdate(luh, lduh, dt,
                                              5,  // getNumberOfVariables(),
                                              4   // getNodesPerCoordinateAxis()
                                              );

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        luh[i],
        ::exahype::tests::testdata::generic_euler::testSolutionUpdate::luh[i],
        eps, i);
  }

  delete[] luh;
  delete[] lduh;
}  // testSolutionUpdate

void GenericEulerKernelTest::testSpaceTimePredictor() {
  cout << "Test space time predictor, ORDER=3, DIM=3" << endl;

  // inputs:
  double *luh = new double[320]();  // nVar * nDOFx * nDOFy * nDOFz
  for (int i = 0; i < 320; i += 5) {
    luh[i] = 1.0;
    luh[i + 4] = 2.5;
  }

  const double dx[3] = {0.5, 0.5, 0.5};
  const double timeStepSize = 1.267423918681417E-002;

  // local:
  double *lQi = new double[1280];  // nVar * nDOFx * nDOFy * nDOFz * nDOFt
  double *lFi = new double[3840];  // nVar * nDOFx * nDOFy * nDOFz * nDOFt * dim

  // outputs:
  double *lQhi = new double[320];  // intentionally left uninitialised
  double *lFhi = new double[960];  // nVar * nDOFx * nDOFy * nDOFz * dim
  double *lQhbnd = new double[480];
  double *lFhbnd = new double[480];

  kernels::aderdg::generic::fortran::spaceTimePredictorNonlinear<testFlux>(
      lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd, luh, dx[0], timeStepSize, 5, 4);

  validateNumericalEqualsWithEps(lQhi[0], 9.9999999999978384E-001, eps);
  validateNumericalEqualsWithEps(lQhi[1], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[2], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[3], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[4], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[5], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[6], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[7], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[8], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[9], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[10], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[11], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[12], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[13], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[14], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[15], 9.9999999999978384E-001, eps);
  validateNumericalEqualsWithEps(lQhi[16], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[17], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[18], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[19], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[20], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[21], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[22], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[23], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[24], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[25], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[26], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[27], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[28], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[29], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[30], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[31], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[32], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[33], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[34], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[35], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[36], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[37], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[38], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[39], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[40], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[41], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[42], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[43], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[44], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[45], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[46], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[47], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[48], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[49], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[50], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[51], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[52], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[53], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[54], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[55], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[56], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[57], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[58], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[59], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[60], 9.9999999999978384E-001, eps);
  validateNumericalEqualsWithEps(lQhi[61], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[62], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[63], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[64], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[65], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[66], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[67], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[68], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[69], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[70], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[71], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[72], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[73], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[74], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[75], 9.9999999999978384E-001, eps);
  validateNumericalEqualsWithEps(lQhi[76], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[77], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[78], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[79], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[80], 9.9999999999978351E-001, eps);
  validateNumericalEqualsWithEps(lQhi[81], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[82], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[83], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[84], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[85], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[86], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[87], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[88], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[89], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[90], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[91], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[92], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[93], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[94], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[95], 9.9999999999978351E-001, eps);
  validateNumericalEqualsWithEps(lQhi[96], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[97], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[98], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[99], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[100], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[101], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[102], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[103], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[104], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[105], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[106], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[107], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[108], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[109], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[110], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[111], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[112], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[113], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[114], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[115], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[116], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[117], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[118], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[119], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[120], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[121], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[122], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[123], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[124], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[125], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[126], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[127], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[128], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[129], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[130], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[131], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[132], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[133], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[134], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[135], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[136], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[137], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[138], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[139], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[140], 9.9999999999978351E-001, eps);
  validateNumericalEqualsWithEps(lQhi[141], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[142], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[143], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[144], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[145], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[146], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[147], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[148], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[149], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[150], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[151], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[152], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[153], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[154], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[155], 9.9999999999978351E-001, eps);
  validateNumericalEqualsWithEps(lQhi[156], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[157], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[158], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[159], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[160], 9.9999999999978351E-001, eps);
  validateNumericalEqualsWithEps(lQhi[161], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[162], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[163], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[164], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[165], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[166], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[167], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[168], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[169], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[170], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[171], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[172], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[173], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[174], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[175], 9.9999999999978351E-001, eps);
  validateNumericalEqualsWithEps(lQhi[176], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[177], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[178], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[179], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[180], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[181], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[182], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[183], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[184], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[185], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[186], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[187], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[188], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[189], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[190], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[191], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[192], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[193], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[194], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[195], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[196], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[197], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[198], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[199], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[200], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[201], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[202], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[203], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[204], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[205], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[206], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[207], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[208], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[209], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[210], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[211], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[212], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[213], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[214], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[215], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[216], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[217], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[218], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[219], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[220], 9.9999999999978351E-001, eps);
  validateNumericalEqualsWithEps(lQhi[221], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[222], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[223], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[224], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[225], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[226], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[227], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[228], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[229], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[230], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[231], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[232], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[233], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[234], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[235], 9.9999999999978351E-001, eps);
  validateNumericalEqualsWithEps(lQhi[236], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[237], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[238], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[239], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[240], 9.9999999999978384E-001, eps);
  validateNumericalEqualsWithEps(lQhi[241], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[242], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[243], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[244], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[245], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[246], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[247], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[248], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[249], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[250], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[251], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[252], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[253], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[254], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[255], 9.9999999999978384E-001, eps);
  validateNumericalEqualsWithEps(lQhi[256], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[257], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[258], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[259], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[260], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[261], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[262], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[263], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[264], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[265], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[266], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[267], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[268], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[269], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[270], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[271], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[272], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[273], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[274], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[275], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[276], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[277], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[278], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[279], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[280], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[281], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[282], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[283], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[284], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[285], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[286], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[287], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[288], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[289], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[290], 9.9999999999978373E-001, eps);
  validateNumericalEqualsWithEps(lQhi[291], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[292], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[293], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[294], 2.4999999999994591E+000, eps);
  validateNumericalEqualsWithEps(lQhi[295], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[296], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[297], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[298], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[299], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[300], 9.9999999999978384E-001, eps);
  validateNumericalEqualsWithEps(lQhi[301], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[302], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[303], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[304], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[305], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[306], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[307], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[308], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[309], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[310], 9.9999999999978362E-001, eps);
  validateNumericalEqualsWithEps(lQhi[311], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[312], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[313], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[314], 2.4999999999994595E+000, eps);
  validateNumericalEqualsWithEps(lQhi[315], 9.9999999999978384E-001, eps);
  validateNumericalEqualsWithEps(lQhi[316], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[317], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[318], 0.0000000000000000E+000, eps);
  validateNumericalEqualsWithEps(lQhi[319], 2.4999999999994595E+000, eps);

  delete[] luh;
  delete[] lQi;
  delete[] lFi;
  delete[] lFhi;
  delete[] lQhi;
  delete[] lQhbnd;
  delete[] lFhbnd;
}  // testSpaceTimePredictor

}  // namespace c
}  // namespace tests
}  // namespace exahype

#endif  // DIMENSIONS==3
