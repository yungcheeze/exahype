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

#include "exahype/tests/kernels/c/GenericEulerKernelTest.h"

#include "../testdata/generic_euler_testdata.h"
#include "kernels/aderdg/generic/Kernels.h"

using std::cout;
using std::endl;

#if DIMENSIONS == 3

namespace exahype {
namespace tests {
namespace c {

void GenericEulerKernelTest::testFlux(const double *const Q, double **F) {
  double *f = F[0];
  double *g = F[1];
  double *h = F[2];

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

void GenericEulerKernelTest::testEigenvalues(const double *const Q,
                                             const int normalNonZeroIndex,
                                             double *lambda) {
  const double GAMMA = 1.4;

  double irho = 1.0 / Q[0];
  double p = (GAMMA - 1) *
             (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * irho);

  double u_n = Q[normalNonZeroIndex + 1] * irho;
  double c = std::sqrt(GAMMA * p * irho);

  lambda[0] = u_n - c;
  lambda[1] = u_n;
  lambda[2] = u_n;
  lambda[3] = u_n;
  lambda[4] = u_n + c;
}

void GenericEulerKernelTest::testNCP(const double *const Q,
                                     const double *const gradQ,
                                     double *BgradQ) {
  // 3D compressible Euler equations
  std::memset(BgradQ, 0, 3 * 5 * sizeof(double));
}  // testNCP

void GenericEulerKernelTest::testMatrixB(const double *const Q,
                                         const int normalNonZero, double *Bn) {
  // 3D compressible Euler equations
  double *B1 = new double[5 * 5];
  double *B2 = new double[5 * 5];
  double *B3 = new double[5 * 5];

  std::memset(B1, 0, 5 * 5 * sizeof(double));
  std::memset(B2, 0, 5 * 5 * sizeof(double));
  std::memset(B3, 0, 5 * 5 * sizeof(double));

  // Bn = B1 if normalNonZero == 0
  //      B2 if normalNonZero == 1
  //      B3 if normalNonZero == 2
  std::memcpy(Bn, (normalNonZero == 0) ? B1 : (normalNonZero == 1) ? B2 : B3,
              5 * 5 * sizeof(double));

  delete[] B1;
  delete[] B2;
  delete[] B3;
}  // testMatrixB

void GenericEulerKernelTest::testPDEFluxes() {
  cout << "Test PDE-related functions, DIM=3" << endl;

  double Q[5] = {1., 0.1, 0.2, 0.3, 3.5};  // pressure = 1.372
  double f[5], g[5], h[5];
  double *F[3] = {f, g, h};

  testFlux(Q, F);

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

void GenericEulerKernelTest::testVolumeIntegralLinear() {
  cout << "Test volume integral linear, ORDER=3, DIM=3" << endl;

  // output:
  double *lduh = new double[320];  // intentionally left uninitialised

  // input:
  const tarch::la::Vector<DIMENSIONS, double> dx(0.5, 0.5,
                                                 0.5);  // mesh spacing
  double *lFhi = new double[960]();  // nVar * nDOFx * nDOFy * nDOFz * dim
  // lFhi = [ lFhi_x  | lFhi_y | lFhi_z ], 320 entries each
  double *lFhi_x = &lFhi[0];
  double *lFhi_y = &lFhi[320];
  double *lFhi_z = &lFhi[640];

  // seed direction
  for (int i = 0; i < 320; i += 5) {
    lFhi_x[i + 1] = 1.;
    lFhi_y[i + 2] = 2.;
    lFhi_z[i + 3] = 3.;
  }

  kernels::aderdg::generic::c::volumeIntegralLinear(lduh, lFhi, dx,
                                                    5,  // numberOfVariables
                                                    4   // basisSize
                                                    );

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lduh[i], ::exahype::tests::testdata::generic_euler::
                     testVolumeIntegralLinear::lduh[i],
        eps, i);
  }

  delete[] lduh;
  delete[] lFhi;
}  // testVolumeIntegralLinear

void GenericEulerKernelTest::testVolumeIntegralNonlinear() {
  cout << "Test volume integral nonlinear, ORDER=3, DIM=3" << endl;

  // output:
  double *lduh = new double[320];  // intentionally left uninitialised

  // input:
  const double dx[3] = {0.05, 0.05, 0.05};  // mesh spacing
  double *lFhi = new double[960]();  // nVar * nDOFx * nDOFy * nDOFz * dim
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
        lduh[i], ::exahype::tests::testdata::generic_euler::
                     testVolumeIntegralNonlinear::lduh[i],
        eps, i);
  }

  delete[] lduh;
  delete[] lFhi;
}  // testVolumeIntegralNonlinear

void GenericEulerKernelTest::testSurfaceIntegralLinear() {
  cout << "Test surface integral linear, ORDER=3, DIM=3" << endl;

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
  kernels::aderdg::generic::c::surfaceIntegralLinear(
      lduh, lFhbnd, dx[0],
      5,   // getNumberOfVariables(),
      4);  // getNodesPerCoordinateAxis()

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lduh[i], ::exahype::tests::testdata::generic_euler::
                     testSurfaceIntegralLinear::lduh[i],
        eps, i);
  }

  delete[] lFhbnd;
  delete[] lduh;
}  // testSurfaceIntegralLinear

void GenericEulerKernelTest::testSurfaceIntegralNonlinear() {
  cout << "Test surface integral nonlinear, ORDER=3, DIM=3" << endl;

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
        lduh[i], ::exahype::tests::testdata::generic_euler::
                     testSurfaceIntegralNonlinear::lduh[i],
        eps, i);
  }

  delete[] lFhbnd;
  delete[] lduh;
}  // testSurfaceIntegralNonlinear

void GenericEulerKernelTest::testRiemannSolverLinear() {
  cout << "Test Riemann solver linear, ORDER=3, DIM=3" << endl;

  // output (intentionally left uninitialised):
  double *FL = new double[4 * 4 * 5];  // nDOF(3) * nDOF(2) * nVar
  double *FR = new double[4 * 4 * 5];  // nDOF(3) * nDOF(2) * nVar

  // inputs:
  // exahype::tests::testdata::generic_euler::testRiemannSolver::QL[80 =
  // nVar * nVar * nDOF]
  // exahype::tests::testdata::generic_euler::testRiemannSolver::QR[80 =
  // nVar * nVar * nDOF]
  const double dt = 1.40831757919882352703e-03;

  kernels::aderdg::generic::c::riemannSolverLinear<testEigenvalues,
                                                   testMatrixB>(
      FL, FR, exahype::tests::testdata::generic_euler::testRiemannSolver::QL,
      exahype::tests::testdata::generic_euler::testRiemannSolver::QR, dt,
      1 /*normalNonZero (only changes result of testEigenvalues, testMatrixB) */,
      5 /* numberOfVariables */, 0 /* numberOfParameters */, 4 /* basisSize */
      );

  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        FL[i], ::exahype::tests::testdata::generic_euler::
                   testRiemannSolverLinear::FL[i],
        eps, i);
  }

  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        FR[i], ::exahype::tests::testdata::generic_euler::
                   testRiemannSolverLinear::FR[i],
        eps, i);
  }

  delete[] FL;
  delete[] FR;
}  // testRiemannSolverLinear

void GenericEulerKernelTest::testRiemannSolverNonlinear() {
  cout << "Test Riemann solver nonlinear, ORDER=3, DIM=3" << endl;

  // inout:
  double *FL = new double[4 * 4 * 5];  // nDOF(3) * nDOF(2) * nVar
  double *FR = new double[4 * 4 * 5];  // nDOF(3) * nDOF(2) * nVar
  for (int i = 0; i < 80; i++) {
    // arbitrary values
    FL[i] = static_cast<double>(i + 1);
    FR[i] = static_cast<double>(i - 1);
  }

  // inputs:
  // exahype::tests::testdata::generic_euler::testRiemannSolver::QL[80 =
  // nVar * nVar * nDOF]
  // exahype::tests::testdata::generic_euler::testRiemannSolver::QR[80 =
  // nVar * nVar * nDOF]
  const double dt = 1.40831757919882352703e-03;

  kernels::aderdg::generic::c::riemannSolverNonlinear<testEigenvalues>(
      FL, FR, ::exahype::tests::testdata::generic_euler::testRiemannSolver::QL,
      ::exahype::tests::testdata::generic_euler::testRiemannSolver::QR, dt,
      1,  // normalNonZero
      5,  // numberOfVariables
      0,  // numberOfParameters
      4   // basisSize
      );

  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        FL[i], ::exahype::tests::testdata::generic_euler::
                   testRiemannSolverNonlinear::FL[i],
        eps, i);
  }

  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        FR[i], ::exahype::tests::testdata::generic_euler::
                   testRiemannSolverNonlinear::FR[i],
        eps, i);
  }

  delete[] FL;
  delete[] FR;
}  // testRiemannSolverNonlinear

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
                                              5,  // getNumberOfVariables()
                                              0,  // getNumberOfParameters()
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

void GenericEulerKernelTest::testSpaceTimePredictorLinear() {
  cout << "Test space time predictor linear, ORDER=3, DIM=3" << endl;

  // inputs:
  // exahype::tests::testdata::generic_euler::testSpaceTimePredictor::luh[320 =
  // nVar * nDOFx * nDOFy * nDOFz]

  const tarch::la::Vector<DIMENSIONS, double> dx(0.5, 0.5, 0.5);
  const double timeStepSize = 1.267423918681417E-002;

  // local:
  double *lQi = new double[1280];  // nVar * nDOFx * nDOFy * nDOFz * nDOFt
  double *lFi = new double[3840];  // nVar * nDOFx * nDOFy * nDOFz * nDOFt * dim

  // outputs:
  double *lQhi = new double[320];    // nVar * nDOFx * nDOFy * nDOFz
                                     // intentionally left uninitialised
  double *lFhi = new double[960];    // nVar * nDOFx * nDOFy * nDOFz * dim
  double *lQhbnd = new double[480];  // nVar * nDOFy * nDOF_z * 6
  double *lFhbnd = new double[480];  // nVar * nDOFy * nDOF_z * 6

  kernels::aderdg::generic::c::spaceTimePredictorLinear<testNCP>(
      lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd,
      exahype::tests::testdata::generic_euler::testSpaceTimePredictor::luh, dx,
      timeStepSize,
      5,  // numberOfVariables
      4   // basisSize
      );

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lQhi[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictorLinear::lQhi[i],
        eps, i);
  }

  for (int i = 0; i < 960; i++) {
    validateNumericalEqualsWithEpsWithParams1(lFhi[i], 0.0, eps, i);
  }

  for (int i = 0; i < 480; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lQhbnd[i], ::exahype::tests::testdata::generic_euler::
                       testSpaceTimePredictorLinear::lQhbnd[i],
        eps, i);
  }

  for (int i = 0; i < 480; i++) {
    validateNumericalEqualsWithEpsWithParams1(lFhbnd[i], 0.0, eps, i);
  }

  delete[] lQi;
  delete[] lFi;
  delete[] lFhi;
  delete[] lQhi;
  delete[] lQhbnd;
  delete[] lFhbnd;
}  // testSpaceTimePredictorLinear

void GenericEulerKernelTest::testSpaceTimePredictorNonlinear() {
  cout << "Test space time predictor nonlinear, ORDER=3, DIM=3" << endl;

  // inputs:
  // exahype::tests::testdata::generic_euler::testSpaceTimePredictor::luh[320 =
  // nVar * nDOFx * nDOFy * nDOFz]

  const tarch::la::Vector<DIMENSIONS, double> dx(0.5, 0.5, 0.5);
  const double timeStepSize = 1.267423918681417E-002;

  // local:
  double *lQi = new double[1280];  // nVar * nDOFx * nDOFy * nDOFz * nDOFt
  double *lFi = new double[3840];  // nVar * nDOFx * nDOFy * nDOFz * nDOFt * dim

  // outputs:
  double *lQhi = new double[320];    // nVar * nDOFx * nDOFy * nDOFz
                                     // intentionally left uninitialised
  double *lFhi = new double[960];    // nVar * nDOFx * nDOFy * nDOFz * dim
  double *lQhbnd = new double[480];  // nVar * nDOFy * nDOF_z * 6
  double *lFhbnd = new double[480];  // nVar * nDOFy * nDOF_z * 6

  kernels::aderdg::generic::c::spaceTimePredictorNonlinear<testFlux>(
      lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd,
      exahype::tests::testdata::generic_euler::testSpaceTimePredictor::luh, dx,
      timeStepSize,
      5,  // numberOfVariables
      4   // basisSize
      );

  for (int i = 0; i < 1280; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lQi[i], ::exahype::tests::testdata::generic_euler::
                    testSpaceTimePredictorNonlinear::lQi[i],
        5e-7, i);
  }

  for (int i = 0; i < 3840; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lFi[i], ::exahype::tests::testdata::generic_euler::
                    testSpaceTimePredictorNonlinear::lFi[i],
        2.1e-5, i);
  }

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lQhi[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictorNonlinear::lQhi[i],
        8e-8, i);
  }

  for (int i = 0; i < 960; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lFhi[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictorNonlinear::lFhi[i],
        3.4e-6, i);
  }

  for (int i = 0; i < 480; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lQhbnd[i], ::exahype::tests::testdata::generic_euler::
                       testSpaceTimePredictorNonlinear::lQhbnd[i],
        1.1e-7, i);
  }

  for (int i = 0; i < 480; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lFhbnd[i], ::exahype::tests::testdata::generic_euler::
                       testSpaceTimePredictorNonlinear::lFhbnd[i],
        4.2e-6, i);
  }

  delete[] lQi;
  delete[] lFi;
  delete[] lFhi;
  delete[] lQhi;
  delete[] lQhbnd;
  delete[] lFhbnd;
}  // testSpaceTimePredictorNonlinear

void GenericEulerKernelTest::testFaceUnknownsProjection() {
  cout << "Test face unknowns projection operators, ORDER=3, DIM=3" << endl;

  const int numberOfVariables = 1;
  const int basisSize = 4;

  // in/out
  double *lQhbndCoarseOut =
      new double[basisSize * basisSize * numberOfVariables];
  double *lFhbndCoarseOut =
      new double[basisSize * basisSize * numberOfVariables];
  double *lQhbndFineOut = new double[basisSize * basisSize * numberOfVariables];
  double *lFhbndFineOut = new double[basisSize * basisSize * numberOfVariables];

  // in:
  double *lQhbndCoarseIn =
      new double[basisSize * basisSize * numberOfVariables];
  double *lFhbndCoarseIn =
      new double[basisSize * basisSize * numberOfVariables];
  double *lQhbndFineIn = new double[basisSize * basisSize * numberOfVariables];
  double *lFhbndFineIn = new double[basisSize * basisSize * numberOfVariables];

  // Initialise to constant value.
  for (int i = 0; i < basisSize * basisSize * numberOfVariables; ++i) {
    lQhbndCoarseIn[i] = 1.0;
    lFhbndCoarseIn[i] = 1.0;
    lQhbndFineIn[i] = 1.0;
    lFhbndFineIn[i] = 1.0;
  }

  // Test the prolongation operator.
  for (int levelCoarse = 0; levelCoarse < 3; ++levelCoarse) {
    for (int levelDelta = 1; levelDelta < 1; ++levelDelta) {
      // todo For a levelDelta >= 2, assertionNumericalEquals
      // fails since the restriction result is not precise enough anymore.
      const int numberOfSubIntervals = tarch::la::aPowI(levelDelta, 3);

      tarch::la::Vector<DIMENSIONS - 1, int> subfaceIndex(0);

      // Test the restriction operator.
      memset(lQhbndCoarseOut, 0,
             basisSize * basisSize * numberOfVariables * sizeof(double));
      memset(lFhbndCoarseOut, 0,
             basisSize * basisSize * numberOfVariables * sizeof(double));
      for (int i2 = 0; i2 < numberOfSubIntervals; ++i2) {
        for (int i1 = 0; i1 < numberOfSubIntervals; ++i1) {
          subfaceIndex[0] = i1;
          subfaceIndex[1] = i2;

          // Prolongate.
          kernels::aderdg::generic::c::faceUnknownsProlongation(
              lQhbndFineOut, lFhbndFineOut, lQhbndCoarseIn, lFhbndCoarseIn,
              levelCoarse, levelCoarse + levelDelta, subfaceIndex,
              numberOfVariables, basisSize);

          // Test prolongated values.
          for (int m = 0; m < basisSize * basisSize * numberOfVariables; ++m) {
            assertionNumericalEquals5(lQhbndFineOut[m], lQhbndFineIn[m], m,
                                      levelCoarse, levelDelta, i1, i2);
            assertionNumericalEquals5(lFhbndFineOut[m], lFhbndFineIn[m], m,
                                      levelCoarse, levelDelta, i1, i2);
          }

          // Restrict.
          kernels::aderdg::generic::c::faceUnknownsRestriction(
              lQhbndCoarseOut, lFhbndCoarseOut, lQhbndFineIn, lFhbndFineIn,
              levelCoarse, levelCoarse + levelDelta, subfaceIndex,
              numberOfVariables, basisSize);
        }
        // Test restricted values.
        for (int m = 0; m < basisSize * basisSize * numberOfVariables; ++m) {
          assertionNumericalEquals3(lQhbndCoarseOut[m], lQhbndCoarseIn[m], m,
                                    levelCoarse, levelDelta);
          assertionNumericalEquals3(lFhbndCoarseOut[m], lFhbndCoarseIn[m], m,
                                    levelCoarse, levelDelta);
        }
      }
    }
  }

  delete[] lQhbndCoarseOut;
  delete[] lFhbndCoarseOut;
  delete[] lQhbndFineOut;
  delete[] lFhbndFineOut;

  delete[] lQhbndCoarseIn;
  delete[] lFhbndCoarseIn;
  delete[] lQhbndFineIn;
  delete[] lFhbndFineIn;
}

void GenericEulerKernelTest::testVolumeUnknownsProjection() {
  cout << "Test volume unknowns projection operators, ORDER=3, DIM=3" << endl;

  const int numberOfVariables = 1;
  const int basisSize = 4;

  // in/out
  double *luhCoarseOut =
      new double[basisSize * basisSize * basisSize * numberOfVariables];
  double *luhFineOut =
      new double[basisSize * basisSize * basisSize * numberOfVariables];

  // in:
  double *luhCoarseIn =
      new double[basisSize * basisSize * basisSize * numberOfVariables];
  double *luhFineIn =
      new double[basisSize * basisSize * basisSize * numberOfVariables];

  // Initialise to constant value.
  for (int i = 0; i < basisSize * basisSize * basisSize * numberOfVariables;
       ++i) {
    luhCoarseIn[i] = 1.0;
    luhFineIn[i] = 1.0;
  }

  // Test the prolongation operator.
  for (int levelCoarse = 0; levelCoarse < 3; ++levelCoarse) {
    for (int levelDelta = 1; levelDelta < 2; ++levelDelta) {
      // todo For a levelDelta >= 2, assertionNumericalEquals
      // fails since the prolongation result is not precise enough anymore.
      const int numberOfSubIntervals = tarch::la::aPowI(levelDelta, 3);

      // Test the restriction operator.
      tarch::la::Vector<DIMENSIONS, int> subcellIndex(0);
      memset(luhCoarseOut, 0, basisSize * basisSize * basisSize *
                                  numberOfVariables * sizeof(double));

      for (int i3 = 0; i3 < numberOfSubIntervals; ++i3) {
        for (int i2 = 0; i2 < numberOfSubIntervals; ++i2) {
          for (int i1 = 0; i1 < numberOfSubIntervals; ++i1) {
            subcellIndex[0] = i1;
            subcellIndex[1] = i2;
            subcellIndex[2] = i3;

            // Prolongate.
            kernels::aderdg::generic::c::volumeUnknownsProlongation(
                luhFineOut, luhCoarseIn, levelCoarse, levelCoarse + levelDelta,
                subcellIndex, numberOfVariables, basisSize);

            // Test prolongated values.
            for (int m = 0;
                 m < basisSize * basisSize * basisSize * numberOfVariables;
                 ++m) {
              /*
              assertionNumericalEquals5(luhFineOut[m], luhFineIn[m], m,
                                        levelCoarse, levelDelta, i1, i2);
                                        */
            }

            // Restrict.
            kernels::aderdg::generic::c::volumeUnknownsRestriction(
                luhCoarseOut, luhFineIn, levelCoarse, levelCoarse + levelDelta,
                subcellIndex, numberOfVariables, basisSize);
          }
        }
      }
      // Test restricted values.
      for (int m = 0; m < basisSize * basisSize * basisSize * numberOfVariables;
           ++m) {
        assertionNumericalEquals3(luhCoarseOut[m], luhCoarseIn[m], m,
                                  levelCoarse, levelDelta);
      }
    }
  }

  delete[] luhCoarseOut;
  delete[] luhFineOut;

  delete[] luhCoarseIn;
  delete[] luhFineIn;
}

}  // namespace c
}  // namespace tests
}  // namespace exahype

#endif  // DIMENSIONS==3
