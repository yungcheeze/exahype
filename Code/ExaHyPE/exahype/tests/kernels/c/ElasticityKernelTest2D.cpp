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

#include "exahype/tests/kernels/c/ElasticityKernelTest.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <numeric>

#include "../testdata/elasticity_testdata.h"
#include "kernels/KernelUtils.h"
#include "kernels/aderdg/generic/Kernels.h"

#if DIMENSIONS == 2

namespace exahype {
namespace tests {
namespace c {

static const int kNumberOfParameters = 3;
static const int kNumberOfVariables = 9 + kNumberOfParameters;
static const int kN = 4;
static const int kBasisSize = kN + 1;

void ElasticityKernelTest::testFlux(const double *Q, double **F) {
  double *f = F[0];
  double *g = F[1];
  // double* h = F[2];

  std::fill(f, f + kNumberOfVariables, 0.0);
  std::fill(g, g + kNumberOfVariables, 0.0);
  // std::fill(h, h + kNumberOfVariables, 0.0);

  double lam = Q[9];          // par(2)
  double mu = Q[10];          // par(2)
  double irho = 1.0 / Q[11];  // 1.0 / par(3)

  f[1 - 1] = -(lam + 2 * mu) * Q[7 - 1];
  f[2 - 1] = -lam * Q[7 - 1];
  f[3 - 1] = -lam * Q[7 - 1];
  f[4 - 1] = -mu * Q[8 - 1];
  f[5 - 1] = 0.0;
  f[6 - 1] = -mu * Q[9 - 1];
  f[7 - 1] = -irho * Q[1 - 1];
  f[8 - 1] = -irho * Q[4 - 1];
  f[9 - 1] = -irho * Q[6 - 1];

  g[1 - 1] = -lam * Q[8 - 1];
  g[2 - 1] = -(lam + 2 * mu) * Q[8 - 1];
  g[3 - 1] = -lam * Q[8 - 1];
  g[4 - 1] = -mu * Q[7 - 1];
  g[5 - 1] = -mu * Q[9 - 1];
  g[6 - 1] = 0.0;
  g[7 - 1] = -irho * Q[4 - 1];
  g[8 - 1] = -irho * Q[2 - 1];
  g[9 - 1] = -irho * Q[5 - 1];

  //  h[1 - 1] = - lam*Q[9 - 1];
  //  h[2 - 1] = - lam*Q[9 - 1];
  //  h[3 - 1] = - (lam+2*mu)*Q[9 - 1];
  //  h[4 - 1] = 0.0;
  //  h[5 - 1] = - mu *Q[8 - 1];
  //  h[6 - 1] = - mu *Q[7 - 1];
  //  h[7 - 1] = - irho *Q[4 - 1];
  //  h[8 - 1] = - irho *Q[2 - 1];
  //  h[9 - 1] = - irho *Q[5 - 1];
}

void ElasticityKernelTest::testSource(const double *Q, double *S) {
  std::fill(S, S + kNumberOfVariables, 0.0);
}

void ElasticityKernelTest::testEigenvalues(const double *const Q,
                                           const int normalNonZeroIndex,
                                           double *lambda) {
  std::fill(lambda, lambda + kNumberOfParameters, 0.0);

  double lam = Q[9];    // par(1)
  double mu = Q[10];    // par(2)
  double rho0 = Q[11];  // par(3)
  double cp = std::sqrt((lam + 2 * mu) / rho0);
  double cs = std::sqrt(mu / rho0);

  lambda[1 - 1] = -cp;
  lambda[2 - 1] = -cs;
  lambda[3 - 1] = -cs;
  lambda[4 - 1] = 0.0;
  lambda[5 - 1] = 0.0;
  lambda[6 - 1] = 0.0;
  lambda[7 - 1] = +cs;
  lambda[8 - 1] = +cs;
  lambda[9 - 1] = +cp;
}

void ElasticityKernelTest::testNCP(const double *const Q,
                                   const double *const gradQ, double *BgradQ) {
  std::fill(BgradQ, BgradQ + kNumberOfVariables * DIMENSIONS, 0.0);

  double lam = Q[9];          // par(1)
  double mu = Q[10];          // par(2)
  double irho = 1.0 / Q[11];  // 1.0 / par(3)

  //  const double *gradQx = gradQ + 0 * kNumberOfVariables;
  //  const double *gradQy = gradQ + 1 * kNumberOfVariables;
  //  const double *gradQz = gradQ + 2 * kNumberOfVariables;

  double *BgradQx = BgradQ + 0 * kNumberOfVariables;
  double *BgradQy = BgradQ + 1 * kNumberOfVariables;
  double *BgradQz = BgradQ + 2 * kNumberOfVariables;

  BgradQx[1 - 1] = -(lam + 2 * mu) * BgradQx[7 - 1];
  BgradQx[2 - 1] = -lam * BgradQx[7 - 1];
  BgradQx[3 - 1] = -lam * BgradQx[7 - 1];
  BgradQx[4 - 1] = -mu * BgradQx[8 - 1];
  BgradQx[5 - 1] = 0.0;
  BgradQx[6 - 1] = -mu * BgradQx[9 - 1];
  BgradQx[7 - 1] = -irho * BgradQx[1 - 1];
  BgradQx[8 - 1] = -irho * BgradQx[4 - 1];
  BgradQx[9 - 1] = -irho * BgradQx[6 - 1];

  BgradQy[1 - 1] = -lam * BgradQy[8 - 1];
  BgradQy[2 - 1] = -(lam + 2 * mu) * BgradQy[8 - 1];
  BgradQy[3 - 1] = -lam * BgradQy[8 - 1];
  BgradQy[4 - 1] = -mu * BgradQy[7 - 1];
  BgradQy[5 - 1] = -mu * BgradQy[9 - 1];
  BgradQy[6 - 1] = 0.0;
  BgradQy[7 - 1] = -irho * BgradQy[4 - 1];
  BgradQy[8 - 1] = -irho * BgradQy[2 - 1];
  BgradQy[9 - 1] = -irho * BgradQy[5 - 1];

  BgradQy[1 - 1] = -lam * BgradQz[9 - 1];
  BgradQy[2 - 1] = -lam * BgradQz[9 - 1];
  BgradQy[3 - 1] = -(lam + 2 * mu) * BgradQz[9 - 1];
  BgradQy[4 - 1] = 0.0;
  BgradQy[5 - 1] = -mu * BgradQz[8 - 1];
  BgradQy[6 - 1] = -mu * BgradQz[7 - 1];
  BgradQy[7 - 1] = -irho * BgradQz[6 - 1];
  BgradQy[8 - 1] = -irho * BgradQz[5 - 1];
  BgradQy[9 - 1] = -irho * BgradQz[3 - 1];

}  // testNCP

void ElasticityKernelTest::testMatrixB(const double *const Q,
                                       const int normalNonZero, double *Bn) {
  std::fill(Bn, Bn + kNumberOfVariables * kNumberOfVariables, 0.0);

  kernels::idx2 idx_Bn(kNumberOfVariables, kNumberOfVariables);

  double lam = Q[9];          // par(1)
  double mu = Q[10];          // par(2)
  double irho = 1.0 / Q[11];  // 1./par(3)

  switch (normalNonZero) {
    case 0:
      Bn[idx_Bn(7 - 1, 1 - 1)] = -(lam + 2 * mu);
      Bn[idx_Bn(7 - 1, 2 - 1)] = -lam;
      Bn[idx_Bn(7 - 1, 3 - 1)] = -lam;
      Bn[idx_Bn(8 - 1, 4 - 1)] = -mu;
      Bn[idx_Bn(9 - 1, 6 - 1)] = -mu;
      Bn[idx_Bn(1 - 1, 7 - 1)] = -irho;
      Bn[idx_Bn(4 - 1, 8 - 1)] = -irho;
      Bn[idx_Bn(6 - 1, 9 - 1)] = -irho;
      break;
    case 1:
      Bn[idx_Bn(8 - 1, 1 - 1)] = -lam;
      Bn[idx_Bn(8 - 1, 2 - 1)] = -(lam + 2 * mu);
      Bn[idx_Bn(8 - 1, 3 - 1)] = -lam;
      Bn[idx_Bn(7 - 1, 4 - 1)] = -mu;
      Bn[idx_Bn(9 - 1, 5 - 1)] = -mu;
      Bn[idx_Bn(4 - 1, 7 - 1)] = -irho;
      Bn[idx_Bn(2 - 1, 8 - 1)] = -irho;
      Bn[idx_Bn(5 - 1, 9 - 1)] = -irho;
      break;
    //    case 2:
    //      Bn[idx_Bn(9 - 1, 1 - 1)] = -lam;
    //      Bn[idx_Bn(9 - 1, 2 - 1)] = -lam;
    //      Bn[idx_Bn(9 - 1, 3 - 1)] = -(lam + 2 * mu);
    //      Bn[idx_Bn(8 - 1, 5 - 1)] = -mu;
    //      Bn[idx_Bn(7 - 1, 6 - 1)] = -mu;
    //      Bn[idx_Bn(6 - 1, 7 - 1)] = -irho;
    //      Bn[idx_Bn(5 - 1, 8 - 1)] = -irho;
    //      Bn[idx_Bn(3 - 1, 9 - 1)] = -irho;
    //      break;
    default:
      assert(false);
      break;
  }
}  // testMatrixB

void ElasticityKernelTest::testRiemannSolverLinear() {
  logInfo("ElasticityKernelTest::testRiemannSolverLinear()",
          "Test Riemann solver linear, ORDER=4, DIM=2");

  double *qL = new double[kBasisSize * kNumberOfVariables];
  double *qR = new double[kBasisSize * kNumberOfVariables];

  kernels::idx2 idx_q(kBasisSize, kNumberOfVariables);
  kernels::idx2 idx_q_in(kBasisSize, kNumberOfVariables - kNumberOfParameters);
  kernels::idx2 idx_param_in(kBasisSize, kNumberOfParameters);

  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kNumberOfVariables - kNumberOfParameters; j++) {
      qL[idx_q(i, j)] = exahype::tests::testdata::elasticity::
          testRiemannSolverLinear::qL_IN[idx_q_in(i, j)];
      qR[idx_q(i, j)] = exahype::tests::testdata::elasticity::
          testRiemannSolverLinear::qR_IN[idx_q_in(i, j)];
    }

    for (int j = 0; j < kNumberOfParameters; j++) {
      qL[idx_q(i, j + kNumberOfVariables - kNumberOfParameters)] =
          exahype::tests::testdata::elasticity::testRiemannSolverLinear::
              paramL_IN[idx_param_in(i, j)];
      qR[idx_q(i, j + kNumberOfVariables - kNumberOfParameters)] =
          exahype::tests::testdata::elasticity::testRiemannSolverLinear::
              paramR_IN[idx_param_in(i, j)];
    }
  }

  double *FL = new double[kBasisSize * kNumberOfVariables];
  double *FR = new double[kBasisSize * kNumberOfVariables];
  std::fill(FL, FL + kBasisSize * kNumberOfVariables,
            std::numeric_limits<double>::quiet_NaN());
  std::fill(FR, FR + kBasisSize * kNumberOfVariables,
            std::numeric_limits<double>::quiet_NaN());

  kernels::idx2 idx_F(kBasisSize, kNumberOfVariables);
  kernels::idx2 idx_F_in(kBasisSize, kNumberOfVariables - kNumberOfParameters);

  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kNumberOfVariables - kNumberOfParameters; j++) {
      FL[idx_F(i, j)] = exahype::tests::testdata::elasticity::
          testRiemannSolverLinear::FL_IN[idx_F_in(i, j)];
      FR[idx_F(i, j)] = exahype::tests::testdata::elasticity::
          testRiemannSolverLinear::FR_IN[idx_F_in(i, j)];
    }
  }

  const double dt = 1.916666666666667E-004;

  kernels::aderdg::generic::c::riemannSolverLinear<testEigenvalues,
                                                   testMatrixB>(
      FL, FR, qL, qR, dt, 1 /* normalNonZero */, kNumberOfVariables,
      kNumberOfParameters, kBasisSize);

  kernels::idx2 idx_F_out(kBasisSize, kNumberOfVariables - kNumberOfParameters);

  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kNumberOfVariables - kNumberOfParameters; j++) {
      validateNumericalEqualsWithEpsWithParams1(
          FL[idx_F(i, j)], exahype::tests::testdata::elasticity::
                               testRiemannSolverLinear::FL_OUT[idx_F_out(i, j)],
          eps, idx_F_out(i, j));
    }
  }

  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kNumberOfVariables - kNumberOfParameters; j++) {
      validateNumericalEqualsWithEpsWithParams1(
          FR[idx_F(i, j)], exahype::tests::testdata::elasticity::
                               testRiemannSolverLinear::FR_OUT[idx_F_out(i, j)],
          eps, idx_F_out(i, j));
    }
  }

  delete[] qL;
  delete[] qR;
  delete[] FL;
  delete[] FR;
}

}  // namespace c
}  // namespace tests
}  // namespace exahype

#endif  // DIMENSIONS==2
