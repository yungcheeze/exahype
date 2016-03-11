#include "EulerFlow/problem/Problem.h"

#include "EulerFlow/Constants.h"

#include "string.h"

#include "cmath"

// UNCOMMENT FOR DEBUGGING PURPOSES
// void exahype::problem::PDEInitialValue2d(const double x,const double y,const
// int nvar,double * value) {
//  for (int n=0; n < 5; n++) {
//    value[n] = 0;
//  }
//  value[0] = 1. + 0.1*x;
//  value[1] = 1. * value[0];
//  value[4] = 1./(GAMMA-1) + 0.5 * value[0] * 1.*1.;
//}

void exahype::problem::PDEInitialValue2d(const double /*in*/ x,
                                         const double /*in*/ y,
                                         double *restrict /*out*/ value) {
  value[0] = 1.;
  value[1] = 0.;
  value[2] = 0.;
  value[3] = 0.;
  value[4] = 1. / (GAMMA - 1) +
             std::exp(-((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) /
                      (0.05 * 0.05)) *
                 1.0e-3;
}

void exahype::problem::PDEFlux(const double *restrict const /*in*/ Q,
                               double *restrict /*out*/ f,
                               double *restrict /*out*/ g) {
  const double irho = 1.0 / Q[0];
  const double p =
      (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);

  /*f[0] = Q[1];
  f[1] = irho*Q[1]*Q[1] + p;
  f[2] = irho*Q[1]*Q[2];
  f[3] = irho*Q[1]*Q[3];
  f[4] = irho*Q[1]*(Q[4]+p);

  g[0] = Q[2];
  g[1] = irho*Q[2]*Q[1];
  g[2] = irho*Q[2]*Q[2] + p;
  g[3] = irho*Q[2]*Q[3];
  g[4] = irho*Q[2]*(Q[4]+p);*/

  f[0] = Q[1];
  for (int i = 1; i < 5; i++) {
    f[i] = irho * Q[1];
  }
  f[1] = f[1] * Q[1] + p;
  f[2] = f[2] * Q[2];
  f[3] = f[3] * Q[3];
  f[4] = f[4] * (Q[4] + p);

  g[0] = Q[2];
  for (int i = 1; i < 5; i++) {
    g[i] = irho * Q[2];
  }
  g[1] = g[1] * Q[1];
  g[2] = g[2] * Q[2] + p;
  g[3] = g[3] * Q[3];
  g[4] = g[4] * (Q[4] + p);
}

void exahype::problem::PDEFlux(const double *restrict const /*in*/ Q,
                               double *restrict /*out*/ f,
                               double *restrict /*out*/ g,
                               double *restrict /*out*/ h) {
  const double irho = 1.0 / Q[0];
  const double p =
      (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);

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

void exahype::problem::PDEEigenvalues(const double *restrict const /*in*/ Q,
                                      const double *restrict const /*in*/ n,
                                      double *restrict /*out*/ lambda) {
  constexpr int d = DIMENSIONS;
  constexpr int nvar = EXAHYPE_NVARS;

  __assume_aligned(lambda, ALIGNMENT);

  double irho = 1.0 / Q[0];
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);

  double u = 0;
  for (int i = 0; i < d; i++) {
    u += n[i] * Q[i + 1] * irho;
  }

  double c = std::sqrt(GAMMA * p * irho);

  lambda[0] = u - c;

  lambda[1] = u;
  lambda[2] = u;
  lambda[3] = u;

  lambda[4] = u + c;
}

double exahype::problem::PDENormalFlux(const double x, const double y,
                                       const double nx, const double ny) {
  const double r = std::sqrt(x * x + y * y);
  const double bx = -y / r;
  const double by = x / r;

  return bx * nx + by * ny;
}

double exahype::problem::DGNormalFlux(const double x, const double y,
                                      const double nx, const double ny) {
  const double pdeNormalFlux = PDENormalFlux(x, y, nx, ny);
  if (pdeNormalFlux > 0) return pdeNormalFlux;
  return 0.;
}

void exahype::problem::DGRiemannSolver(const double x, const double y,
                                       const double nx, const double ny,
                                       double *selfFlux,
                                       double *neighbourFlux) {
  *selfFlux = DGNormalFlux(x, y, nx, ny);
  *neighbourFlux = DGNormalFlux(x, y, -nx, -ny);
}
