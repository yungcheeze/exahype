#include "kernels/aderdg/default/PDEFluxes.h"

#include <cmath>

#define GAMMA 1.4

void exahype::pde::PDEFlux2d(const double * const Q,const int nvar,double * f,double * g) {
  double irho = 1.0/Q[0];
  double p    = (GAMMA-1)*( Q[4] - 0.5* (Q[1]*Q[1] + Q[2]*Q[2]) * irho );

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

void exahype::pde::PDEEigenvalues2d(const double * const Q,const int nvar,const double * const n,const int d,double * lambda) {
  double irho = 1.0/Q[0];
  double p    = (GAMMA-1)*( Q[4] - 0.5* (Q[1]*Q[1] + Q[2]*Q[2]) * irho );

  double u = 0;
  for (int i=0; i < d; i++) {
    u += n[i] * Q[i+1] * irho;
  }

  double c = std::sqrt(GAMMA * p * irho);

  lambda[0] = u-c;

  lambda[1] = u;
  lambda[2] = u;
  lambda[3] = u;

  lambda[4] = u+c;
}
