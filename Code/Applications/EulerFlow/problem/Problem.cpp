#include "EulerFlow/problem/Problem.h"

#include "cmath"

// UNCOMMENT FOR DEBUGGING PURPOSES
//void exahype::problem::PDEInitialValue2d(const double x,const double y,const int nvar,double * value) {
//  for (int n=0; n < 5; n++) {
//    value[n] = 0;
//  }
//  value[0] = 1. + 0.1*x;
//  value[1] = 1. * value[0];
//  value[4] = 1./(GAMMA-1) + 0.5 * value[0] * 1.*1.;
//}

void exahype::problem::PDEInitialValue2d(const double x,const double y,const int nvar,double * value) {
  for (int n=0; n < 5; n++) {
    value[n] = 0;
  }
  value[0] = 1.;
  value[4] = 1./(GAMMA-1) + std::exp(-((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))/(0.05*0.05)) * 1.0e-3;
}


void exahype::problem::PDEFlux(const double * const Q,const int nvar,double * f,double * g) {
  double irho = 1.0/Q[0];
  double p = (GAMMA-1)*( Q[4] - 0.5* (Q[1]*Q[1] + Q[2]*Q[2]) * irho );

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

void exahype::problem::PDEEigenvalues(const double * const Q,const int nvar,const double * const n,const int d,double * lambda) {
  double irho = 1.0/Q[0];
  double p = (GAMMA-1)*( Q[4] - 0.5* (Q[1]*Q[1] + Q[2]*Q[2]) * irho );

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


double
exahype::problem::PDENormalFlux(const double x,const double y,const double nx,const double ny) {
  const double r          = std::sqrt(x*x+y*y);
  const double bx         = -y/r;
  const double by         = x/r;

  return bx*nx + by*ny;
}

double
exahype::problem::DGNormalFlux(const double x,const double y,const double nx,const double ny) {
  const double pdeNormalFlux = PDENormalFlux(x,y,nx,ny);
  if (pdeNormalFlux > 0)
    return pdeNormalFlux;
  return 0.;
}

void
exahype::problem::DGRiemannSolver(
         const double x,const double y,
         const double nx,const double ny,
         double *selfFlux,double *neighbourFlux) {
  *selfFlux      = DGNormalFlux(x,y, nx, ny);
  *neighbourFlux = DGNormalFlux(x,y,-nx,-ny);
}
