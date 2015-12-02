#include "EulerFlow3d/problem/Problem.h"

#include "cmath"

void exahype::problem::PDEInitialValue2d(const double x,const double y,const int nvar,double * value) {
  for (int n=0; n < 5; n++) {
    value[n] = 0;
  }
  value[0] = 1. + 0.1*x;
  value[1] = 1. * value[0];
  value[4] = 1./(GAMMA-1) + 0.5 * value[0] * 1.*1.;
}

double
exahype::problem::PDEInflow(const double x,const double y) {
  return 0;
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
