#include "EulerFlow3d/problem/Problem.h"

#include "cmath"

void exahype::problem::PDEInitialValue2d(const double x,const double y,const int numberOfVariables,double * value) {
  for (int n=0; n < 5; n++) {
    value[n] = 0;
  }
  value[0] = 1. + 0.1*x + 0.2*y;
  value[4] = 1;

//  for (int n=0; n < 5; n++) {
//    value[n] = n+1;
//  }
}

double
exahype::problem::PDEInflow(const double x,const double y) {
  return 0;
}

double
exahype::problem::PDEVolumeFlux(const double x, const double y) {
  return 0.;
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
