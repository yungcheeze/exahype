#include "EulerFlow3d/problem/Problem.h"

#include "cmath"

double
exahype::problem::PDEInitialValue(const double x,const double y) {
  if ((x-0.5)*(x-0.5) + (y-0.25)*(y-0.25)  < 0.01) {
    return 1.;
  }
  return 0.;
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
