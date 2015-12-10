#include "EulerFlow3d/dg/ADERDG.h"

#include "math.h"
#include "stdlib.h"

#include "EulerFlow3d/Constants.h"
#include "EulerFlow3d/problem/Problem.h"

// 3D specialisation
template <>
double exahype::dg::stableTimeStepSize<3>(
    const double * const luh,
    const double * const dx,
    double * lambda,
    const int nvar,
    const int basisSize) {
  constexpr int dim = 3;

  // todo insert your code here
}

// 2D specialisation
template <>
double exahype::dg::stableTimeStepSize<2>(
    const double * const luh,
    const double * const dx,
    double * lambda,
    const int nvar,
    const int basisSize) {
  constexpr int dim = 2;

  // todo insert your code here

  const double normal[dim][dim]= {
      { 1., 0.},
      { 0., 1.},
  };

  double dt=1e20;
  for(int ii=0; ii < basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {
      const int nodeIndex     = ii + basisSize * jj;
      const int dofStartIndex = nodeIndex * nvar;

      double denominator=0.0;
      for (int d=0; d<dim; d++) {
        problem::PDEEigenvalues(&luh[dofStartIndex],nvar,normal[d],dim,lambda);

        double maxEigenvalue=0.0;
        for (int ivar=0; ivar<nvar; ivar++) {
          maxEigenvalue = std::max(fabs(lambda[ivar]),maxEigenvalue);
        }
        denominator += maxEigenvalue/dx[d];
      }
      dt = std::min(dt,EXAHYPE_CFL_FACTOR*PNPM[basisSize-1]/denominator); // order N = basisSize-1
    }
  }
  return dt;
}
