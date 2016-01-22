#include "exahype/aderdg/ADERDG.h"

#include "math.h"
#include "stdlib.h"

#include "limits"

#include "exahype/Constants.h"

#include "kernels/aderdg/default/PDEFluxes.h"

// 3D specialisation
template <>
double exahype::aderdg::stableTimeStepSize<3>(
    const double * const luh,
    const double * const dx,
    double * lambda,
    const int nvar,
    const int basisSize) {
  constexpr int DIM3=3;
  // todo insert your code here
}

// 2D specialisation
template <>
double exahype::aderdg::stableTimeStepSize<2>(
    const double * const luh,
    const double * const dx,
    double * lambda,
    const int nvar,
    const int basisSize) {
  constexpr int DIM2=2;

  const double normal[DIM2][DIM2]= {
      { 1., 0.},
      { 0., 1.},
  };

  double dt=std::numeric_limits<double>::max();
  for(int ii=0; ii < basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {
      const int nodeIndex     = ii + basisSize * jj;
      const int dofStartIndex = nodeIndex * nvar;

      double denominator=0.0;
      for (int d=0; d<DIM2; d++) {
        pde::PDEEigenvalues2d(&luh[dofStartIndex],nvar,normal[d],DIM2,lambda);

        double maxEigenvalue=0.0;
        for (int ivar=0; ivar<nvar; ivar++) {
          maxEigenvalue = std::max(fabs(lambda[ivar]),maxEigenvalue);
        }
        denominator += maxEigenvalue/dx[d];
      }
      dt = std::min(dt,EXAHYPE_CFL_FACTOR*PNPM[basisSize-1]/denominator); // order N = basisSize-1
    }
  }
  //return dt;
  return 0.1; // @Dominic Warum?
}
