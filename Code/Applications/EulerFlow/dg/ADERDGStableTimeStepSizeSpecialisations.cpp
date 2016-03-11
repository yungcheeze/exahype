#include "EulerFlow/dg/ADERDG.h"

#include "math.h"
#include "stdlib.h"

#include "limits"

#include "EulerFlow/Constants.h"
#include "EulerFlow/problem/Problem.h"

// 3D specialisation
template <>
double exahype::dg::stableTimeStepSize<3>(const double* restrict const luh,
                                          const double* restrict const dx,
                                          double* restrict lambda) {
  constexpr int dim = DIMENSIONS;  // 3
  constexpr int nvar = EXAHYPE_NVARS;
  constexpr int basisSize = EXAHYPE_ORDER + 1;

  // todo insert your code here
}

// 2D specialisation
template <>
double exahype::dg::stableTimeStepSize<2>(const double* restrict const luh,
                                          const double* restrict const dx,
                                          double* restrict lambda) {
  constexpr int dim = DIMENSIONS;  // 2
  constexpr int nvar = EXAHYPE_NVARS;
  constexpr int basisSize = EXAHYPE_ORDER + 1;

  __assume_aligned(lambda, ALIGNMENT);

  const double normal[dim][dim] = {
      {1., 0.}, {0., 1.},
  };

  double dt = std::numeric_limits<double>::max();
  for (int ii = 0; ii < basisSize; ii++) {
    for (int jj = 0; jj < basisSize; jj++) {
      const int nodeIndex = ii + basisSize * jj;
      const int dofStartIndex = nodeIndex * nvar;

      double denominator = 0.0;
      for (int d = 0; d < dim; d++) {
        problem::PDEEigenvalues(&luh[dofStartIndex], normal[d], lambda);

        double maxEigenvalue = 0.0;
        for (int ivar = 0; ivar < nvar; ivar++) {
          maxEigenvalue = std::max(fabs(lambda[ivar]), maxEigenvalue);
        }
        denominator += maxEigenvalue / dx[d];
      }
      dt = std::min(dt, EXAHYPE_CFL_FACTOR * PNPM[basisSize - 1] /
                            denominator);  // order N = basisSize-1
    }
  }
  return dt;
}
