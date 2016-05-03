#include <tarch/la/Vector.h>

#include "kernels/aderdg/generic/Kernels.h"

#if DIMENSIONS == 3

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

void surfaceIntegralNonlinear(double *lduh, const double *const lFbnd,
                              const tarch::la::Vector<DIMENSIONS, double> &dx,
                              const int numberOfVariables,
                              const int basisSize) {
  const int basisSize2 = basisSize * basisSize;
  const int order = basisSize - 1;

  // x faces
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      const double weight = kernels::gaussLegendreWeights[order][i] *
                            kernels::gaussLegendreWeights[order][j];
      const double updateSize = weight / dx[0];

      // TODO(gyro): Benchmark l,k vs. k,l vs. blocking
      for (int k = 0; k < numberOfVariables; k++) {
        // left flux minus right flux
        for (int l = 0; l < basisSize; l++) {
          lduh[i * basisSize2 * numberOfVariables +
               j * basisSize * numberOfVariables + l * numberOfVariables + k] -=
              (lFbnd[1 * basisSize2 * numberOfVariables +
                     i * basisSize * numberOfVariables + j * numberOfVariables +
                     k] *
                   kernels::FRCoeff[order][l] -
               lFbnd[0 * basisSize2 * numberOfVariables +
                     i * basisSize * numberOfVariables + j * numberOfVariables +
                     k] *
                   kernels::FLCoeff[order][l]) *
              updateSize;
        }
      }
    }
  }

  // y faces
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      const double weight = kernels::gaussLegendreWeights[order][i] *
                            kernels::gaussLegendreWeights[order][j];
      const double updateSize = weight / dx[1];

      for (int k = 0; k < numberOfVariables; k++) {
        // back flux minus front flux
        for (int l = 0; l < basisSize; l++) {
          lduh[i * basisSize2 * numberOfVariables +
               l * basisSize * numberOfVariables + j * numberOfVariables + k] -=
              (lFbnd[3 * basisSize2 * numberOfVariables +
                     i * basisSize * numberOfVariables + j * numberOfVariables +
                     k] *
                   kernels::FRCoeff[order][l] -
               lFbnd[2 * basisSize2 * numberOfVariables +
                     i * basisSize * numberOfVariables + j * numberOfVariables +
                     k] *
                   kernels::FLCoeff[order][l]) *
              updateSize;
        }
      }
    }
  }

  // z faces
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      const double weight = kernels::gaussLegendreWeights[order][i] *
                            kernels::gaussLegendreWeights[order][j];
      const double updateSize = weight / dx[2];

      for (int k = 0; k < numberOfVariables; k++) {
        // bottom flux minus top flux
        for (int l = 0; l < basisSize; l++) {
          lduh[l * basisSize2 * numberOfVariables +
               i * basisSize * numberOfVariables + j * numberOfVariables + k] -=
              (lFbnd[5 * basisSize2 * numberOfVariables +
                     i * basisSize * numberOfVariables + j * numberOfVariables +
                     k] *
                   kernels::FRCoeff[order][l] -
               lFbnd[4 * basisSize2 * numberOfVariables +
                     i * basisSize * numberOfVariables + j * numberOfVariables +
                     k] *
                   kernels::FLCoeff[order][l]) *
              updateSize;
        }
      }
    }
  }
}

}  // namespace c
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels

#endif  // DIMENSIONS == 3
