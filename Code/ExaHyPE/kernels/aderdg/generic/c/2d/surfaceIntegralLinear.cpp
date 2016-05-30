#include "../../Kernels.h"

#if DIMENSIONS == 2

// WARNING: Untested!

void surfaceIntegralLinear(double* lduh, const double* const lFbnd,
                           const tarch::la::Vector<DIMENSIONS, double>& dx,
                           const int numberOfVariables, const int basisSize) {
  const int order = basisSize - 1;

  // x faces
  for (int i = 0; i < basisSize; i++) {
    const double weight = kernels::gaussLegendreWeights[order][i];
    const double updateSize = weight / dx[0];

    for (int j = 0; j < basisSize; j++) {
      // left flux minus right flux
      for (int k = 0; k < numberOfVariables; k++) {
        lduh[i * basisSize * numberOfVariables + j * numberOfVariables + k] -=
            (lFbnd[1 * basisSize * numberOfVariables + i * numberOfVariables +
                   k] *
                 kernels::FRCoeff[order][j] +
             lFbnd[0 * basisSize * numberOfVariables + i * numberOfVariables +
                   k] *
                 kernels::FLCoeff[order][j]) *
            updateSize;
      }
    }
  }

  // y faces
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      const double weight = kernels::gaussLegendreWeights[order][i] *
                            kernels::gaussLegendreWeights[order][j];
      const double updateSize = weight / dx[1];

      // back flux minus front flux
      for (int k = 0; k < numberOfVariables; k++) {
        lduh[i * basisSize * numberOfVariables + j * numberOfVariables + k] -=
            (lFbnd[3 * basisSize * numberOfVariables + j * numberOfVariables +
                   k] *
                 kernels::FRCoeff[order][i] +
             lFbnd[2 * basisSize * numberOfVariables + j * numberOfVariables +
                   k] *
                 kernels::FLCoeff[order][i]) *
            updateSize;
      }
    }
  }
}

#endif  // DIMENSIONS == 2
