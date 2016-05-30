#include "../../Kernels.h"

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

#if DIMENSIONS == 2

// WARNING: Untested!

void volumeIntegralLinear(double* lduh, const double* const lFhi,
                          const tarch::la::Vector<DIMENSIONS, double>& dx,
                          const int numberOfVariables, const int basisSize) {
  // for linear non-conservative PDE, the volume integral is trivial, since it
  // only involves the element mass matrix, which later will cancel

  const int basisSize2 = basisSize * basisSize;
  const int order = basisSize - 1;

  std::memset(lduh, 0, basisSize2 * numberOfVariables * sizeof(double));

  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      double weight = kernels::gaussLegendreWeights[order][i] *
                      kernels::gaussLegendreWeights[order][i] *
                      kernels::gaussLegendreWeights[order][j];

      // Fortran: lduh(:,k,j) = -SUM(lFhi(:,k,j,1:nDim), dim = 4) * weight
      for (int k = 0; k < numberOfVariables; k++) {
        for (int l = 0; l < 2; l++) {
          lduh[i * basisSize * numberOfVariables + j * numberOfVariables + k] -=
              weight * lFhi[l * basisSize2 * numberOfVariables +

                            i * basisSize * numberOfVariables +
                            j * numberOfVariables + k];
        }
      }
    }
  }
}

#endif

}  // namespace c
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels
