/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released unter the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#include "../../Kernels.h"

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

#if DIMENSIONS == 2

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
                      kernels::gaussLegendreWeights[order][j];

      // Fortran: lduh(:,j,i) = -SUM(lFhi(:,j,i,1:nDim), dim = 4) * weight
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
