/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#include "../../Kernels.h"

#if DIMENSIONS == 2

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

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
      const double weight = kernels::gaussLegendreWeights[order][j];
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

}  // namespace c
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels

#endif  // DIMENSIONS == 2
