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
 
#include <tarch/la/Vector.h>

#include "kernels/aderdg/generic/Kernels.h"

#if DIMENSIONS == 3

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

void surfaceIntegralLinear(double *lduh, const double *const lFbnd,
                           const tarch::la::Vector<DIMENSIONS, double> &dx,
                           const int numberOfVariables, const int basisSize) {
  const int basisSize2 = basisSize * basisSize;
  const int order = basisSize - 1;

  // x faces
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      const double weight = kernels::gaussLegendreWeights[order][i] *
                            kernels::gaussLegendreWeights[order][j];
      const double updateSize = weight / dx[0];

      for (int k = 0; k < basisSize; k++) {
        // left flux minus right flux
        for (int l = 0; l < numberOfVariables; l++) {
          lduh[i * basisSize2 * numberOfVariables +
               j * basisSize * numberOfVariables + k * numberOfVariables + l] -=
              (lFbnd[1 * basisSize2 * numberOfVariables +
                     i * basisSize * numberOfVariables + j * numberOfVariables +
                     l] *
                   kernels::FRCoeff[order][k] +
               lFbnd[0 * basisSize2 * numberOfVariables +
                     i * basisSize * numberOfVariables + j * numberOfVariables +
                     l] *
                   kernels::FLCoeff[order][k]) *
              updateSize;
        }
      }
    }
  }

  // y faces
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      for (int k = 0; k < basisSize; k++) {
        const double weight = kernels::gaussLegendreWeights[order][i] *
                              kernels::gaussLegendreWeights[order][k];
        const double updateSize = weight / dx[1];

        // back flux minus front flux
        for (int l = 0; l < numberOfVariables; l++) {
          lduh[i * basisSize2 * numberOfVariables +
               j * basisSize * numberOfVariables + k * numberOfVariables + l] -=
              (lFbnd[3 * basisSize2 * numberOfVariables +
                     i * basisSize * numberOfVariables + k * numberOfVariables +
                     l] *
                   kernels::FRCoeff[order][j] +
               lFbnd[2 * basisSize2 * numberOfVariables +
                     i * basisSize * numberOfVariables + k * numberOfVariables +
                     l] *
                   kernels::FLCoeff[order][j]) *
              updateSize;
        }
      }
    }
  }

  // z faces
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      for (int k = 0; k < basisSize; k++) {
        const double weight = kernels::gaussLegendreWeights[order][j] *
                              kernels::gaussLegendreWeights[order][k];
        const double updateSize = weight / dx[2];

        // bottom flux minus top flux
        for (int l = 0; l < numberOfVariables; l++) {
          lduh[i * basisSize2 * numberOfVariables +
               j * basisSize * numberOfVariables + k * numberOfVariables + l] -=
              (lFbnd[5 * basisSize2 * numberOfVariables +
                     j * basisSize * numberOfVariables + k * numberOfVariables +
                     l] *
                   kernels::FRCoeff[order][i] +
               lFbnd[4 * basisSize2 * numberOfVariables +
                     j * basisSize * numberOfVariables + k * numberOfVariables +
                     l] *
                   kernels::FLCoeff[order][i]) *
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
