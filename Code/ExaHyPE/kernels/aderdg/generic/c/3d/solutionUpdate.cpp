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

#include "../../../../GaussLegendreQuadrature.h"

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

#if DIMENSIONS == 3

void solutionUpdate(double* luh, const double* const lduh, const double dt,
                    const int numberOfVariables, const int basisSize) {
  const int order = basisSize - 1;

  for (int i = 0; i < basisSize; i++) {
    const int i_offset = i * basisSize * basisSize * numberOfVariables;

    for (int j = 0; j < basisSize; j++) {
      const int j_offset = j * basisSize * numberOfVariables;

      for (int k = 0; k < basisSize; k++) {
        const int k_offset = k * numberOfVariables;

        const double weight = kernels::gaussLegendreWeights[order][i] *
                              kernels::gaussLegendreWeights[order][j] *
                              kernels::gaussLegendreWeights[order][k];
        const double updateSize = dt / weight;

        for (int l = 0; l < numberOfVariables; l++) {
          const int offset = i_offset + j_offset + k_offset + l;
          luh[offset] += lduh[offset] * updateSize;
        }
      }
    }
  }
}

#endif  // DIMENSIONS == 3

}  // namespace c
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels
