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
 
#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

using std::endl;
using std::cout;

#if DIMENSIONS == 2

extern "C" {
void elementupdate_(double *luh, double *lduh, double *dt);
}
void kernels::aderdg::generic::c::solutionUpdate(double *luh,
                                                 const double *const lduh,
                                                 const double dt,
                                                 const int numberOfVariables,
                                                 const int basisSize) {
  const int order = basisSize - 1;

  for (int ii = 0; ii < basisSize; ii++) {
    for (int jj = 0; jj < basisSize; jj++) {
      const int nodeIndex = jj + basisSize * ii;
      const int dofStartIndex = nodeIndex * numberOfVariables;

      const double weight = kernels::gaussLegendreWeights[order][ii] *
                            kernels::gaussLegendreWeights[order][jj];
      const double updateSize = dt / weight;

      for (int ivar = 0; ivar < numberOfVariables; ivar++) {
        luh[dofStartIndex + ivar] += lduh[dofStartIndex + ivar] * updateSize;
      }
    }
  }
}

#endif  // DIMENSIONS == 2
