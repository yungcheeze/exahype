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

#ifndef _EXAHYPE_KERNELS_LIMITER_GENERIC_H_
#define _EXAHYPE_KERNELS_ADERDG_GENERIC_H_

#include <algorithm>
#include <stdexcept>
#include <stdlib.h>

#include "../../LimiterProjectionMatrices.h"
#include "../../GaussLegendreQuadrature.h"
#include "../../KernelUtils.h"

#include "peano/utils/Globals.h"

namespace kernels {
namespace limiter {
namespace generic {
namespace c {

int getLimBasisSize(const int basisSize);

double* getGaussLobattoData(const double* const luh, const int numberOfVariables, const int basisSize); //for test only

template <bool anticipateNewDGSolution>
void compareWithADERDGSolutionAtGaussLobattoNodes(const double* const luh, const int numberOfVariables, const int basisSize, double* const min, double* const max, const double* const lduh, const double dt);

void projectOnFVLimiterSpace(const double* const luh, const int numberOfVariables, const int basisSize, const int basisSizeLim, double* const lim);

void findCellLocalLimMinAndMax(const double* const lim, const int numberOfVariables, const int basisSize, double* const localMin, double* const localMax); //deprecated
void findCellLocalMinAndMax(const double* const luh, const int numberOfVariables, const int basisSize, double* const localMin, double* const localMax);
bool isTroubledCell(const double* const luh, const int numberOfVariables, const int basisSize, const double* const troubledMin, const double* const troubledMax);

void projectOnADERDGSpace(const double* const lim, const int numberOfVariables, const int basisSizeLim, const int basisSize, double* const luh);

inline double anticipateLuh(const double* const luh, const double* const lduh, const double dt, const int order, const int idx, const int x, const int y, const int z) {
  return (luh[idx] + kernels::gaussLegendreWeights[order][x] * kernels::gaussLegendreWeights[order][y]
#if DIMENSIONS == 3
                                                  * kernels::gaussLegendreWeights[order][z]
#endif
                                                  /dt * lduh[idx]);
}

} // namespace c
} // namespace generic
} // namespace limiter
} // namespace kernel


#endif //_EXAHYPE_KERNELS_LIMITER_GENERIC_H_
