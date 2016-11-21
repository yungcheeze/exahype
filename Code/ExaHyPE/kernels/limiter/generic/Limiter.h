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

/**
 * \brief Projection ADERDG -> FV
 *
 * Projects the ADERDG solution onto
 * the finite volumes limiter space.
 *
 * \param[in] basisSize The size of the ADER-DG basis per coordinate axis (order+1).
 * \param[in] ghostLayerWidth The ghost layer width in cells of the finite volumes patch
 */
void projectOnFVLimiterSpace(const double* const luh, const int numberOfVariables, const int basisSize, const int ghostLayerWidth, double* const lim);
/**
 * \brief Projection FV -> ADERDG
 *
 * Projects the finite volumes limiter solution onto
 * the DG space.
 *
 * \param[in] basisSize The size of the ADER-DG basis per coordinate axis (order+1)
 * \param[in] ghostLayerWidth The ghost layer width in cells of the finite volumes patch.
 */
void projectOnDGSpace(const double* const lim, const int numberOfVariables, const int basisSize, const int ghostLayerWidth, double* const luh);

// Get the local min/max from the DG and Gauss Lobatto nodes
void findCellLocalMinAndMax(const double* const luh, const int numberOfVariables, const int basisSize, double* const localMinPerVariables, double* const localMaxPerVariables);

//Test if the DG solution is troubled
bool isTroubledCell(const double* const luh, const int numberOfVariables, const int basisSize, const double* const boundaryMinPerVariables, const double* const boundaryMaxPerVariables);

// TODO(Dominic): @JM: We have to do a rollback in every neighbour cell of the troubled cells. Furthermore, the
// troubled cells are not that many compared to the non-troubled ones. Thus, I decided to get
// rid of the solution anticipation. I am sorry for the confusion.
//bool isTroubledCell(const double* const luh, const double* const lduh, const double dt, const int numberOfVariables, const int basisSize, const double* const boundaryMinPerVariables, const double* const boundaryMaxPerVariables);

//************************
//*** Helper functions ***
//************************

//inline double anticipateLuh(const double* const luh, const double* const lduh, const double dt, const int order, const int idx, const int x, const int y, const int z) {
//  double weight =
//  #if DIMENSIONS == 3
//  kernels::gaussLegendreWeights[order][z] *
//  #endif
//  kernels::gaussLegendreWeights[order][y] * kernels::gaussLegendreWeights[order][x];
//
//  return (luh[idx] + dt / weight * lduh[idx]); // TODO(Dominic): The compiler might not able to optimise for the dt=0 case.
//}

inline int getBasisSizeLim(const int basisSize) {
  return 2*(basisSize-1)+1;
}

//*************************
//*** Private functions ***
//*************************

//Projection ADERDG -> Gauss-Lobatto, for test only
double* getGaussLobattoData(const double* const luh, const int numberOfVariables, const int basisSize); 

void compareWithADERDGSolutionAtGaussLobattoNodes(const double* const luh, const int numberOfVariables, const int basisSize, double* const min, double* const max);

} // namespace c
} // namespace generic
} // namespace limiter
} // namespace kernel


#endif //_EXAHYPE_KERNELS_LIMITER_GENERIC_H_
