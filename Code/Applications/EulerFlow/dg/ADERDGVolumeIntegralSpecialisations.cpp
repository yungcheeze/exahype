#include "EulerFlow/dg/ADERDG.h"

#include "string.h"

#include "tarch/la/ScalarOperations.h"

#include "EulerFlow/quad/GaussLegendre.h"

#include "EulerFlow/dg/DGMatrices.h"

// 3D specialisation
template <>
void exahype::dg::volumeIntegral<3>(
    double* lduh,
    const double * const lFhi,
    const double * const dx
) {
  constexpr int dim         = DIMENSIONS;     // 3
  constexpr int dimTimesTwo = (2*DIMENSIONS); // 6
  constexpr int nvar = EXAHYPE_NVARS;
  constexpr int basisSize = EXAHYPE_ORDER+1;

  // todo insert your code here
  // todo insert your code here
}

// 2D specialisation
template <>
void exahype::dg::volumeIntegral<2>(
    double* lduh,
    const double * const lFhi,
    const double * const dx
) {
  constexpr int dim = DIMENSIONS; // 2
  constexpr int nvar = EXAHYPE_NVARS;
  constexpr int basisSize = EXAHYPE_ORDER+1;
  constexpr int numberOfDof = nvar * power(basisSize,dim);

  const double * f;
  const double * g;

  memset(lduh,0,sizeof(double) * numberOfDof);

  // Compute the "derivatives" (contributions of the stiffness matrix)
  // x direction (independent from the y and z derivatives)
  for (int ii=0; ii<basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {
      const int nodeIndex     = ii + basisSize * jj;
      const int dofStartIndex = nodeIndex * nvar;

      double weight = exahype::quad::gaussLegendreWeights[jj];

      for(int mm=0; mm < basisSize; mm++) {
        const int mmNodeIndex         = mm + basisSize * jj;
        const int mmDofStartIndex     = mmNodeIndex * nvar;
        const int mmFluxDofStartIndex = mmDofStartIndex * dim;

        f = &lFhi[mmFluxDofStartIndex];

        for(int ivar=0; ivar < nvar; ivar++) {
          lduh[dofStartIndex+ivar] +=  weight/dx[0] * dg::Kxi[ii][mm] * f[ivar];
        }
      }
    }
  }
  // Above seems okay!

  // Compute the "derivatives" (contributions of the stiffness matrix)
  /// y direction (independent from the y and z derivatives)
  for (int ii=0; ii<basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {
      const int nodeIndex     = ii + basisSize * jj;
      const int dofStartIndex = nodeIndex * nvar;

      double weight = exahype::quad::gaussLegendreWeights[ii];

      for(int mm=0; mm < basisSize; mm++) {
        const int mmNodeIndex         = ii + basisSize * mm;
        const int mmDofStartIndex     = mmNodeIndex * nvar;
        const int mmFluxDofStartIndex = mmDofStartIndex * dim;

        g = &lFhi[mmFluxDofStartIndex+nvar];

        double * du = &lduh[0];
        for(int ivar=0; ivar < nvar; ivar++) {
          du[dofStartIndex+ivar] +=  weight/dx[1] * dg::Kxi[jj][mm] * g[ivar];
        }
      }
    }
  }
  // Above seems okay!
}
