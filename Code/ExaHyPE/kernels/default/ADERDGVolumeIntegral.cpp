#include "kernels/compressibleeuler/ADERDG.h"

#include "string.h"

#include "tarch/la/ScalarOperations.h"

#include "EulerFlow/quad/GaussLegendre.h"

#include "kernels/compressibleeuler/DGMatrices.h"

// 3D specialisation
template <>
void exahype::dg::volumeIntegral<3>(
    double* /*out*/ lduh,
    const double * const lFhi,
    const double * const dx
) {
  constexpr int dim         = DIMENSIONS;     // 3
  constexpr int dimTimesTwo = (2*DIMENSIONS); // 6
  constexpr int nvar = EXAHYPE_NVARS;
  constexpr int basisSize = EXAHYPE_ORDER+1;

  // todo insert your code here
}

// 2D specialisation
template <>
void exahype::dg::volumeIntegral<2>(
    double* /*out*/ lduh,
    const double * const lFhi,
    const double * const dx
) {
  constexpr int dim = DIMENSIONS;             // 2
  constexpr int nvar = EXAHYPE_NVARS;
  constexpr int basisSize = EXAHYPE_ORDER+1;
  constexpr int numberOfDof = nvar * power(basisSize,dim);

  const double * f;
  const double * g;

  memset(lduh,0,sizeof(double) * numberOfDof);

  // access lduh(nDOF[2] x nDOF[1] x nvar) in the usual 3D array manner
  typedef double tensor_t[basisSize][nvar];
  tensor_t *lduh3D = (tensor_t *)lduh;

  // Compute the "derivatives" (contributions of the stiffness matrix)
  // x direction (independent from the y and z derivatives)
  for (int ii=0; ii<basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {

      double weight = exahype::quad::gaussLegendreWeights[ii];

      for(int mm=0; mm < basisSize; mm++) {
        const int mmNodeIndex         = mm + basisSize * ii;
        const int mmDofStartIndex     = mmNodeIndex * nvar;
        const int mmFluxDofStartIndex = mmDofStartIndex * dim;

        f = &lFhi[mmFluxDofStartIndex];

        for(int ivar=0; ivar < nvar; ivar++) {
          lduh3D[ii][jj][ivar] += weight/dx[0] * dg::Kxi[jj][mm] * f[ivar];
        }
      }
    }
  }

  // Above seems okay!

  // Compute the "derivatives" (contributions of the stiffness matrix)
  // y direction (independent from the y and z derivatives)
  for (int ii=0; ii<basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {

      double weight = exahype::quad::gaussLegendreWeights[jj];

      for(int mm=0; mm < basisSize; mm++) {
        const int mmNodeIndex         = jj + basisSize * mm;
        const int mmDofStartIndex     = mmNodeIndex * nvar;
        const int mmFluxDofStartIndex = mmDofStartIndex * dim;

        g = &lFhi[mmFluxDofStartIndex+nvar];

        for(int ivar=0; ivar < nvar; ivar++) {
          lduh3D[ii][jj][ivar] += weight/dx[1] * dg::Kxi[ii][mm] * g[ivar];
        }
      }
    }
  }
  // Above seems okay!
}
