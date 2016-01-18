#include "EulerFlow/dg/ADERDG.h"

#include "string.h"

#include "tarch/la/ScalarOperations.h"

#include "EulerFlow/quad/GaussLegendre.h"

#include "EulerFlow/dg/DGMatrices.h"

// 3D specialisation
template <>
void exahype::dg::volumeIntegral<3>(
    double* /*out*/ lduh,
    const double * const /*in*/ lFhi,
    const double * const /*in*/ dx
) {
  constexpr int dim         = DIMENSIONS;     // 3
  constexpr int dimTimesTwo = (2*DIMENSIONS); // 6
  constexpr int nvar        = EXAHYPE_NVARS;
  constexpr int basisSize   = EXAHYPE_ORDER+1;

  // todo insert your code here
}

// 2D specialisation
template <>
void exahype::dg::volumeIntegral<2>(
    double* /*out*/ lduh,
    const double * const /*in*/ lFhi,
    const double * const /*in*/ dx
) {
  constexpr int dim         = DIMENSIONS;                 // 2
  constexpr int nvar        = EXAHYPE_NVARS;
  constexpr int basisSize   = EXAHYPE_ORDER+1;
  constexpr int numberOfDof = nvar * power(basisSize,dim);

  // memory layout of lFhi:
  // lFhi = [ lFhi_x | lFhi_y ] ordered as follows
  // (a) lFhi_x[nDOF_y][nDOF_x][nVar]
  // (b) lFhi_y[nDOF_y][nDOF_x][nVar]
  // let's not bother with offsets and define separate flux matrices
  const double * lFhi_x = &lFhi[0];            // f flux
  const double * lFhi_y = &lFhi[numberOfDof];  // g flux

  memset(lduh,0,sizeof(double) * numberOfDof);

  // access lduh(nDOF[2] x nDOF[1] x nvar) in the usual 3D array manner
  typedef double tensor_t[basisSize][nvar];
  tensor_t *lduh3D = (tensor_t *)lduh;

  // Compute the "derivatives" (contributions of the stiffness matrix)
  // x direction (independent from the y and z derivatives)
  for (int ii=0; ii<basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {

      double weight = exahype::quad::gaussLegendreWeights[ii];

      // MATMUL: Kxi * lFhi_x
      for(int mm=0; mm < basisSize; mm++) {
        const int mmNodeIndex         = mm + basisSize * ii;
        const int mmDofStartIndex     = mmNodeIndex * nvar;

        for(int ivar=0; ivar < nvar; ivar++) {
          lduh3D[ii][jj][ivar] += weight/dx[0] * dg::Kxi[jj][mm] * lFhi_x[mmDofStartIndex+ivar];
        }
      }
    }
  }


  // Compute the "derivatives" (contributions of the stiffness matrix)
  // y direction (independent from the y and z derivatives)
  for (int ii=0; ii<basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {

      double weight = exahype::quad::gaussLegendreWeights[jj];

      // MATMUL: Kxi * lFhi_y
      for(int mm=0; mm < basisSize; mm++) {
        const int mmNodeIndex         = jj + basisSize * mm;
        const int mmDofStartIndex     = mmNodeIndex * nvar;

        for(int ivar=0; ivar < nvar; ivar++) {
          lduh3D[ii][jj][ivar] += weight/dx[1] * dg::Kxi[ii][mm] * lFhi_y[mmDofStartIndex+ivar];
        }
      }
    }
  }

}
