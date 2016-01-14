#include "EulerFlow/dg/ADERDG.h"

#include "string.h"

#include "tarch/la/ScalarOperations.h"

#include "EulerFlow/quad/GaussLegendre.h"

#include "EulerFlow/dg/DGMatrices.h"

// 3D specialisation
template <>
void exahype::dg::volumeIntegral<3>(
    double* /*out*/ lduh,
    const double * const lFhi,
    const double * const dx
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
    const double * const lFhi,
    const double * const dx
) {
  constexpr int dim         = DIMENSIONS;                 // 2
  constexpr int nvar        = EXAHYPE_NVARS;
  constexpr int basisSize   = EXAHYPE_ORDER+1;
  constexpr int numberOfDof = nvar * power(basisSize,dim);

  const double * f;
  const double * g;
  const double * lFhi_x;
  const double * lFhi_y;

  // lFhi's memory layout:
  // lFhi = [ lFhi_x | lFhi_y ] ordered as follows
  // (a) lFhi_x[nDOF_y][nDOF_x][nVar]
  // (b) lFhi_y[nDOF_y][nDOF_x][nVar]

  memset(lduh,0,sizeof(double) * numberOfDof);
  lFhi_x = &lFhi[0];
  lFhi_y = &lFhi[numberOfDof];

  // access lduh(nDOF[2] x nDOF[1] x nvar) in the usual 3D array manner
  typedef double tensor_t[basisSize][nvar];
  tensor_t *lduh3D = (tensor_t *)lduh;

  // Compute the "derivatives" (contributions of the stiffness matrix)
  // x direction (independent from the y and z derivatives)

  for (int ii=0; ii<basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {

      double weight = exahype::quad::gaussLegendreWeights[ii];

      // MATMUL: lFhi_x * Kxi
      for(int mm=0; mm < basisSize; mm++) {
        const int mmNodeIndex         = mm + basisSize * ii;
        const int mmDofStartIndex     = mmNodeIndex * nvar;

        //f = &lFhi[mmFluxDofStartIndex]; // TODO
        //f = &lFhi[mmDofStartIndex];

        for(int ivar=0; ivar < nvar; ivar++) {
          //lduh3D[ii][jj][ivar] += weight/dx[0] * dg::Kxi[jj][mm] * f[ivar];
          lduh3D[ii][jj][ivar] += weight/dx[0] * dg::Kxi[jj][mm] * lFhi_x[mmDofStartIndex+ivar];
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

      // MATMUL: lFhi_y * Kxi
      for(int mm=0; mm < basisSize; mm++) {
        const int mmNodeIndex         = jj + basisSize * mm;
        const int mmDofStartIndex     = mmNodeIndex * nvar;

        //g = &lFhi[mmFluxDofStartIndex+nvar];  // TODO
        //g = &lFhi[numberOfDof+mmDofStartIndex];

        for(int ivar=0; ivar < nvar; ivar++) {
          //lduh3D[ii][jj][ivar] += weight/dx[1] * dg::Kxi[ii][mm] * g[ivar]; // lFhi_y * Kxi
          lduh3D[ii][jj][ivar] += weight/dx[1] * dg::Kxi[ii][mm] * lFhi_y[mmDofStartIndex+ivar];
        }
      }
    }
  }
  // Above seems okay!
}
