#include "EulerFlow/dg/ADERDG.h"

#include "EulerFlow/quad/GaussLegendre.h"

#include "EulerFlow/dg/DGMatrices.h"

// 3D
void exahype::dg::surfaceIntegral(
    double * restrict /*inout*/ lduh,
    const double * restrict const /*in*/ dx,
    const double * restrict const /*in*/ FLeft,
    const double * restrict const /*in*/ FRight,
    const double * restrict const /*in*/ FFront,
    const double * restrict const /*in*/ FBack,
    const double * restrict const /*in*/ FBottom,
    const double * restrict const /*in*/ FTop
) {
  constexpr int nvar        = EXAHYPE_NVARS;
  constexpr int basisSize   = EXAHYPE_ORDER+1;

  // todo insert your code here
}

// 2D
void exahype::dg::surfaceIntegral(
    double * restrict /*inout*/ lduh,
    const double * restrict const /*in*/ dx,
    const double * restrict const /*in*/ FLeft,
    const double * restrict const /*in*/ FRight,
    const double * restrict const /*in*/ FFront,
    const double * restrict const /*in*/ FBack
) {
  constexpr int nvar        = EXAHYPE_NVARS;
  constexpr int basisSize   = EXAHYPE_ORDER+1;

  // access lduh(nDOF[2] x nDOF[1] x nvar) in the usual 3D array manner
  typedef double tensor_t[basisSize][nvar];
  tensor_t *lduh3D = (tensor_t *)lduh;

  // x direction (independent from the y and z)
  for (int jj=0; jj<basisSize; jj++) {
    const int nodeIndex     = jj;
    const int dofStartIndex = nodeIndex * nvar;

    double weight =  quad::gaussLegendreWeights[jj];

    for(int mm=0; mm < basisSize; mm++) {
      for(int ivar=0; ivar < nvar; ivar++) {
        lduh3D[jj][mm][ivar]
            -=  weight/dx[0] * ( dg::FRCoeff[mm] * FRight[dofStartIndex+ivar] - dg::FLCoeff[mm] * FLeft[dofStartIndex+ivar] );
      }
    }
  }

  // Above seems okay!

  // y direction (independent from the y and z)
  for (int jj=0; jj<basisSize; jj++) {
    const int nodeIndex     = jj;
    const int dofStartIndex = nodeIndex * nvar;

    double weight =  quad::gaussLegendreWeights[jj];

    for(int mm=0; mm < basisSize; mm++) {
      for(int ivar=0; ivar < nvar; ivar++) {
        lduh3D[mm][jj][ivar]
           -=  weight/dx[1] * ( dg::FRCoeff[mm] * FBack[dofStartIndex+ivar] - dg::FLCoeff[mm] * FFront[dofStartIndex+ivar] );
      }
    }
  }

  // Above seems okay!
}
