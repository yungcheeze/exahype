#include "EulerFlow3d/dg/ADERDG.h"

#include "EulerFlow3d/quad/GaussLegendre.h"

#include "EulerFlow3d/dg/DGMatrices.h"

// 3D
void exahype::dg::surfaceIntegral(
    double * lduh,
    const double * const dx,
    const int nvar,
    const int basisSize,
    const double * const FLeft,
    const double * const FRight,
    const double * const FFront,
    const double * const FBack,
    const double * const FBottom,
    const double * const FTop
) {
  // todo insert your code here
}

// 2D
void exahype::dg::surfaceIntegral(
    double * lduh,
    const double * const dx,
    const int nvar,
    const int basisSize,
    const double * const FLeft,
    const double * const FRight,
    const double * const FFront,
    const double * const FBack
) {
  // x direction (independent from the y and z)
  for (int jj=0; jj<basisSize; jj++) {
    const int nodeIndex     = jj;
    const int dofStartIndex = nodeIndex * nvar;

    double weight =  quad::gaussLegendreWeights[basisSize-1][jj];

    for(int mm=0; mm < basisSize; mm++) {
      const int mmNodeIndex         = mm + basisSize * jj;
      const int mmDofStartIndex     = mmNodeIndex * nvar;

      for(int ivar=0; ivar < nvar; ivar++) {
        lduh[mmDofStartIndex+ivar]
            -=  weight/dx[0] * ( dg::FRCoeff[mm] * FRight[dofStartIndex+ivar] - dg::FLCoeff[mm] * FLeft[dofStartIndex+ivar] ); // todo FL/RCoeff is hard coded
      }
    }
  }
  // Above seems okay!

  // y direction (independent from the y and z)
  for (int jj=0; jj<basisSize; jj++) {
    const int nodeIndex     = jj;
    const int dofStartIndex = nodeIndex * nvar;

    double weight =  quad::gaussLegendreWeights[basisSize-1][jj];

    for(int mm=0; mm < basisSize; mm++) {
      const int mmNodeIndex         = jj + basisSize * mm;
      const int mmDofStartIndex     = mmNodeIndex * nvar;

      for(int ivar=0; ivar < nvar; ivar++) {
        lduh[mmDofStartIndex+ivar]
           -=  weight/dx[1] * ( dg::FRCoeff[mm] * FBack[dofStartIndex+ivar] - dg::FLCoeff[mm] * FFront[dofStartIndex+ivar] ); // todo FL/RCoeff is hard coded
      }
    }
  }
  // Above seems okay!
}
