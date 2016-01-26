/*
#include "exahype/aderdg/ADERDG.h"

#include "kernels/quad/GaussLegendre.h"
#include "kernels/aderdg/DGMatrices.h"

// 3D
void exahype::aderdg::surfaceIntegral(
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
void exahype::aderdg::surfaceIntegral(
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

    double weight =  quad::gaussLegendreWeights[jj];

    for(int mm=0; mm < basisSize; mm++) {
      const int mmNodeIndex         = mm + basisSize * jj;
      const int mmDofStartIndex     = mmNodeIndex * nvar;

      for(int ivar=0; ivar < nvar; ivar++) {
        lduh[mmDofStartIndex+ivar]
            -=  weight/dx[0] * ( aderdg::FRCoeff[mm] * FRight[dofStartIndex+ivar] - aderdg::FLCoeff[mm] * FLeft[dofStartIndex+ivar] ); // todo FL/RCoeff is hard coded
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
      const int mmNodeIndex         = jj + basisSize * mm;
      const int mmDofStartIndex     = mmNodeIndex * nvar;

      for(int ivar=0; ivar < nvar; ivar++) {
        lduh[mmDofStartIndex+ivar]
           -=  weight/dx[1] * ( aderdg::FRCoeff[mm] * FBack[dofStartIndex+ivar] - aderdg::FLCoeff[mm] * FFront[dofStartIndex+ivar] ); // todo FL/RCoeff is hard coded
      }
    }
  }
  // Above seems okay!
}
*/
