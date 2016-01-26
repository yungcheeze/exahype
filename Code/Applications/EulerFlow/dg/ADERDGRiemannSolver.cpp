#include "EulerFlow/dg/ADERDG.h"

#include "string.h"

#include "tarch/la/ScalarOperations.h"

#include "EulerFlow/quad/GaussLegendre.h"

#include "EulerFlow/dg/DGMatrices.h"

#include "EulerFlow/problem/Problem.h"

// explicit specialisation
template <>
void exahype::dg::solveRiemannProblem<3>(
            double * restrict /*out*/ FL,
            double * restrict /*out*/ FR,
            const double * restrict const /*in*/ QL,
            const double * restrict const /*in*/ QR,
            double * restrict /*local*/ QavL,
            double * restrict /*local*/ QavR,
            double * restrict /*local*/ lambdaL,
            double * restrict /*local*/ lambdaR,
            const double /*unused*/ dt,
            const double /*unused*/ hFace,
            const double * restrict const /*in*/ n
) {
  constexpr int dim       = DIMENSIONS;       // 3
  constexpr int nvar      = EXAHYPE_NVARS;
  constexpr int basisSize = EXAHYPE_ORDER+1;

  // todo insert your code here
}

template <>
void exahype::dg::solveRiemannProblem<2>(
            double * restrict /*out*/ FL,
            double * restrict /*out*/ FR,
            const double * restrict const /*in*/ QL,
            const double * restrict const /*in*/ QR,
            double * restrict /*local*/ QavL,
            double * restrict /*local*/ QavR,
            double * restrict /*local*/ lambdaL,
            double * restrict /*local*/ lambdaR,
            const double /*unused*/ dt,
            const double /*unused*/ hFace,
            const double * restrict const /*in*/ n
) {
  constexpr int dim       = DIMENSIONS;       // 2
  constexpr int nvar      = EXAHYPE_NVARS;
  constexpr int basisSize = EXAHYPE_ORDER+1;

  __assume_aligned(QavL, ALIGNMENT);
  __assume_aligned(QavR, ALIGNMENT);
  __assume_aligned(lambdaL, ALIGNMENT);
  __assume_aligned(lambdaR, ALIGNMENT);

  // Compute the average states from the left and the right, which we need to compute the numerical dissipation

  memset((double *) QavL,0,nvar * sizeof(double));
  memset((double *) QavR,0,nvar * sizeof(double));

  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    const int nodeIndex     = ii;
    const int dofStartIndex = nodeIndex * nvar;

    double weight =  quad::gaussLegendreWeights[ii];

    for(int ivar=0; ivar < nvar; ivar++) {
      QavL[ivar] +=  weight * QL[dofStartIndex+ivar];
      QavR[ivar] +=  weight * QR[dofStartIndex+ivar];
    }
  }
  //
  // Here, we implement a very simple Rusanov scheme with scalar dissipation (smax*Id).
  // We can change this into a more sophisticated Osher or HLLEM Riemann solver whenever needed.
  //
  exahype::problem::PDEEigenvalues(QavL,n,lambdaL);
  exahype::problem::PDEEigenvalues(QavR,n,lambdaR);

  double sMax = 0;
  for(int ivar=0; ivar < nvar; ivar++) {
    sMax = std::max(sMax,std::max(fabs(lambdaL[ivar]),fabs(lambdaR[ivar])));
  }
  //
  // We now compute the numerical flux. Note that the scheme is at the moment written in
  // CONSERVATION FORM => no fluctuations, but real fluxes.
  // Later, this will be converted into the left and right fluctuations.
  //
  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    const int nodeIndex     = ii;
    const int dofStartIndex = nodeIndex * nvar;

    for(int ivar=0; ivar < nvar; ivar++) {
      FL[dofStartIndex+ivar] = 0.5 * (FL[dofStartIndex+ivar] + FR[dofStartIndex+ivar])
                              -0.5 * sMax *  (QR[dofStartIndex+ivar] - QL[dofStartIndex+ivar]);

      FR[dofStartIndex+ivar] = FL[dofStartIndex+ivar];
    }
  }
}
