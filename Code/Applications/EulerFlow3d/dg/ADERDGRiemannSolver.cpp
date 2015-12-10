#include "EulerFlow3d/dg/ADERDG.h"

#include "string.h"

#include "tarch/la/ScalarOperations.h"

#include "EulerFlow3d/quad/GaussLegendre.h"

#include "EulerFlow3d/dg/DGMatrices.h"

#include "EulerFlow3d/problem/Problem.h"

// explicit specialisation
template <>
void exahype::dg::solveRiemannProblem<3>(
            double * FL,
            double * FR,
            const double * const QL,
            const double * const QR,
            double * QavL,
            double * QavR,
            double * lambdaL,
            double * lambdaR,
            const double dt,
            const double hFace,
            const double * const n,
            const int nvar,
            const int basisSize
) {
  constexpr int dim = 2;
  // todo insert your code here
}

template <>
void exahype::dg::solveRiemannProblem<2>(
            double * FL,
            double * FR,
            const double * const QL,
            const double * const QR,
            double * QavL,
            double * QavR,
            double * lambdaL,
            double * lambdaR,
            const double dt,
            const double hFace,
            const double * const n,
            const int nvar,
            const int basisSize
) {
  constexpr int dim = 2;

  // Compute the average states from the left and the right, which we need to compute the numerical dissipation
  double sMax = 0;
  memset((double *) QavL,0,nvar * sizeof(double));
  memset((double *) QavR,0,nvar * sizeof(double));

  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    const int nodeIndex     = ii;
    const int dofStartIndex = nodeIndex * nvar;

    double weight =  quad::gaussLegendreWeights[basisSize-1][ii];

    for(int ivar=0; ivar < nvar; ivar++) {
      QavL[ivar] +=  weight * QL[dofStartIndex+ivar];
      QavR[ivar] +=  weight * QR[dofStartIndex+ivar];
    }
  }
  //
  // Here, we implement a very simple Rusanov scheme with scalar dissipation (smax*Id).
  // We can change this into a more sophisticated Osher or HLLEM Riemann solver whenever needed.
  //
  exahype::problem::PDEEigenvalues(QavL,nvar,n,dim,lambdaL);
  exahype::problem::PDEEigenvalues(QavR,nvar,n,dim,lambdaR);

  sMax = 0;
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
