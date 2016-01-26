/*
#include "exahype/aderdg/ADERDG.h"

#include "string.h"

#include "tarch/la/ScalarOperations.h"

#include "kernels/quad/GaussLegendre.h"
#include "kernels/aderdg/DGMatrices.h"
#include "Kernels.h"


// explicit specialisation
template <>
void exahype::aderdg::riemannSolver<3>(
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
            const double * const n
) {
  constexpr int dim       = DIMENSIONS;       // 3
  constexpr int nvar      = EXAHYPE_NVARS;
  constexpr int basisSize = EXAHYPE_ORDER+1;

  // todo insert your code here
}

template <>
void exahype::aderdg::riemannSolver<2>(
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
            const double * const n
) {
  constexpr int DIM2       = 2;             // 2
  constexpr int nvar      = EXAHYPE_NVARS;
  constexpr int basisSize = EXAHYPE_ORDER+1;

  // Compute the average states from the left and the right, which we need to compute the numerical dissipation
  double sMax = 0;
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
  exahype::pde::PDEEigenvalues2d(QavL,nvar,n,DIM2,lambdaL);
  exahype::pde::PDEEigenvalues2d(QavR,nvar,n,DIM2,lambdaR);

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
*/
