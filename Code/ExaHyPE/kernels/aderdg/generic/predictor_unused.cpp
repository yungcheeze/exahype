#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"


void kernels::aderdg::generic::predictor(
    double* lQhi,
    double* lFhi,
    const double* const lQi,
    const double* const lFi,
    const double predictorTimeStepSize,
    const int numberOfVariables,
    const int basisSize
) {
  /////////////////////////////////////////////////
  // Post processing of the predictor:
  // Immediately compute the time-averaged space-time polynomials
  /////////////////////////////////////////////////
  constexpr int dim         = 2;         // 2
  //constexpr int dimTimesTwo = (2*2);     // 4
  constexpr int BASISSIZE   = (3+1); // basisSize=order+1
  constexpr int NVAR        = 5;     // numberOfVariables
  constexpr int numberOfDof           = NVAR * BASISSIZE * BASISSIZE;/*tarch::la::aPowI(dim,BASISSIZE);*/
  constexpr int numberOfFluxDof  = numberOfDof * dim;

  memset((double *) lQhi,0,sizeof(double) * numberOfDof);
  memset((double *) lFhi,0,sizeof(double) * numberOfFluxDof);

  // memory layout of lFhi:
  // lFhi = [ lFhi_x | lFhi_y ] ordered as
  // (a) lFhi_x[nDOF_y][nDOF_x][nVar]
  // (b) lFhi_y[nDOF_x][nDOF_y][nVar]
  // Note the order of lFhi_y. Rationale is that matrix multiplications
  // then no longer have strided access pattern
  //
  // For 3D the variables should be as follows.
  // lFhi = [ lFhi_x | lFhi_y | lFhi_z ] where
  // (a) lFhi_x[nDOF_z][nDOF_y][nDOF_x][nVar]
  // (b) lFhi_y[nDOF_z][nDOF_x][nDOF_y][nVar]
  // (c) lFhi_z[nDOF_x][nDOF_y][nDOF_z][nVar]
  //
  // Above seems to work!
  double * lFhi_x = &lFhi[0];
  double * lFhi_y = &lFhi[numberOfDof];

  const double*    Q;
  const double*    f;
  const double*    g;

  for (int ii=0; ii<BASISSIZE; ii++) { // loop over dof
    for (int jj=0; jj<BASISSIZE; jj++) {
      const int nodeIndex         = jj + BASISSIZE * ii;
      const int dofStartIndex     = nodeIndex * NVAR;
      //const int fluxDofStartIndex = dim * dofStartIndex;

      for (int ll=0; ll<BASISSIZE; ll++) { // loop over dof
        const int spaceTimeNodeIndex         = nodeIndex  + BASISSIZE * BASISSIZE * ll;
        const int spaceTimeDofStartIndex     = spaceTimeNodeIndex * NVAR;
        const int spaceTimeFluxDofStartIndex = spaceTimeDofStartIndex * dim;

        Q = &lQi[spaceTimeDofStartIndex];

        f = &lFi[spaceTimeFluxDofStartIndex     ];
        g = &lFi[spaceTimeFluxDofStartIndex+NVAR];

        double weight = kernels::gaussLegendreWeights[BASISSIZE-1][ll];

        //double * temp = &(lQhi[dofStartIndex]);
        for(int ivar=0; ivar < NVAR; ivar++) {
          lQhi[dofStartIndex+ivar] += weight * Q[ivar];

          lFhi_x[dofStartIndex+ivar] += weight * f[ivar];   // lFhi_x[nDOF_y][nDOF_x][nVar]
          //lFhi_y[dofStartIndex+ivar] += weight * g[ivar]; // lFhi_y[nDOF_y][nDOF_x][nVar] (former version)
          lFhi_y[(ii + BASISSIZE * jj) * NVAR+ivar] += weight * g[ivar];  // lFhi_y(DOFx,DOFy,nVar)
        }
      }
    }
  }
}


