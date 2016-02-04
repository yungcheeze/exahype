#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"


void kernels::aderdg::generic::solutionUpdate(double * luh,
                                              const double * const lduh,
                                              const tarch::la::Vector<DIMENSIONS,double>&  dx,
                                              const double dt,
                                              const int numberOfVariables,
                                              const int basisSize) {
  const int order = basisSize-1;

  for (int ii=0; ii<basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {
      const int nodeIndex     = jj + basisSize * ii;
      const int dofStartIndex = nodeIndex * numberOfVariables;

      const double weight     =  kernels::gaussLegendreWeights[order][ii] * kernels::gaussLegendreWeights[order][jj];
      const double updateSize = dt/weight;

      for(int ivar=0; ivar < numberOfVariables; ivar++) {
        luh[dofStartIndex+ivar] +=  lduh[dofStartIndex+ivar]*updateSize;
      }
    }
  }
}

void kernels::aderdg::generic::volumeIntegral(
    double * lduh,
    const double * const lFhi,
    const tarch::la::Vector<DIMENSIONS,double>& dx,
    const int numberOfVariables,
    const int basisSize
) {
  // todo Angelika
  // Please remove the typedefs in generic kernels again since numberOf(...)Dof is not
  // a compile time variable anymore
  constexpr int numberOfDof = 5 * (3+1)*(3+1); //tarch::la::aPowI(3+1,DIMENSIONS);
  constexpr int order       = 3;

  // memory layout of lFhi:
  // lFhi = [ lFhi_x | lFhi_y ] ordered as follows
  // (a) lFhi_x[nDOF_y][nDOF_x][nVar]
  // (b) lFhi_y[nDOF_x][nDOF_y][nVar]
  // Note the variable order of lFhi_y
  // let's not bother with offsets and define separate flux matrices
  const double * lFhi_x = &lFhi[0];            // f flux
  const double * lFhi_y = &lFhi[numberOfDof];  // g flux

  memset(lduh,0,sizeof(double) * numberOfDof);

  // access lduh(nDOF[2] x nDOF[1] x numberOfVariables) in the usual 3D array manner
  // @todo Angelika
  typedef double tensor_t[3+1][5];
  tensor_t *lduh3D = (tensor_t *)lduh;

  // Compute the "derivatives" (contributions of the stiffness matrix)
  // x direction (independent from the y and z derivatives)
  for (int ii=0; ii<basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {

      double weight = kernels::gaussLegendreWeights[order][ii];

      // MATMUL: Kxi * lFhi_x
      for(int mm=0; mm < basisSize; mm++) {
        const int mmNodeIndex         = mm + basisSize * ii;
        const int mmDofStartIndex     = mmNodeIndex * numberOfVariables;

        for(int ivar=0; ivar < numberOfVariables; ivar++) {
          lduh3D[ii][jj][ivar] += weight/dx[0] * kernels::Kxi[order][jj][mm] * lFhi_x[mmDofStartIndex+ivar];
        }
      }
    }
  }


  // Compute the "derivatives" (contributions of the stiffness matrix)
  // y direction (independent from the y and z derivatives)
  for (int ii=0; ii<basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {

      double weight = kernels::gaussLegendreWeights[order][jj];

      // MATMUL: Kxi * lFhi_y
      for(int mm=0; mm < basisSize; mm++) {
        // without reordering:
        //const int mmNodeIndex         = jj + basisSize * mm;
        //const int mmDofStartIndex     = mmNodeIndex * numberOfVariables;
        const int mmNodeIndex         = mm + basisSize * jj;
        const int mmDofStartIndex     = mmNodeIndex * numberOfVariables;

        // now we benefit from the reordering of lFhi_y
        for(int ivar=0; ivar < numberOfVariables; ivar++) {
          lduh3D[ii][jj][ivar] += weight/dx[1] * kernels::Kxi[order][ii][mm] * lFhi_y[mmDofStartIndex+ivar];
        }
      }
    }
  }
}

void kernels::aderdg::generic::surfaceIntegral(
    double * lduh,
    const double * const lFhbnd,
    const tarch::la::Vector<DIMENSIONS,double>& dx,
    const int numberOfVariables,
    const int basisSize
){
  // @todo Angelika
  // Please remove the typedefs in generic kernels again since numberOf(...)Dof is not
  // a compile time variable anymore

  // @todo 03/02/16:Dominic Etienne Charrier
  // Fixed bug: Number of variables was missing
  constexpr int numberOfFaceDof = 5 * (3+1);//numberOfVariables * tarch::la::aPowI(DIMENSIONS-1,basisSize);
  constexpr int order           = 3;

  const double * FLeft  = &(lFhbnd[EXAHYPE_FACE_LEFT  * numberOfFaceDof]);
  const double * FRight = &(lFhbnd[EXAHYPE_FACE_RIGHT * numberOfFaceDof]);
  const double * FFront = &(lFhbnd[EXAHYPE_FACE_FRONT * numberOfFaceDof]);
  const double * FBack  = &(lFhbnd[EXAHYPE_FACE_BACK  * numberOfFaceDof]);

  // access lduh(nDOF[2] x nDOF[1] x numberOfVariables) in the usual 3D array manner
  typedef double tensor_t[3+1][5];
  tensor_t *lduh3D = (tensor_t *)lduh;

  // x direction (independent from the y and z)
  for (int jj=0; jj<basisSize; jj++) {
    const int nodeIndex     = jj;
    const int dofStartIndex = nodeIndex * numberOfVariables;

    double weight =  kernels::gaussLegendreWeights[order][jj];

    for(int mm=0; mm < basisSize; mm++) {
      for(int ivar=0; ivar < numberOfVariables; ivar++) {
        lduh3D[jj][mm][ivar]
                       -=  weight/dx[0] * ( kernels::FRCoeff[order][mm] * FRight[dofStartIndex+ivar] - kernels::FLCoeff[order][mm] * FLeft[dofStartIndex+ivar] );
      }
    }
  }

  // Above seems okay!

  // y direction (independent from the y and z)
  for (int jj=0; jj<basisSize; jj++) {
    const int nodeIndex     = jj;
    const int dofStartIndex = nodeIndex * numberOfVariables;

    double weight =  kernels::gaussLegendreWeights[order][jj];

    for(int mm=0; mm < basisSize; mm++) {
      for(int ivar=0; ivar < numberOfVariables; ivar++) {
        lduh3D[mm][jj][ivar]
                       -=  weight/dx[1] * ( kernels::FRCoeff[order][mm] * FBack[dofStartIndex+ivar] - kernels::FLCoeff[order][mm] * FFront[dofStartIndex+ivar] );
      }
    }
  }
}

