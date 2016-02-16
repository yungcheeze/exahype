#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"


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

  //  surfaceIntegralXDirection( lduh,FLeft, dx[0],0,+1.,numberOfVariables,basisSize );
  //  surfaceIntegralXDirection( lduh,FRight,dx[0],1,-1.,numberOfVariables,basisSize );

  // x direction (independent from the y and z)
  for (int jj=0; jj<basisSize; jj++) {
    const int nodeIndex     = jj;
    const int dofStartIndex = nodeIndex * numberOfVariables;

    double weight =  kernels::gaussLegendreWeights[order][jj];

    for(int mm=0; mm < basisSize; mm++) {
      for(int ivar=0; ivar < numberOfVariables; ivar++) {
        int dimension     = 1;
        double updateSign = -1.;
        const double * F2   = &FRight[0];
        lduh3D[jj][mm][ivar] += updateSign * weight/dx[0] * ( kernels::FCoeff[order][dimension][mm] * F2 [dofStartIndex+ivar] );

        dimension  =  0;
        updateSign = +1.;
        F2     = &FLeft[0];
        lduh3D[jj][mm][ivar] += updateSign * weight/dx[0] * ( kernels::FCoeff[order][dimension][mm] * F2 [dofStartIndex+ivar] );
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
        int dimension     =  1;
        double updateSign = -1.;
        const double * F2   = &FFront[0];
        lduh3D[mm][jj][ivar] += updateSign * weight/dx[1] * ( kernels::FCoeff[order][dimension][mm] * F2 [dofStartIndex+ivar] );

        dimension  =  0;
        updateSign = +1.;
        F2         = &FFront[0];
        lduh3D[mm][jj][ivar] += updateSign * weight/dx[1] * ( kernels::FCoeff[order][dimension][mm] * F2 [dofStartIndex+ivar] );
      }
    }
  }
}

void kernels::aderdg::generic::surfaceIntegral2(
    double * lduh,
    const double * const lFhbnd,
    const tarch::la::Vector<DIMENSIONS,double>& dx,
    const int numberOfVariables,
    const int basisSize
){
  constexpr int numberOfFaceDof = 5 * (3+1);//numberOfVariables * tarch::la::aPowI(DIMENSIONS-1,basisSize);
  const double * FLeft  = &(lFhbnd[EXAHYPE_FACE_LEFT  * numberOfFaceDof]);
  const double * FRight = &(lFhbnd[EXAHYPE_FACE_RIGHT * numberOfFaceDof]);
  const double * FFront = &(lFhbnd[EXAHYPE_FACE_FRONT * numberOfFaceDof]);
  const double * FBack  = &(lFhbnd[EXAHYPE_FACE_BACK  * numberOfFaceDof]);

  surfaceIntegralXDirection( lduh,FLeft, dx[0],0,+1.,numberOfVariables,basisSize );
  surfaceIntegralXDirection( lduh,FRight,dx[0],1,-1.,numberOfVariables,basisSize );
  surfaceIntegralYDirection( lduh,FFront,dx[1],0,+1.,numberOfVariables,basisSize );
  surfaceIntegralYDirection( lduh,FBack, dx[1],1,-1.,numberOfVariables,basisSize );
}

void kernels::aderdg::generic::surfaceIntegralXDirection(
    double * lduh,
    const double * const lFhbnd,
    const double area,
    const int facePosition,  // 0 for "left" face, 1 for "right" face.
    const double updateSign, // -1 for "left" face, 1 for "right" face.
    const int numberOfVariables,
    const int basisSize
) {
  // @todo Angelika
  // Please remove the typedefs in generic kernels again since numberOf(...)Dof is not
  // a compile time variable anymore
  constexpr int numberOfFaceDof = 5 * (3+1);//numberOfVariables * tarch::la::aPowI(DIMENSIONS-1,basisSize);
  constexpr int order           = 3;

  // access lduh(nDOF[2] x nDOF[1] x numberOfVariables) in the usual 3D array manner
  typedef double tensor_t[3+1][5];
  tensor_t *lduh3D = (tensor_t *)lduh;

  // x direction (independent from the y and z)
  for (int jj=0; jj<basisSize; jj++) {
    const int nodeIndex     = jj;
    const int dofStartIndex = nodeIndex * numberOfVariables;

    double weight = kernels::gaussLegendreWeights[order][jj];

    for(int mm=0; mm < basisSize; mm++) {
      for(int ivar=0; ivar < numberOfVariables; ivar++) {
        lduh3D[jj][mm][ivar] // direction dependent!
                       += updateSign * weight/area * ( kernels::FCoeff[order][facePosition][mm] * lFhbnd [dofStartIndex+ivar] );
      }
    }
  }
}

void kernels::aderdg::generic::surfaceIntegralYDirection(
    double * lduh,
    const double * const lFhbnd,
    const double area,
    const int facePosition,  // 0 for "left" face, 1 for "right" face.
    const double updateSign, // -1 for "left" face, 1 for "right" face.
    const int numberOfVariables,
    const int basisSize
) {
  // @todo Angelika
  // Please remove the typedefs in generic kernels again since numberOf(...)Dof is not
  // a compile time variable anymore
  constexpr int numberOfFaceDof = 5 * (3+1);//numberOfVariables * tarch::la::aPowI(DIMENSIONS-1,basisSize);
  constexpr int order           = 3;

  // access lduh(nDOF[2] x nDOF[1] x numberOfVariables) in the usual 3D array manner
  typedef double tensor_t[3+1][5];
  tensor_t *lduh3D = (tensor_t *)lduh;

  // y direction (independent from the y and z)
  for (int jj=0; jj<basisSize; jj++) {
    const int nodeIndex     = jj;
    const int dofStartIndex = nodeIndex * numberOfVariables;

    double weight =  kernels::gaussLegendreWeights[order][jj];

    for(int mm=0; mm < basisSize; mm++) {
      for(int ivar=0; ivar < numberOfVariables; ivar++) {
        lduh3D[mm][jj][ivar] // direction dependent!
                       += updateSign * weight/area * ( kernels::FCoeff[order][facePosition][mm] * lFhbnd[dofStartIndex+ivar] );
      }
    }
  }
}


