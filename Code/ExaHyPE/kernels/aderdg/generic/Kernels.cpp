#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"


void kernels::aderdg::generic::solutionUpdate(double * luh,
                                              const double * const lduh,
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

//void kernels::aderdg::generic::surfaceIntegral(
//    double * lduh,
//    const double * const lFhbnd,
//    const tarch::la::Vector<DIMENSIONS,double>& dx,
//    const int numberOfVariables,
//    const int basisSize
//){
//  // @todo Angelika
//  // Please remove the typedefs in generic kernels again since numberOf(...)Dof is not
//  // a compile time variable anymore
//
//  // @todo 03/02/16:Dominic Etienne Charrier
//  // Fixed bug: Number of variables was missing
//  constexpr int numberOfFaceDof = 5 * (3+1);//numberOfVariables * tarch::la::aPowI(DIMENSIONS-1,basisSize);
//  constexpr int order           = 3;
//
//  const double * FLeft  = &(lFhbnd[EXAHYPE_FACE_LEFT  * numberOfFaceDof]);
//  const double * FRight = &(lFhbnd[EXAHYPE_FACE_RIGHT * numberOfFaceDof]);
//  const double * FFront = &(lFhbnd[EXAHYPE_FACE_FRONT * numberOfFaceDof]);
//  const double * FBack  = &(lFhbnd[EXAHYPE_FACE_BACK  * numberOfFaceDof]);
//
//  // access lduh(nDOF[2] x nDOF[1] x numberOfVariables) in the usual 3D array manner
//  typedef double tensor_t[3+1][5];
//  tensor_t *lduh3D = (tensor_t *)lduh;
//
//  // x direction (independent from the y and z)
//  for (int jj=0; jj<basisSize; jj++) {
//    const int nodeIndex     = jj;
//    const int dofStartIndex = nodeIndex * numberOfVariables;
//
//    double weight =  kernels::gaussLegendreWeights[order][jj];
//
//    for(int mm=0; mm < basisSize; mm++) {
//      for(int ivar=0; ivar < numberOfVariables; ivar++) {
//        lduh3D[jj][mm][ivar]
//                       -=  weight/dx[0] * ( kernels::FRCoeff[order][mm] * FRight[dofStartIndex+ivar] - kernels::FLCoeff[order][mm] * FLeft[dofStartIndex+ivar] );
//      }
//    }
//  }
//
//  // Above seems okay!
//
//  // y direction (independent from the y and z)
//  for (int jj=0; jj<basisSize; jj++) {
//    const int nodeIndex     = jj;
//    const int dofStartIndex = nodeIndex * numberOfVariables;
//
//    double weight =  kernels::gaussLegendreWeights[order][jj];
//
//    for(int mm=0; mm < basisSize; mm++) {
//      for(int ivar=0; ivar < numberOfVariables; ivar++) {
//        lduh3D[mm][jj][ivar]
//                       -=  weight/dx[1] * ( kernels::FRCoeff[order][mm] * FBack[dofStartIndex+ivar] - kernels::FLCoeff[order][mm] * FFront[dofStartIndex+ivar] );
//      }
//    }
//  }
//}

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
  constexpr int dim         = DIMENSIONS;         // 2
  constexpr int dimTimesTwo = (2*DIMENSIONS);     // 4
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

void kernels::aderdg::generic::extrapolatedPredictor(
    double* lQhbnd,
    double* lFhbnd,
    const double* const lQhi,
    const double* const lFhi,
    const double predictorTimeStepSize,
    const int numberOfVariables,
    const int basisSize
){
  constexpr int BASISSIZE       = (3+1); // basisSize=order+1
  constexpr int NVAR            = 5;     // numberOfVariables
  constexpr int numberOfFaceDof = NVAR * BASISSIZE;// tarch::la::aPowI(dim-1,BASISSIZE);

  extrapolatedPredictorXDirection(&lQhbnd[EXAHYPE_FACE_LEFT * numberOfFaceDof],
                                  &lFhbnd[EXAHYPE_FACE_LEFT * numberOfFaceDof],
                                  lQhi,lFhi,
                                  0, // face position
                                  predictorTimeStepSize, // evaluationTimeStepSize
                                  predictorTimeStepSize,
                                  numberOfVariables,basisSize);

  extrapolatedPredictorXDirection(&lQhbnd[EXAHYPE_FACE_RIGHT * numberOfFaceDof],
                                  &lFhbnd[EXAHYPE_FACE_RIGHT * numberOfFaceDof],
                                  lQhi,lFhi,
                                  1, // face position
                                  predictorTimeStepSize, // evaluationTimeStepSize
                                  predictorTimeStepSize,
                                  numberOfVariables,basisSize);

  extrapolatedPredictorYDirection(&lQhbnd[EXAHYPE_FACE_FRONT * numberOfFaceDof],
                                  &lFhbnd[EXAHYPE_FACE_FRONT * numberOfFaceDof],
                                  lQhi,lFhi,
                                  0, // face position
                                  predictorTimeStepSize, // evaluationTimeStepSize
                                  predictorTimeStepSize,
                                  numberOfVariables,basisSize);

  extrapolatedPredictorYDirection(&lQhbnd[EXAHYPE_FACE_BACK * numberOfFaceDof],
                                  &lFhbnd[EXAHYPE_FACE_BACK * numberOfFaceDof],
                                  lQhi,lFhi,
                                  1, // face position
                                  predictorTimeStepSize, // evaluationTimeStepSize
                                  predictorTimeStepSize,
                                  numberOfVariables,basisSize);
}

void kernels::aderdg::generic::extrapolatedPredictorXDirection(
    double* lQhbnd,
    double* lFhbnd,
    const double* const lQhi,
    const double* const lFhi,
    const int facePosition,
    const double evaluationTimeStepSize,
    const double predictorTimeStepSize,
    const int numberOfVariables,
    const int basisSize
){
  constexpr int direction = 0;
  constexpr int dim                = DIMENSIONS;         // 2
  constexpr int dimTimesTwo        = (2*DIMENSIONS);     // 4
  constexpr int BASISSIZE          = (3+1); // basisSize=order+1
  constexpr int NVAR               = 5;     // numberOfVariables

  /////////////////////////////////////////////////
  // Compute the bounday-extrapolated values for Q and F*n
  /////////////////////////////////////////////////
  constexpr int numberOfDof     = NVAR * BASISSIZE * BASISSIZE;/*tarch::la::aPowI(dim,BASISSIZE);*/
  constexpr int numberOfFaceDof = NVAR * BASISSIZE;// tarch::la::aPowI(dim-1,BASISSIZE);

  // todo 10/02/16:Dominic Etienne Charrier
  // REMOVED:
  /*
      memset((double *) &lQhbnd[0],0,sizeof(double) * numberOfFaceDof * dimTimesTwo);
      memset((double *) &lFhbnd[0],0,sizeof(double) * numberOfFaceDof * dimTimesTwo);
   */
  // ADDED:
  memset(lQhbnd,0,sizeof(double) * numberOfFaceDof);
  memset(lFhbnd,0,sizeof(double) * numberOfFaceDof);

  const double* Q;
  const double* lFhi_n = &lFhi[direction * numberOfDof];

  for (int jj=0; jj<BASISSIZE; jj++) {
    const int nodeIndex      = jj;
    const int dofStartIndex  = nodeIndex * NVAR;

    for (int mm=0; mm<BASISSIZE; mm++) { // loop over dof
      const int mmNodeIndex         = mm + BASISSIZE * jj; // dimension depending line
      const int mmDofStartIndex     = mmNodeIndex * NVAR;

      Q = &lQhi[mmDofStartIndex];

      for(int ivar=0; ivar < NVAR; ivar++) {

        // todo 10/02/16:Dominic Etienne Charrier
        // REMOVED:
        /*
        lQhbnd[dofStartIndexL+ivar] += kernels::FLCoeff[BASISSIZE-1][mm] * Q[ivar]; // lQhbnd(Facei,nDOF,nVar)
        lQhbnd[dofStartIndexR+ivar] += kernels::FRCoeff[BASISSIZE-1][mm] * Q[ivar]; // lQhbnd(Facei,nDOF,nVar)

        lFhbnd[dofStartIndexL+ivar] += kernels::FLCoeff[BASISSIZE-1][mm] * lFhi_x[mmDofStartIndex+ivar]; // lFhi_x * FLCoeff
        lFhbnd[dofStartIndexR+ivar] += kernels::FRCoeff[BASISSIZE-1][mm] * lFhi_x[mmDofStartIndex+ivar]; // lFhi_x * FRCoeff
         */
        // todo 10/02/16:Dominic Etienne Charrier
        // ADDED:
        lQhbnd[dofStartIndex+ivar] += kernels::FCoeff[BASISSIZE-1][facePosition][mm] * Q     [ivar]; // lQhbnd(Facei,nDOF,nVar)
        lFhbnd[dofStartIndex+ivar] += kernels::FCoeff[BASISSIZE-1][facePosition][mm] * lFhi_n[mmDofStartIndex+ivar]; // lFhi_x * FLCoeff
      }
    }
    continue;
  }
}

void kernels::aderdg::generic::extrapolatedPredictorYDirection(
    double* lQhbnd,
    double* lFhbnd,
    const double* const lQhi,
    const double* const lFhi,
    const int facePosition,
    const double evaluationTimeStepSize,
    const double predictorTimeStepSize,
    const int numberOfVariables,
    const int basisSize
){
  constexpr int direction       = 1;
  constexpr int dim             = DIMENSIONS;         // 2
  constexpr int dimTimesTwo     = (2*DIMENSIONS);     // 4
  constexpr int BASISSIZE       = (3+1); // basisSize=order+1
  constexpr int NVAR            = 5;     // numberOfVariables

  /////////////////////////////////////////////////
  // Compute the bounday-extrapolated values for Q and F*n
  /////////////////////////////////////////////////
  constexpr int numberOfDof     = NVAR * BASISSIZE * BASISSIZE;/*tarch::la::aPowI(dim,BASISSIZE);*/
  constexpr int numberOfFaceDof = NVAR * BASISSIZE;// tarch::la::aPowI(dim-1,BASISSIZE);

  // todo 10/02/16:Dominic Etienne Charrier
  // REMOVED:
  /*
      memset((double *) &lQhbnd[0],0,sizeof(double) * numberOfFaceDof * dimTimesTwo);
      memset((double *) &lFhbnd[0],0,sizeof(double) * numberOfFaceDof * dimTimesTwo);
   */
  // ADDED:
  memset(lQhbnd,0,sizeof(double) * numberOfFaceDof);
  memset(lFhbnd,0,sizeof(double) * numberOfFaceDof);

  const double* Q;
  const double* lFhi_n = &lFhi[direction * numberOfDof];

  for (int jj=0; jj<BASISSIZE; jj++) {
    const int nodeIndex      = jj;
    const int dofStartIndex  = nodeIndex * NVAR;

    for (int mm=0; mm<BASISSIZE; mm++) { // loop over dof
      const int mmNodeIndex         = mm * BASISSIZE + jj; // face depending line
      const int mmDofStartIndex     = mmNodeIndex * NVAR;

      Q = &lQhi[mmDofStartIndex];

      for(int ivar=0; ivar < NVAR; ivar++) {

        // todo 10/02/16:Dominic Etienne Charrier
        // REMOVED:
        /*
        lQhbnd[dofStartIndexL+ivar] += kernels::FLCoeff[BASISSIZE-1][mm] * Q[ivar]; // lQhbnd(Facei,nDOF,nVar)
        lQhbnd[dofStartIndexR+ivar] += kernels::FRCoeff[BASISSIZE-1][mm] * Q[ivar]; // lQhbnd(Facei,nDOF,nVar)

        lFhbnd[dofStartIndexL+ivar] += kernels::FLCoeff[BASISSIZE-1][mm] * lFhi_x[mmDofStartIndex+ivar]; // lFhi_x * FLCoeff
        lFhbnd[dofStartIndexR+ivar] += kernels::FRCoeff[BASISSIZE-1][mm] * lFhi_x[mmDofStartIndex+ivar]; // lFhi_x * FRCoeff
         */
        // todo 10/02/16:Dominic Etienne Charrier
        // ADDED:
        lQhbnd[dofStartIndex+ivar] += kernels::FCoeff[BASISSIZE-1][facePosition][mm] * Q     [ivar]; // lQhbnd(Facei,nDOF,nVar)
        lFhbnd[dofStartIndex+ivar] += kernels::FCoeff[BASISSIZE-1][facePosition][mm] * lFhi_n[mmDofStartIndex+ivar]; // lFhi_x * FLCoeff
      }
    }
    continue;
  }
}
