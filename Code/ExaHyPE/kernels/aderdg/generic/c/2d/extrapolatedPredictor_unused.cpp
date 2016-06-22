/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 

#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

#if DIMENSIONS == 2

void kernels::aderdg::generic::c::extrapolatedPredictor(
    double* lQhbnd, double* lFhbnd, const double* const lQhi,
    const double* const lFhi, const double predictorTimeStepSize,
    const int numberOfVariables, const int basisSize) {
  constexpr int BASISSIZE = (3 + 1);  // basisSize=order+1
  constexpr int NVAR = 5;             // numberOfVariables
  constexpr int numberOfFaceDof =
      NVAR * BASISSIZE;  // tarch::la::aPowI(dim-1,BASISSIZE);

  extrapolatedPredictorXDirection(
      &lQhbnd[EXAHYPE_FACE_LEFT * numberOfFaceDof],
      &lFhbnd[EXAHYPE_FACE_LEFT * numberOfFaceDof], lQhi, lFhi,
      0,                      // face position
      predictorTimeStepSize,  // evaluationTimeStepSize
      predictorTimeStepSize, numberOfVariables, basisSize);

  extrapolatedPredictorXDirection(
      &lQhbnd[EXAHYPE_FACE_RIGHT * numberOfFaceDof],
      &lFhbnd[EXAHYPE_FACE_RIGHT * numberOfFaceDof], lQhi, lFhi,
      1,                      // face position
      predictorTimeStepSize,  // evaluationTimeStepSize
      predictorTimeStepSize, numberOfVariables, basisSize);

  extrapolatedPredictorYDirection(
      &lQhbnd[EXAHYPE_FACE_FRONT * numberOfFaceDof],
      &lFhbnd[EXAHYPE_FACE_FRONT * numberOfFaceDof], lQhi, lFhi,
      0,                      // face position
      predictorTimeStepSize,  // evaluationTimeStepSize
      predictorTimeStepSize, numberOfVariables, basisSize);

  extrapolatedPredictorYDirection(
      &lQhbnd[EXAHYPE_FACE_BACK * numberOfFaceDof],
      &lFhbnd[EXAHYPE_FACE_BACK * numberOfFaceDof], lQhi, lFhi,
      1,                      // face position
      predictorTimeStepSize,  // evaluationTimeStepSize
      predictorTimeStepSize, numberOfVariables, basisSize);
}

void kernels::aderdg::generic::c::extrapolatedPredictorXDirection(
    double* lQhbnd, double* lFhbnd, const double* const lQhi,
    const double* const lFhi, const int facePosition,
    const double evaluationTimeStepSize, const double predictorTimeStepSize,
    const int numberOfVariables, const int basisSize) {
  constexpr int direction = 0;
  // constexpr int dim                = 2;         // 2
  // constexpr int dimTimesTwo        = (2*2);     // 4
  constexpr int BASISSIZE = (3 + 1);  // basisSize=order+1
  constexpr int NVAR = 5;             // numberOfVariables

  /////////////////////////////////////////////////
  // Compute the bounday-extrapolated values for Q and F*n
  /////////////////////////////////////////////////
  constexpr int numberOfDof =
      NVAR * BASISSIZE * BASISSIZE; /*tarch::la::aPowI(dim,BASISSIZE);*/
  constexpr int numberOfFaceDof =
      NVAR * BASISSIZE;  // tarch::la::aPowI(dim-1,BASISSIZE);

  // todo 10/02/16:Dominic Etienne Charrier
  // REMOVED:
  /*
      memset((double *) &lQhbnd[0],0,sizeof(double) * numberOfFaceDof *
     dimTimesTwo);
      memset((double *) &lFhbnd[0],0,sizeof(double) * numberOfFaceDof *
     dimTimesTwo);
   */
  // ADDED:
  memset(lQhbnd, 0, sizeof(double) * numberOfFaceDof);
  memset(lFhbnd, 0, sizeof(double) * numberOfFaceDof);

  const double* Q;
  const double* lFhi_n = &lFhi[direction * numberOfDof];

  for (int jj = 0; jj < BASISSIZE; jj++) {
    const int nodeIndex = jj;
    const int dofStartIndex = nodeIndex * NVAR;

    for (int mm = 0; mm < BASISSIZE; mm++) {        // loop over dof
      const int mmNodeIndex = mm + BASISSIZE * jj;  // dimension depending line
      const int mmDofStartIndex = mmNodeIndex * NVAR;

      Q = &lQhi[mmDofStartIndex];

      for (int ivar = 0; ivar < NVAR; ivar++) {
        // todo 10/02/16:Dominic Etienne Charrier
        // REMOVED:
        /*
        lQhbnd[dofStartIndexL+ivar] += kernels::FLCoeff[BASISSIZE-1][mm] *
        Q[ivar]; // lQhbnd(Facei,nDOF,nVar)
        lQhbnd[dofStartIndexR+ivar] += kernels::FRCoeff[BASISSIZE-1][mm] *
        Q[ivar]; // lQhbnd(Facei,nDOF,nVar)

        lFhbnd[dofStartIndexL+ivar] += kernels::FLCoeff[BASISSIZE-1][mm] *
        lFhi_x[mmDofStartIndex+ivar]; // lFhi_x * FLCoeff
        lFhbnd[dofStartIndexR+ivar] += kernels::FRCoeff[BASISSIZE-1][mm] *
        lFhi_x[mmDofStartIndex+ivar]; // lFhi_x * FRCoeff
         */
        // todo 10/02/16:Dominic Etienne Charrier
        // ADDED:
        lQhbnd[dofStartIndex + ivar] +=
            kernels::FCoeff[BASISSIZE - 1][facePosition][mm] *
            Q[ivar];  // lQhbnd(Facei,nDOF,nVar)
        lFhbnd[dofStartIndex + ivar] +=
            kernels::FCoeff[BASISSIZE - 1][facePosition][mm] *
            lFhi_n[mmDofStartIndex + ivar];  // lFhi_x * FLCoeff
      }
    }
    continue;
  }
}

void kernels::aderdg::generic::c::extrapolatedPredictorYDirection(
    double* lQhbnd, double* lFhbnd, const double* const lQhi,
    const double* const lFhi, const int facePosition,
    const double evaluationTimeStepSize, const double predictorTimeStepSize,
    const int numberOfVariables, const int basisSize) {
  constexpr int direction = 1;
  // constexpr int dim             = 2;         // 2
  // constexpr int dimTimesTwo     = (2*2);     // 4
  constexpr int BASISSIZE = (3 + 1);  // basisSize=order+1
  constexpr int NVAR = 5;             // numberOfVariables

  /////////////////////////////////////////////////
  // Compute the bounday-extrapolated values for Q and F*n
  /////////////////////////////////////////////////
  constexpr int numberOfDof =
      NVAR * BASISSIZE * BASISSIZE; /*tarch::la::aPowI(dim,BASISSIZE);*/
  constexpr int numberOfFaceDof =
      NVAR * BASISSIZE;  // tarch::la::aPowI(dim-1,BASISSIZE);

  // todo 10/02/16:Dominic Etienne Charrier
  // REMOVED:
  /*
      memset((double *) &lQhbnd[0],0,sizeof(double) * numberOfFaceDof *
     dimTimesTwo);
      memset((double *) &lFhbnd[0],0,sizeof(double) * numberOfFaceDof *
     dimTimesTwo);
   */
  // ADDED:
  memset(lQhbnd, 0, sizeof(double) * numberOfFaceDof);
  memset(lFhbnd, 0, sizeof(double) * numberOfFaceDof);

  const double* Q;
  const double* lFhi_n = &lFhi[direction * numberOfDof];

  for (int jj = 0; jj < BASISSIZE; jj++) {
    const int nodeIndex = jj;
    const int dofStartIndex = nodeIndex * NVAR;

    for (int mm = 0; mm < BASISSIZE; mm++) {        // loop over dof
      const int mmNodeIndex = mm * BASISSIZE + jj;  // face depending line
      const int mmDofStartIndex = mmNodeIndex * NVAR;

      Q = &lQhi[mmDofStartIndex];

      for (int ivar = 0; ivar < NVAR; ivar++) {
        // todo 10/02/16:Dominic Etienne Charrier
        // REMOVED:
        /*
        lQhbnd[dofStartIndexL+ivar] += kernels::FLCoeff[BASISSIZE-1][mm] *
        Q[ivar]; // lQhbnd(Facei,nDOF,nVar)
        lQhbnd[dofStartIndexR+ivar] += kernels::FRCoeff[BASISSIZE-1][mm] *
        Q[ivar]; // lQhbnd(Facei,nDOF,nVar)

        lFhbnd[dofStartIndexL+ivar] += kernels::FLCoeff[BASISSIZE-1][mm] *
        lFhi_x[mmDofStartIndex+ivar]; // lFhi_x * FLCoeff
        lFhbnd[dofStartIndexR+ivar] += kernels::FRCoeff[BASISSIZE-1][mm] *
        lFhi_x[mmDofStartIndex+ivar]; // lFhi_x * FRCoeff
         */
        // todo 10/02/16:Dominic Etienne Charrier
        // ADDED:
        lQhbnd[dofStartIndex + ivar] +=
            kernels::FCoeff[BASISSIZE - 1][facePosition][mm] *
            Q[ivar];  // lQhbnd(Facei,nDOF,nVar)
        lFhbnd[dofStartIndex + ivar] +=
            kernels::FCoeff[BASISSIZE - 1][facePosition][mm] *
            lFhi_n[mmDofStartIndex + ivar];  // lFhi_x * FLCoeff
      }
    }
    continue;
  }
}
#endif 