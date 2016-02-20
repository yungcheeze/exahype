#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"


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
#if DIMENSIONS == 2
  constexpr int numberOfDof = 5 * (3+1)*(3+1); //tarch::la::aPowI(3+1,DIMENSIONS);
#else
  constexpr int numberOfDof = 5 * (3+1)*(3+1)*(3+1); //tarch::la::aPowI(3+1,DIMENSIONS);
#endif  
  constexpr int order       = 3;

#if DIMENSIONS == 2
  // memory layout of lFhi:
  // lFhi = [ lFhi_x | lFhi_y ] ordered as follows
  // (a) lFhi_x[nDOF_y][nDOF_x][nVar]
  // (b) lFhi_y[nDOF_x][nDOF_y][nVar]
  // Note the variable order of lFhi_y
  // let's not bother with offsets and define separate flux matrices
  const double * lFhi_x = &lFhi[0];            // f flux
  const double * lFhi_y = &lFhi[numberOfDof];  // g flux

  memset(lduh,0,sizeof(double) * numberOfDof);
#else
  const double * lFhi_x = &lFhi[0            ];  // f flux
  const double * lFhi_y = &lFhi[numberOfDof  ];  // g flux
  const double * lFhi_z = &lFhi[numberOfDof*2];  // h flux

  memset(lduh,0,sizeof(double) * numberOfDof);
#endif  


#if DIMENSIONS == 2
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
#else
  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    for (int jj=0; jj<basisSize; jj++) {
      for (int kk=0; kk<basisSize; kk++) { // loop over dof

      double weight = kernels::gaussLegendreWeights[basisSize-1][jj] *
                        kernels::gaussLegendreWeights[basisSize-1][kk];


        for(int mm=0; mm < basisSize; mm++) {
          const int mmNodeIndex         = mm + basisSize * jj + basisSize * basisSize * kk;
          const int mmDofStartIndex     = mmNodeIndex * numberOfVariables;

          for(int ivar=0; ivar < numberOfVariables; ivar++) {
            lduh[mmDofStartIndex+ivar] += weight/dx[0] * kernels::Kxi[basisSize-1][ii][mm] * lFhi_x[mmNodeIndex+ivar];
            //todo VV
            lduh[mmDofStartIndex+ivar] += weight/dx[0] * kernels::Kxi[basisSize-1][ii][mm] * lFhi_x[mmNodeIndex+ivar];
          }
        }
      }
    }
  }
  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    for (int jj=0; jj<basisSize; jj++) {
      for (int kk=0; kk<basisSize; kk++) { // loop over dof

      double weight = kernels::gaussLegendreWeights[basisSize-1][jj] *
                        kernels::gaussLegendreWeights[basisSize-1][kk];


        for(int mm=0; mm < basisSize; mm++) {
          const int mmNodeIndex         = ii + basisSize * mm + basisSize * basisSize * kk;
          const int mmDofStartIndex     = mmNodeIndex * numberOfVariables;

          for(int ivar=0; ivar < numberOfVariables; ivar++) {
            lduh[mmDofStartIndex+ivar] += weight/dx[1] * kernels::Kxi[basisSize-1][jj][mm] * lFhi_y[mmNodeIndex+ivar];
            //todo VV
            lduh[mmDofStartIndex+ivar] += weight/dx[1] * kernels::Kxi[basisSize-1][jj][mm] * lFhi_y[mmNodeIndex+ivar];
          }
        }
      }
    }
  }
  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    for (int jj=0; jj<basisSize; jj++) {
      for (int kk=0; kk<basisSize; kk++) { // loop over dof

      double weight = kernels::gaussLegendreWeights[basisSize-1][jj] *
                        kernels::gaussLegendreWeights[basisSize-1][kk];


        for(int mm=0; mm < basisSize; mm++) {
          const int mmNodeIndex         = ii + basisSize * jj + basisSize * basisSize * mm;
          const int mmDofStartIndex     = mmNodeIndex * numberOfVariables;

          // TODO please check this again
          // something is definitely wrong here. lduh and lFhi_z are accessed in the same way
          // but are not stored in the same way
          // lFhi_z[(ii*basisSize*basisSize+basisSize*jj+mm)*numberOfVariables+ivar] ???
          for(int ivar=0; ivar < numberOfVariables; ivar++) {
            lduh[mmDofStartIndex+ivar] += weight/dx[2] * kernels::Kxi[basisSize-1][kk][mm] * lFhi_z[mmNodeIndex+ivar];
            //todo VV
            lduh[mmDofStartIndex+ivar] += weight/dx[2] * kernels::Kxi[basisSize-1][kk][mm] * lFhi_z[mmNodeIndex+ivar];
          }
        }
      }
    }
  }
 
#endif  
}

