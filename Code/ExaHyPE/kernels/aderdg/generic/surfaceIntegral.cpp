#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

using std::endl;
using std::cout;

extern "C" 
{
  void adersurfaceintegral_(double *lduh, double *lFhi, double *dx);
}


void kernels::aderdg::generic::surfaceIntegral(
    double * lduh,
    const double * const lFbnd,
    const tarch::la::Vector<DIMENSIONS,double>& dx,
    const int numberOfVariables,
    const int basisSize
){
#if DIMENSIONS == 2
  constexpr int numberOfFaceDof = 5 * (3+1);//numberOfVariables * tarch::la::aPowI(DIMENSIONS-1,basisSize);
  const double * FLeft  = &(lFbnd[EXAHYPE_FACE_LEFT  * numberOfFaceDof]);
  const double * FRight = &(lFbnd[EXAHYPE_FACE_RIGHT * numberOfFaceDof]);
  const double * FFront = &(lFbnd[EXAHYPE_FACE_FRONT * numberOfFaceDof]);
  const double * FBack  = &(lFbnd[EXAHYPE_FACE_BACK  * numberOfFaceDof]);

  surfaceIntegralXDirection( lduh,FLeft,  dx[0],0,+1.,numberOfVariables,basisSize );
  surfaceIntegralXDirection( lduh,FRight, dx[0],1,-1.,numberOfVariables,basisSize );
  surfaceIntegralYDirection( lduh,FFront, dx[1],0,+1.,numberOfVariables,basisSize );
  surfaceIntegralYDirection( lduh,FBack,  dx[1],1,-1.,numberOfVariables,basisSize );
  
#elif DIMENSIONS == 3

    double* lFbndFortran = new double[numberOfVariables*6*basisSize*basisSize];
    for(int i=0; i < numberOfVariables*6*basisSize*basisSize; i++){
      lFbndFortran[i] = -123.45;
    }
    
    // Permutation of lFbnd and lFbnd

    // for(int i=0; i < numberOfVariables*6*basisSize*basisSize; i++){
      // lFbndFortran[i] = i;
    // }
    
    for (int bb=0; bb<basisSize; bb++) {  // loop over dof
      for (int aa=0; aa<basisSize; aa++) {
        for(int ivar=0; ivar < numberOfVariables; ivar++) {
          for(int face=0; face < 6; face++) {
            lFbndFortran[f2p4(ivar, face, aa, bb)] = lFbnd[p2f4(ivar, face, aa, bb)];
          }
        }
      }
    }

    
    
    
    // std::ofstream ofs;
    // ofs.open ("boutput_lFbnd.txt", std::ofstream::out);
    // for (int ii=0; ii<numberOfVariables*6*basisSize*basisSize; ii++) {
      // ofs << lFbnd[ii] << "\n";
    // }
    // ofs.close();

    // ofs.open ("boutput_lFbndFortran.txt", std::ofstream::out);
    // for (int ii=0; ii<numberOfVariables*6*basisSize*basisSize; ii++) {
      // ofs << lFbndFortran[ii] << "\n";
    // }
    // ofs.close();
    

    double dxTemp[3] = {dx[0], dx[1], dx[2]};
    
    adersurfaceintegral_(lduh, lFbndFortran, &dxTemp[0]);

    // for(int i=0; i < numberOfVariables*basisSize*basisSize*basisSize; i++){
      // cout << lduh[i] << endl;
    // }
    // exit(0);
  
#endif
}

#if DIMENSIONS == 2
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
  for (int ii=0; ii<basisSize; ii++) {
    const int nodeIndex     = ii; 
    const int dofStartIndex = nodeIndex * numberOfVariables;

    double weight =  kernels::gaussLegendreWeights[order][ii];

    for(int mm=0; mm < basisSize; mm++) {
      for(int ivar=0; ivar < numberOfVariables; ivar++) {
        lduh3D[mm][ii][ivar] // direction dependent! 
                       += updateSign * weight/area * ( kernels::FCoeff[order][facePosition][mm] * lFhbnd[dofStartIndex+ivar] );
      }
    }
  }
}
#endif

