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

#include "../../Limiter.h"
 
#if DIMENSIONS == 3 
 
namespace kernels {
namespace limiter {
namespace generic {
namespace c {

        
//TODO JMG
double* getGaussLobattoData(const double* const luh, const int numberOfVariables, const int basisSize, int& basisSizeLob) {
  basisSizeLob = 0;
  
  return new double[1];
}

/*
REAL :: limz(nVar,nDOF(1),nDOF(2),nSubLimV(3)), limy(nVar,nDOF(1),nSubLimV(2),nSubLimV(3))
IF (nDim == 3) THEN
        DO jjj = 1, nDOF(2)  ! y
            DO iii = 1, nDOF(1) ! x
                limz(:,iii,jjj,:) = MATMUL( luh(:,iii,jjj,:), TRANSPOSE(uh2lim) ) 
            ENDDO
        ENDDO
    ELSE
        limz = 0. 
        DO jjj = 1, nDOF(2)  ! y
            DO iii = 1, nDOF(1) ! x
                limz(:,iii,jjj,1) = luh(:,iii,jjj,1) 
            ENDDO
        ENDDO
    ENDIF
    ! mapping of uh to the sub-cell limiter along y direction
    IF( nDim >= 2 ) THEN
        limy = 0. 
        DO jjj = 1, nSubLimV(3)  ! z 
            DO iii = 1, nDOF(1) ! x
                limy(:,iii,:,jjj) = MATMUL( limz(:,iii,:,jjj), TRANSPOSE(uh2lim) ) 
            ENDDO
        ENDDO
    ELSE
        limy = 0. 
        DO iii = 1, nDOF(1) ! x
            limy(:,iii,1,1) = luh(:,iii,1,1) 
        ENDDO
    ENDIF    
    ! mapping of uh to the sub-cell limiter along x direction 
    lim = 0. 
    DO jjj = 1, nSubLimV(3)  ! z 
        DO iii = 1, nSubLimV(2) ! y
            lim(:,:,iii,jjj) = MATMUL( limy(:,:,iii,jjj), TRANSPOSE(uh2lim) )  
        ENDDO
    ENDDO   
*/
double* getFVMData(const double* const luh, const int numberOfVariables, const int basisSize, int& basisSizeLim) {
  
  basisSizeLim = getLimBasisSize(basisSize);
  
  idx4 idxLuh(basisSize, basisSize, basisSize, numberOfVariables);
  idx2 idxConv(basisSize, basisSizeLim);
  
  double* lim = new double[basisSizeLim*basisSizeLim*basisSizeLim*numberOfVariables]; //Fortran ref: lim(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3))
  idx4 idxLim(basisSizeLim, basisSizeLim, basisSizeLim, numberOfVariables);
  
  double* tmpZ = new double[basisSizeLim*basisSize*basisSize*numberOfVariables]; //Fortran ref: limz(nVar,nDOF(1),nDOF(2),nSubLimV(3))
  idx4 idxZ(basisSizeLim, basisSize, basisSize, numberOfVariables);
  double* tmpY = new double[basisSizeLim*basisSizeLim*basisSize*numberOfVariables]; //Fortran ref: limy(nVar,nDOF(1),nSubLimV(2),nSubLimV(3))
  idx4 idxY(basisSizeLim, basisSizeLim, basisSize, numberOfVariables);
  int x,y,z,v,k;
  
  for(y=0; y<basisSize; y++) {
    for(x=0; x<basisSize; x++) {
      //limz(:,iii,jjj,:) = MATMUL( luh(:,iii,jjj,:), TRANSPOSE(uh2lim) ) 
      for(z=0; z<basisSizeLim; z++) {
        for(v=0; v<numberOfVariables; v++) {
          tmpZ[idxZ(z,y,x,v)] = 0;
          for(k=0; k<basisSize; k++) {
            tmpZ[idxZ(z,y,x,v)] += luh[idxLuh(k,y,x,v)] * uh2lim[idxConv(k,z)];
          }
        }
      }
    }
  }
  
  for(z=0; z<basisSizeLim; z++) {
    for(x=0; x<basisSize; x++) {
      //limy(:,iii,:,jjj) = MATMUL( limz(:,iii,:,jjj), TRANSPOSE(uh2lim) ) 
      for(y=0; y<basisSizeLim; y++) {
        for(v=0; v<numberOfVariables; v++) {
          tmpY[idxY(z,y,x,v)] = 0;
          for(k=0; k<basisSize; k++) {
            tmpY[idxY(z,y,x,v)] += tmpZ[idxZ(z,k,x,v)] * uh2lim[idxConv(k,y)];
          }
        }
      }
    }
  }
  
  for(z=0; z<basisSizeLim; z++) {
    for(y=0; y<basisSizeLim; y++) {
      //lim(:,:,iii,jjj) = MATMUL( limy(:,:,iii,jjj), TRANSPOSE(uh2lim) )
      for(x=0; x<basisSizeLim; x++) {
        for(v=0; v<numberOfVariables; v++) {
          lim[idxLim(z,y,x,v)] = 0;
          for(k=0; k<basisSize; k++) {
            lim[idxLim(z,y,x,v)] += tmpY[idxY(z,y,k,v)] * uh2lim[idxConv(k,x)];
          }
        }
      }
    }
  }
  delete[] tmpZ;
  delete[] tmpY;
  
  return lim;
}


} // namespace c
} // namespace generic
} // namespace limiter
} // namespace kernel

#endif //DIMENSIONS == 3