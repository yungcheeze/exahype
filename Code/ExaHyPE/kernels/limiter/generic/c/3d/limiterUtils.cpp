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
  
double* uh2lim;
double* uh2lob;
double* lim2uh;

int getLimBasisSize(const int basisSize) {
  return 2*(basisSize-1) +1;
}

// exist already ?
// xin = xiGPN = gaussLegendreNodes[basisSize-1]
void BaseFunc1D(double* phi, double xi, const int basisSize) {
  int i,j,m;
  for(i=0; i<basisSize; i++) {
    phi[i] = 1;
  }
  for(m=0; m<basisSize; m++) {
    for(j=0; j<basisSize; j++) {
      if(j == m) continue;
      phi[m] = phi[m]*(xi- gaussLegendreNodes[basisSize-1][j])/(gaussLegendreNodes[basisSize-1][m]-gaussLegendreNodes[basisSize-1][j]);
    }
  }
}


void initProjectionMatrices(const int basisSize) {
  int i,j,k;
  double* phi = new double[basisSize];
  
  {
    int basisSizeLim = getLimBasisSize(basisSize);
    idx2 idx(basisSize, basisSizeLim);
    uh2lim = new double[basisSize*basisSizeLim](); //initialized at 0, Fortran ref: uh2lim(nSubLim,N+1)
    const double dxi = 1.0 / basisSizeLim;
    double xLeft, xRight, xi;
    for(i=0; i<basisSizeLim; i++) {
      xLeft = i*dxi;
      xRight = xLeft+dxi;
      for(j=0; j<basisSize; j++) {
        xi = xLeft + dxi*gaussLegendreNodes[basisSize-1][j];
        BaseFunc1D(phi, xi, basisSize);
        for(k=0; k<basisSize; k++) { //
          uh2lim[idx(k,i)] += gaussLegendreWeights[basisSize-1][j]*phi[k];
        }
      }
    }
  }
  
  //TODO JMG init uh2lob and lim2uh
  
  delete[] phi;
}

void freeProjectionMatrices(const int basisSize) {
  delete[] uh2lim;
}
        
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

/**
 * localMin, localMax are double[numberOfVariables]
 */
void findCellLocallocalMinlocalMax(const double* const luh, const int numberOfVariables, const int basisSize, double* localMin, double* localMax) {      

  int index, ii, iVar;
  
  // initialize and process luh
  index = 0;
  for(int iVar = 0; iVar < numberOfVariables; iVar++) {
    localMin[iVar] = luh[index];
    localMax[iVar] = luh[index];   
    index++;
  }
  for(ii = 1; ii < basisSize*basisSize*basisSize; ii++) {
    for(int iVar = 0; iVar < numberOfVariables; iVar++) {
      if(luh[index] < localMin[iVar]) {
        localMin[iVar] = luh[index];
      } else if(luh[index] > localMax[iVar]) {
        localMax[iVar] = luh[index];
      }    
      index++;
    }
  }
 
  // process lob
  int basisSizeLob = 0;
  double* lob = getGaussLobattoData(luh, numberOfVariables, basisSize, basisSizeLob);
  index = 0;
  for(ii = 0; ii < basisSizeLob*basisSizeLob*basisSizeLob; ii++) {
    for(int iVar = 0; iVar < numberOfVariables; iVar++) {
      if(lob[index] < localMin[iVar]) {
        localMin[iVar] = lob[index];
      } else if(lob[index] > localMax[iVar]) {
        localMax[iVar] = lob[index];
      }    
      index++;
    }
  }
  delete[] lob;
  
  // process lim
  index = 0;
  int basisSizeLim = 0;
  double* lim = getFVMData(luh, numberOfVariables, basisSize, basisSizeLim);
  for(ii = 0; ii < basisSizeLim*basisSizeLim*basisSizeLim; ii++) {
    for(int iVar = 0; iVar < numberOfVariables; iVar++) {
      if(lim[index] < localMin[iVar]) {
        localMin[iVar] = lim[index];
      } else if(lim[index] > localMax[iVar]) {
        localMax[iVar] = lim[index];
      }    
      index++;
    }
  }
  delete[] lim;

}


bool isTroubledCell(const double* const luh, const int numberOfVariables, const int basisSize, const double* const troubledMin, const double* const troubledMax) {
  
  double minMarginOfError = 0.0001;
  double diffScaling = 0.001;
  
  double* localMin = new double[numberOfVariables];
  double* localMax = new double[numberOfVariables];
  findCellLocallocalMinlocalMax(luh, numberOfVariables, basisSize, localMin, localMax);
  
  double ldiff;

  for(int iVar = 0; iVar < numberOfVariables; iVar++) {
    ldiff = std::max((troubledMax[iVar] - troubledMin[iVar]) * diffScaling, minMarginOfError);
    if((localMin[iVar] < (troubledMin[iVar] - ldiff)) || (localMax[iVar] > (troubledMax[iVar] + ldiff))) {
      return true;
    }
  }
  
  //TODO JMG (todo or not needed???) check PDE positivity and lim data not NAN 
  
  return false;
}

} // namespace c
} // namespace generic
} // namespace limiter
} // namespace kernel

#endif //DIMENSIONS == 3