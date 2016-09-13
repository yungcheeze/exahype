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

//Fortran (Limiter.f90): GetLobattoData
double* getGaussLobattoData(const double* const luh, const int numberOfVariables, const int basisSize) {
  
  idx4 idx(basisSize, basisSize, basisSize, numberOfVariables);
  idx2 idxConv(basisSize, basisSize);
  
  double* lob = new double[basisSize*basisSize*basisSize*numberOfVariables]; //Fortran ref: lob(nVar,nDOF(1),nDOF(2),nDOF(3))
  
  double* tmpZ = new double[basisSize*basisSize*basisSize*numberOfVariables]; //Fortran ref: lobz(nVar,nDOF(1),nDOF(2),nDOF(3))
  double* tmpY = new double[basisSize*basisSize*basisSize*numberOfVariables]; //Fortran ref: loby(nVar,nDOF(1),nDOF(2),nDOF(3))
  int x,y,z,v,k;
  
  for(y=0; y<basisSize; y++) {
    for(x=0; x<basisSize; x++) {
      //lobz(:,iii,jjj,:) = MATMUL( luh(:,iii,jjj,:), TRANSPOSE(uh2lob) )
      for(z=0; z<basisSize; z++) {
        for(v=0; v<numberOfVariables; v++) {
          tmpZ[idx(z,y,x,v)] = 0;
          for(k=0; k<basisSize; k++) {
            tmpZ[idx(z,y,x,v)] += luh[idx(k,y,x,v)] * uh2lob[basisSize-1][idxConv(k,z)];
          }
        }
      }
    }
  }
  
  for(z=0; z<basisSize; z++) {
    for(x=0; x<basisSize; x++) {
      //loby(:,iii,:,jjj) = MATMUL( lobz(:,iii,:,jjj), TRANSPOSE(uh2lob) ) 
      for(y=0; y<basisSize; y++) {
        for(v=0; v<numberOfVariables; v++) {
          tmpY[idx(z,y,x,v)] = 0;
          for(k=0; k<basisSize; k++) {
            tmpY[idx(z,y,x,v)] += tmpZ[idx(z,k,x,v)] * uh2lob[basisSize-1][idxConv(k,y)];
          }
        }
      }
    }
  }
  
  for(z=0; z<basisSize; z++) {
    for(y=0; y<basisSize; y++) {
      //lob(:,:,iii,jjj) = MATMUL( loby(:,:,iii,jjj), TRANSPOSE(uh2lob) )
      for(x=0; x<basisSize; x++) {
        for(v=0; v<numberOfVariables; v++) {
          lob[idx(z,y,x,v)] = 0;
          for(k=0; k<basisSize; k++) {
            lob[idx(z,y,x,v)] += tmpY[idx(z,y,k,v)] * uh2lob[basisSize-1][idxConv(k,x)];
          }
        }
      }
    }
  }
  
  delete[] tmpZ;
  delete[] tmpY;
  
  return lob;
}

//Fortran (Limiter.f90): GetSubcellData
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
            tmpZ[idxZ(z,y,x,v)] += luh[idxLuh(k,y,x,v)] * uh2lim[basisSize-1][idxConv(k,z)];
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
            tmpY[idxY(z,y,x,v)] += tmpZ[idxZ(z,k,x,v)] * uh2lim[basisSize-1][idxConv(k,y)];
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
            lim[idxLim(z,y,x,v)] += tmpY[idxY(z,y,k,v)] * uh2lim[basisSize-1][idxConv(k,x)];
          }
        }
      }
    }
  }
  
  delete[] tmpZ;
  delete[] tmpY;
  
  return lim;
}

void updateSubcellWithLimiterData(const double* const lim, const int numberOfVariables, const int basisSizeLim, const int basisSize, double* const luh) {
  
  idx4 idxLuh(basisSize, basisSize, basisSize, numberOfVariables);
  idx4 idxLim(basisSizeLim, basisSizeLim, basisSizeLim, numberOfVariables);
  idx2 idxConv(basisSizeLim, basisSize);
  
  double* tmpZ = new double[basisSize*basisSizeLim*basisSizeLim*numberOfVariables]; //Fortran ref: luhz(nVar,nSubLimV(1),nSubLimV(2),N+1)
  idx4 idxZ(basisSize, basisSizeLim, basisSizeLim, numberOfVariables);
  double* tmpY = new double[basisSize*basisSize*basisSizeLim*numberOfVariables]; //Fortran ref: luhy(nVar,nSubLimV(1),N+1,N+1)
  idx4 idxY(basisSize, basisSize, basisSizeLim, numberOfVariables);
  int x,y,z,v,k;
  
  for(y=0; y<basisSizeLim; y++) {
    for(x=0; x<basisSizeLim; x++) {
      //luhz(:,iii,jjj,:) = MATMUL( lim(:,iii,jjj,:), TRANSPOSE(lim2uh) ) 
      for(z=0; z<basisSize; z++) {
        for(v=0; v<numberOfVariables; v++) {
          tmpZ[idxZ(z,y,x,v)] = 0;
          for(k=0; k<basisSizeLim; k++) {
            tmpZ[idxZ(z,y,x,v)] += lim[idxLim(k,y,x,v)] * lim2uh[basisSize-1][idxConv(k,z)];
          }
          
        }
      }
    }
  }
  
  for(z=0; z<basisSize; z++) {
    for(x=0; x<basisSizeLim; x++) {
      //luhy(:,iii,:,jjj) = MATMUL( luhz(:,iii,:,jjj), TRANSPOSE(lim2uh) )
      for(y=0; y<basisSize; y++) {
        for(v=0; v<numberOfVariables; v++) {
          tmpY[idxY(z,y,x,v)] = 0;
          for(k=0; k<basisSizeLim; k++) {
            tmpY[idxY(z,y,x,v)] += tmpZ[idxZ(z,k,x,v)] * lim2uh[basisSize-1][idxConv(k,y)];
          }
        }
      }
    }
  }
  
  for(z=0; z<basisSize; z++) {
    for(y=0; y<basisSize; y++) {
      //luh(:,:,iii,jjj) = MATMUL( luhy(:,:,iii,jjj), TRANSPOSE(lim2uh) )  
      for(x=0; x<basisSize; x++) {
        for(v=0; v<numberOfVariables; v++) {
          luh[idxLuh(z,y,x,v)] = 0;
          for(k=0; k<basisSizeLim; k++) {
            luh[idxLuh(z,y,x,v)] += tmpY[idxY(z,y,k,v)] * lim2uh[basisSize-1][idxConv(k,x)];
          }
        }
      }
    }
  }
  
  delete[] tmpZ;
  delete[] tmpY;
  
} 


} // namespace c
} // namespace generic
} // namespace limiter
} // namespace kernel

#endif //DIMENSIONS == 3