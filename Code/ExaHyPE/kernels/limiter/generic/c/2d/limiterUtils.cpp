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
 
#include "peano/utils/Loop.h"

#if DIMENSIONS == 2
 
namespace kernels {
namespace limiter {
namespace generic {
namespace c {

<<<<<<< HEAD
//Fortran (Limiter.f90): GetLobattoData
double* getGaussLobattoData(const double* const luh, const int numberOfVariables, const int basisSize) {
  
  idx3 idx(basisSize, basisSize, numberOfVariables);
  idx2 idxConv(basisSize, basisSize);
  
  double* lob = new double[basisSize*basisSize*numberOfVariables]; //Fortran ref: lob(nVar,nDOF(1),nDOF(2),nDOF(3))
  
  double* tmpY = new double[basisSize*basisSize*numberOfVariables]; //Fortran ref: loby(nVar,nDOF(1),nDOF(2),nDOF(3))
  int x,y,v,k;
  
  // TODO(JM): Compute
  //
  //    lob_yxv = \sum_j \sum_i luh_jiv * uh2lob_jy * uh2lob_ix ;
  //
  // It's simply a projection with respect to two directions.
  // Then compare directly with localMin and localMax (add them to the argument list).
  // No need for storing things.

  for(x=0; x<basisSize; x++) {
    //loby(:,iii,:,jjj) = MATMUL( lobz(:,iii,:,jjj), TRANSPOSE(uh2lob) ) 
    for(y=0; y<basisSize; y++) {
      for(v=0; v<numberOfVariables; v++) {
        tmpY[idx(y,x,v)] = 0;
        for(k=0; k<basisSize; k++) {
          tmpY[idx(y,x,v)] += luh[idx(k,x,v)] * uh2lob[basisSize-1][idxConv(k,y)];
        }
      }
    }
  }

  for(y=0; y<basisSize; y++) {
    //lob(:,:,iii,jjj) = MATMUL( loby(:,:,iii,jjj), TRANSPOSE(uh2lob) )
    for(x=0; x<basisSize; x++) {
      for(v=0; v<numberOfVariables; v++) {
        lob[idx(y,x,v)] = 0;
        for(k=0; k<basisSize; k++) {
          lob[idx(y,x,v)] += tmpY[idx(y,k,v)] * uh2lob[basisSize-1][idxConv(k,x)];
        }
      }
    }
  }
  
  delete[] tmpY;
  
  return lob;
}

=======
/*
>>>>>>> Limiter Kernel - Refactor the projections to remove tmp arrays and make them dimension agnostic
//Fortran (Limiter.f90): GetSubcellData
// Allocate lim via
// double* lim = new double[basisSizeLim*basisSizeLim*numberOfVariables]; //Fortran ref: lim(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3))
void getFVMData(const double* const luh, const int numberOfVariables, const int basisSize, double* lim) {
  const int basisSizeLim = getLimBasisSize(basisSize);
  
  idx3 idxLuh(basisSize, basisSize, numberOfVariables);
  idx2 idxConv(basisSize, basisSizeLim);

  idx3 idxLim(basisSizeLim, basisSizeLim, numberOfVariables);
  
  double* tmpY = new double[basisSizeLim*basisSize*numberOfVariables]; //Fortran ref: limy(nVar,nDOF(1),nSubLimV(2),nSubLimV(3))
  idx3 idxY(basisSizeLim, basisSize, numberOfVariables);
  int x,y,v,k;
  
  // TODO(JM): Compute
  //
  //    lim_yxv = \sum_j \sum_i luh_jiv * uh2lim_jy * uh2lim_ix ;
  //
  // It's simply a projection with respect to two directions.
  // Then compare directly with localMin and localMax (add them to the argument list).
  // No need for storing things.
  // This looks pretty similar to the other projections,
  // doesn't it?
  //
  // If you use peano's dfor loop, you will get it done
  // in a dimension independent fashion.

  for(x=0; x<basisSize; x++) {
    //limy(:,iii,:,jjj) = MATMUL( limz(:,iii,:,jjj), TRANSPOSE(uh2lim) ) 
    for(y=0; y<basisSizeLim; y++) {
      for(v=0; v<numberOfVariables; v++) {
        tmpY[idxY(y,x,v)] = 0;
        for(k=0; k<basisSize; k++) {
          tmpY[idxY(y,x,v)] += luh[idxLuh(k,x,v)] * uh2lim[basisSize-1][idxConv(k,y)];
        }
      }
    }
  }

  for(y=0; y<basisSizeLim; y++) {
    //lim(:,:,iii,jjj) = MATMUL( limy(:,:,iii,jjj), TRANSPOSE(uh2lim) )
    for(x=0; x<basisSizeLim; x++) {
      for(v=0; v<numberOfVariables; v++) {
        lim[idxLim(y,x,v)] = 0;
        for(k=0; k<basisSize; k++) {
          lim[idxLim(y,x,v)] += tmpY[idxY(y,k,v)] * uh2lim[basisSize-1][idxConv(k,x)];
        }
      }
    }
  }

  delete[] tmpY;
}


void updateSubcellWithLimiterData(const double* const lim, const int numberOfVariables, const int basisSizeLim, const int basisSize, double* const luh) {
  
  idx3 idxLuh(basisSize, basisSize, numberOfVariables);
  idx3 idxLim(basisSizeLim, basisSizeLim, numberOfVariables);
  idx2 idxConv(basisSizeLim, basisSize);
  
  const int basisSize2 = basisSize*basisSize;

  double* tmpY = new double[basisSize2*basisSizeLim*numberOfVariables]; //Fortran ref: luhy(nVar,nSubLimV(1),N+1,N+1)
  idx3 idxY(basisSize, basisSizeLim, numberOfVariables);
  int x,y,v,k;
  

  for(x=0; x<basisSizeLim; x++) {
    //luhy(:,iii,:,jjj) = MATMUL( luhz(:,iii,:,jjj), TRANSPOSE(lim2uh) )
    for(y=0; y<basisSize; y++) {
      for(v=0; v<numberOfVariables; v++) {
        tmpY[idxY(y,x,v)] = 0;
        for(k=0; k<basisSizeLim; k++) {
          tmpY[idxY(y,x,v)] += lim[idxLim(k,x,v)] * lim2uh[basisSize-1][idxConv(k,y)];
        }
      }
    }
  }

  for(y=0; y<basisSize; y++) {
    //luh(:,:,iii,jjj) = MATMUL( luhy(:,:,iii,jjj), TRANSPOSE(lim2uh) )  
    for(x=0; x<basisSize; x++) {
      for(v=0; v<numberOfVariables; v++) {
        luh[idxLuh(y,x,v)] = 0;
        for(k=0; k<basisSizeLim; k++) {
          luh[idxLuh(y,x,v)] += tmpY[idxY(y,k,v)] * lim2uh[basisSize-1][idxConv(k,x)];
        }
      }
    }
  }
  
  delete[] tmpY;
  
} 
*/

} // namespace c
} // namespace generic
} // namespace limiter
} // namespace kernel

#endif //DIMENSIONS == 2
