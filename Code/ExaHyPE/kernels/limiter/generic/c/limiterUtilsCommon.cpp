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

#include "../Limiter.h"

#include "tarch/la/ScalarOperations.h"

namespace kernels {
namespace limiter {
namespace generic {
namespace c {

/**
 * localMin, localMax are double[numberOfVariables]
 */
void findCellLocalMinAndMax(const double* const luh, const int numberOfVariables, const int basisSize, double* const localMin, double* const localMax) {
  int index, ii, iVar, iiEnd;
  
  // initialize and process luh
  index = 0;
  iiEnd =  basisSize*basisSize;
#if DIMENSIONS == 3
     iiEnd *= basisSize;
#endif
  for(int iVar = 0; iVar < numberOfVariables; iVar++) {
    localMin[iVar] = luh[index];
    localMax[iVar] = luh[index];   
    index++;
  }
  for(ii = 1; ii < iiEnd; ii++) {
    for(iVar = 0; iVar < numberOfVariables; iVar++) {
      localMin[iVar] = std::min ( localMin[iVar], luh[index] );
      localMax[iVar] = std::max ( localMax[iVar], luh[index] );
      index++;
    }
  }
 
  //process lob
  compareWithADERDGSolutionAtGaussLobattoNodes<false>(luh, numberOfVariables, basisSize, localMin, localMax, nullptr, 0.0);
}


/**
 * localMin, localMax are double[numberOfVariables]
 */
void findCellLocalMinAndMaxWithAnticipatedDG(const double* const luh, const double* const lduh, const double dt, const int numberOfVariables, const int basisSize, double* const localMin, double* const localMax) {
  int x, y, z, iVar; 
  const int order = basisSize-1;
  
  // initialize and process luh
#if DIMENSIONS == 3
  const int basisSize3D = basisSize;
#else
  const int basisSize3D = 1;  
#endif

  idx4 idx(basisSize3D, basisSize, basisSize, numberOfVariables);

  for(int iVar = 0; iVar < numberOfVariables; iVar++) {
    localMin[iVar] = anticipateLuh(luh, lduh, dt, order, iVar, 0, 0, 0);
    localMax[iVar] = anticipateLuh(luh, lduh, dt, order, iVar, 0, 0, 0);   
  }
  for(z=0; z<basisSize3D; z++) {
    for(y=0; y<basisSize; y++) {
      for(x=0; x<basisSize; x++) {
        for(iVar = 0; iVar < numberOfVariables; iVar++) {
          localMin[iVar] = std::min ( localMin[iVar], anticipateLuh(luh, lduh, dt, order, idx(z, y, x, iVar), x, y, z) );
          localMax[iVar] = std::max ( localMax[iVar], anticipateLuh(luh, lduh, dt, order, idx(z, y, x, iVar), x, y, z) );
        }
      }
    }
  }
 
  //process lob
  compareWithADERDGSolutionAtGaussLobattoNodes<true>(luh, numberOfVariables, basisSize, localMin, localMax, lduh, dt);
}

double getMin(const double* const solutionMin, int iVar, int numberOfVariables) {
  double result = std::numeric_limits<double>::max();
  for (int i=0; i<DIMENSIONS_TIMES_TWO; i+=numberOfVariables) {
    result = std::min( result, solutionMin[i+iVar] );
  }
  return result;
}
double getMax(const double* const solutionMax, int iVar, int numberOfVariables) {
  double result = -std::numeric_limits<double>::max();
  for (int i=0; i<DIMENSIONS_TIMES_TWO; i+=numberOfVariables) {
    result = std::max( result, solutionMax[i+iVar] );
  }
  return result;
}

bool isTroubledCell(const double* const luh, const double* const lduh, const double dt, const int numberOfVariables, const int basisSize, const double* const troubledMin, const double* const troubledMax) {
  double minMarginOfError = 0.0001;
  double diffScaling      = 0.001;

  double* localMin = new double[numberOfVariables];
  double* localMax = new double[numberOfVariables];

  findCellLocalMinAndMaxWithAnticipatedDG(luh, lduh, dt, numberOfVariables, basisSize, localMin, localMax);

  for(int iVar = 0; iVar < numberOfVariables; iVar++) {
    double previousSolutionMax = getMax(troubledMax,iVar,numberOfVariables);
    double previousSolutionMin = getMin(troubledMin,iVar,numberOfVariables);

    double ldiff = (previousSolutionMax - previousSolutionMin) * diffScaling;
    assertion1(tarch::la::greaterEquals(ldiff,0.0),ldiff);
    ldiff = std::max( ldiff, minMarginOfError );

    if((localMin[iVar] < (previousSolutionMin - ldiff)) ||
       (localMax[iVar] > (previousSolutionMax + ldiff))) {
      delete[] localMin;
      delete[] localMax;
      return true;
    }
  }
  
  delete[] localMin;
  delete[] localMax;
  
  //TODO JMG (todo or not needed???) check PDE positivity and lim data not NAN 
  
  return false;
}

//Used only for test purpose to check the projection since expandBoundariesWithGaussLobatto uses the same loop
double* getGaussLobattoData(const double* const luh, const int numberOfVariables, const int basisSize) {

#if DIMENSIONS == 3
  const int basisSize3D = basisSize;
#else
  const int basisSize3D = 1;  
#endif

  idx4 idx(basisSize3D, basisSize, basisSize, numberOfVariables);
  idx2 idxConv(basisSize, basisSize);
  
  double* lob = new double[basisSize3D*basisSize*basisSize*numberOfVariables]; //Fortran ref: lob(nVar,nDOF(1),nDOF(2),nDOF(3))
  
  int x,y,z,v,ix,iy,iz;
  
  for(z=0; z<basisSize3D; z++) {
    for(y=0; y<basisSize; y++) {
      for(x=0; x<basisSize; x++) {
        for(v=0; v<numberOfVariables; v++) {
          lob[idx(z,y,x,v)] = 0;
          for(iz=0; iz<basisSize3D; iz++) {
            for(iy=0; iy<basisSize; iy++) {
              for(ix=0; ix<basisSize;ix++) {
                lob[idx(z,y,x,v)] += luh[idx(iz,iy,ix,v)] 
#if DIMENSIONS == 3
                                        * uh2lob[basisSize-1][idxConv(iz,z)] 
#endif
                                        * uh2lob[basisSize-1][idxConv(iy,y)] * uh2lob[basisSize-1][idxConv(ix,x)];
              }
            }
          }
        }
      }
    }
  }
  
  return lob;
}

/**
 * Auxilliary function to findMinMax
 * Project to GaussLobatto and modify the min/max if required
 *
 */
template <bool anticipateNewDGSolution>
void compareWithADERDGSolutionAtGaussLobattoNodes(const double* const luh, const int numberOfVariables, const int basisSize, double* const min, double* const max, const double* const lduh, const double dt) {

#if DIMENSIONS == 3
  const int basisSize3D = basisSize;
#else
  const int basisSize3D = 1;  
#endif
const int order = basisSize-1;

  idx4 idx(basisSize3D, basisSize, basisSize, numberOfVariables);
  idx2 idxConv(basisSize, basisSize);
  
  double lobValue;
  int x,y,z,v,ix,iy,iz;
  
  for(z=0; z<basisSize3D; z++) {
    for(y=0; y<basisSize; y++) {
      for(x=0; x<basisSize; x++) {
        for(v=0; v<numberOfVariables; v++) {
          lobValue = 0;
          for(iz=0; iz<basisSize3D; iz++) {
            for(iy=0; iy<basisSize; iy++) {
              for(ix=0; ix<basisSize;ix++) {
                if(anticipateNewDGSolution) {
                  lobValue += anticipateLuh(luh, lduh, dt, order, idx(iz,iy,ix,v), ix, iy, iz)
#if DIMENSIONS == 3
                                        * uh2lob[order][idxConv(iz,z)] 
#endif
                                        * uh2lob[order][idxConv(iy,y)] * uh2lob[order][idxConv(ix,x)];
                } else {
                  lobValue += luh[idx(iz,iy,ix,v)] 
#if DIMENSIONS == 3
                                        * uh2lob[order][idxConv(iz,z)] 
#endif
                                        * uh2lob[order][idxConv(iy,y)] * uh2lob[order][idxConv(ix,x)];
                }
              }
            }
          }
          if(min[v] > lobValue) {
            min[v] = lobValue;
          } else if(max[v] < lobValue) {
            max[v] = lobValue;
          }
        }
      }
    }
  }

}

//Fortran (Limiter.f90): GetSubcellData
// Allocate lim memory via
// double* lim = new double[basisSizeLim*basisSizeLim*basisSizeLim*numberOfVariables]; //Fortran ref: lim(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3))
void projectOnFVLimiterSpace(const double* const luh, const int numberOfVariables, const int basisSize, const int basisSizeLim, double* const lim) {
  
#if DIMENSIONS == 3
  const int basisSize3D = basisSize;
  const int basisSizeLim3D = basisSizeLim;
#else
  const int basisSize3D = 1;  
  const int basisSizeLim3D = 1;
#endif

  idx4 idxLuh(basisSize3D, basisSize, basisSize, numberOfVariables);
  idx4 idxLim(basisSizeLim3D, basisSizeLim, basisSizeLim, numberOfVariables);
  idx2 idxConv(basisSize, basisSizeLim); 
  
  int x,y,z,v,ix,iy,iz; 
  
  for(z=0; z<basisSizeLim3D; z++) {
    for(y=0; y<basisSizeLim; y++) {
      //lim(:,:,iii,jjj) = MATMUL( limy(:,:,iii,jjj), TRANSPOSE(uh2lim) )
      for(x=0; x<basisSizeLim; x++) {
        for(v=0; v<numberOfVariables; v++) {
          lim[idxLim(z,y,x,v)] = 0;
          for(iz=0; iz<basisSize3D; iz++) {
            for(iy=0; iy<basisSize; iy++) {
              for(ix=0; ix<basisSize; ix++) {
                lim[idxLim(z,y,x,v)] += luh[idxLuh(iz,iy,ix,v)] 
#if DIMENSIONS == 3
                                          * uh2lim[basisSize-1][idxConv(iz,z)]
#endif
                                          * uh2lim[basisSize-1][idxConv(iy,y)] * uh2lim[basisSize-1][idxConv(ix,x)];
              }
            }
          }
        }
      }
    }
  }
  
}

void projectOnADERDGSpace(const double* const lim, const int numberOfVariables, const int basisSizeLim, const int basisSize, double* const luh) {
  
#if DIMENSIONS == 3
  const int basisSize3D = basisSize;
  const int basisSizeLim3D = basisSizeLim;
#else
  const int basisSize3D = 1;  
  const int basisSizeLim3D = 1;
#endif
  
  idx4 idxLuh(basisSize3D, basisSize, basisSize, numberOfVariables);
  idx4 idxLim(basisSizeLim3D, basisSizeLim, basisSizeLim, numberOfVariables);
  idx2 idxConv(basisSizeLim, basisSize);
  
  int x,y,z,v,ix,iy,iz;

  for(z=0; z<basisSize3D; z++) {
    for(y=0; y<basisSize; y++) {
      //luh(:,:,iii,jjj) = MATMUL( luhy(:,:,iii,jjj), TRANSPOSE(lim2uh) )  
      for(x=0; x<basisSize; x++) {
        for(v=0; v<numberOfVariables; v++) {
          luh[idxLuh(z,y,x,v)] = 0;
          for(iz=0; iz<basisSizeLim3D; iz++) {
            for(iy=0; iy<basisSizeLim; iy++) {
              for(ix=0; ix<basisSizeLim; ix++) {
                luh[idxLuh(z,y,x,v)] += lim[idxLim(iz,iy,ix,v)] 
#if DIMENSIONS == 3
                                          * lim2uh[basisSize-1][idxConv(iz,z)] 
#endif
                                          * lim2uh[basisSize-1][idxConv(iy,y)] * lim2uh[basisSize-1][idxConv(ix,x)];
              }
            }
          }
        }
      }
    }
  }
  
} 

} // namespace c
} // namespace generic
} // namespace limiter
} // namespace kernel
