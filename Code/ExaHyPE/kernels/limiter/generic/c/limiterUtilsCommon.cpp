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

namespace kernels {
namespace limiter {
namespace generic {
namespace c {

// TODO JMG get true basisSize from spec file
int getLimBasisSize(const int basisSize) {
  return 2*(basisSize-1) +1;
}

/**
 * localMin, localMax are double[numberOfVariables]
 */
void findCellLocalMinAndMax(const double* const luh, const double* const lim, const int numberOfVariables, const int basisSize, double* localMin, double* localMax) {
  int index, ii, iVar, iiEnd;
  
  // initialize and process luh
  index = 0;
  iiEnd =  basisSize*basisSize;
  if(DIMENSIONS == 3)
     iiEnd *= basisSize;
  for(int iVar = 0; iVar < numberOfVariables; iVar++) {
    localMin[iVar] = luh[index];
    localMax[iVar] = luh[index];   
    index++;
  }
  for(ii = 1; ii < iiEnd; ii++) {
    for(iVar = 0; iVar < numberOfVariables; iVar++) {
      if(luh[index] < localMin[iVar]) {
        localMin[iVar] = luh[index];
      } else if(luh[index] > localMax[iVar]) {
        localMax[iVar] = luh[index];
      }    
      index++;
    }
  }
 
  // process lob
  double* lob = getGaussLobattoData(luh, numberOfVariables, basisSize);
  index = 0;
  iiEnd =  basisSize*basisSize;
  if(DIMENSIONS == 3)
     iiEnd *= basisSize;
  for(ii = 0; ii < iiEnd; ii++) {
    for(iVar = 0; iVar < numberOfVariables; iVar++) {
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
  const int basisSizeLim = getLimBasisSize(basisSize);
  index = 0;
  iiEnd =  basisSizeLim*basisSizeLim;
  if(DIMENSIONS == 3)
     iiEnd *= basisSizeLim;
  for(ii = 0; ii < iiEnd; ii++) {
    for(iVar = 0; iVar < numberOfVariables; iVar++) {
      if(lim[index] < localMin[iVar]) {
        localMin[iVar] = lim[index];
      } else if(lim[index] > localMax[iVar]) {
        localMax[iVar] = lim[index];
      }    
      index++;
    }
  }
}

double getMin(const double* const minOfNeighbours, int iVar, int numberOfVariables) {
  double min = std::numeric_limits<double>::max();
  for (int i=0; i<DIMENSIONS_TIMES_TWO; i+=numberOfVariables) {
    min = std::min( min, minOfNeighbours[i+iVar] );
  }
  return min;
}
double getMax(const double* const maxOfNeighbours, int iVar, int numberOfVariables) {
  double max = std::numeric_limits<double>::min();
  for (int i=0; i<DIMENSIONS_TIMES_TWO; i+=numberOfVariables) {
    max = std::max( max, maxOfNeighbours[i+iVar] );
  }
  return max;
}

bool isTroubledCell(const double* const luh, const int numberOfVariables, const int basisSize, const double* const troubledMin, const double* const troubledMax) {
  
  double minMarginOfError = 0.0001;
  double diffScaling      = 0.001;
  
  double* localMin = new double[numberOfVariables];
  double* localMax = new double[numberOfVariables];

  const int basisSizeLim = getLimBasisSize(basisSize);
  double* lim = new double[basisSizeLim*basisSizeLim*numberOfVariables]; //Fortran ref: lim(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3))
  getFVMData(luh, numberOfVariables, basisSize, lim);
  findCellLocalMinAndMax(luh, lim, numberOfVariables, basisSize, localMin, localMax);
  
  double ldiff;

  for(int iVar = 0; iVar < numberOfVariables; iVar++) {
    double maxMaxOfNeighbours = getMin(troubledMax,iVar,numberOfVariables);
    double minMinOfNeighbours = getMax(troubledMin,iVar,numberOfVariables);

    ldiff = std::max((maxMaxOfNeighbours - minMinOfNeighbours) * diffScaling, minMarginOfError);
    if((localMin[iVar] < (minMinOfNeighbours - ldiff)) ||
       (localMax[iVar] > (maxMaxOfNeighbours + ldiff))) {
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
