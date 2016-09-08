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
void findCellLocallocalMinlocalMax(const double* const luh, const int numberOfVariables, const int basisSize, double* localMin, double* localMax) {      

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
  int basisSizeLim = 0;
  double* lim = getFVMData(luh, numberOfVariables, basisSize, basisSizeLim);
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
