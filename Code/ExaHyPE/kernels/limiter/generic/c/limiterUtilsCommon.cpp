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
  
double* uh2lim;
double* uh2lob;
double* lim2uh;

// TODO JMG get true basisSize from spec file
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
  iiEnd =  basisSizeLob*basisSizeLob;
  if(DIMENSIONS == 3)
     iiEnd *= basisSizeLob;
  for(ii = 0; ii < iiEnd; ii++) {
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
  
  int basisSizeLim = 0;
  double* lim = getFVMData(luh, numberOfVariables, basisSize, basisSizeLim);
  index = 0;
  iiEnd =  basisSizeLim*basisSizeLim;
  if(DIMENSIONS == 3)
     iiEnd *= basisSizeLim;
  for(ii = 0; ii < iiEnd; ii++) {
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
