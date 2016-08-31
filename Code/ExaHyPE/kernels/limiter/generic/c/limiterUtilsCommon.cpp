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

tarch::logging::Log _log("");

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

double* matrixInverse(int n, double* a) {
  
  double* ia = new double[n*n];
  double* c = new double[n*n*2];

  idx2 idx(n,n);
  idx2 idxC(n,2*n);
  
  int i,j,k,ml,mlV;
  double tmp, piv;
  
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      c[idxC(i,j)] = a[idx(j,i)];
    }
    for(; j<2*n; j++) {
      c[idxC(i,j)] = 0.;
    }
    c[idxC(i,i+n)] = 1.;
  }

  //Forward elimination and row swapping (if necessary)
  for(i=0; i<n; i++) {
    ml = i;
    mlV = abs(c[idxC(i,i)]);
    for(j=i+1; j<n; j++) {
      if(abs(c[idxC(j,i)]) > mlV) {
        ml = j;
        mlV = c[idxC(j,i)];
      }
    }
    for(k=0; k<2*n; k++) {
      tmp = c[idxC(ml,k)];
      c[idxC(ml,k)] = c[idxC(i,k)];
      c[idxC(i,k)] = tmp;
    }
    if(c[idxC(i,i)] == 0) {
      logError("matrixInverse()", "Matrix is singular" );
      return nullptr;
    }
    piv = 1. / c[idxC(i,i)];
    for(k=0; k<2*n; k++) {
      c[idxC(i,k)] *= piv;
    }
    for(j=i+1; j<n; j++) {
      tmp = c[idxC(j,i)];
      for(k=0; k<2*n; k++) {
        c[idxC(j,k)] -= tmp*c[idxC(i,k)];
      }
    }
  }
  
  //Back substitution
  for(i=n-1; i>=0; i--) {
    for(j=i-1; j>=0; j--) {
      tmp = c[idxC(j,i)];
      for(k=0; k<2*n; k++) {
        c[idxC(j,k)] -= tmp*c[idxC(i,k)];
      }
    }
  }
  
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      ia[idx(j,i)] = c[idxC(i,j+n)];
    }
  }
  
  return ia;
}


void initProjectionMatrices(const int basisSize) {
  int i,j,k;
  double* phi = new double[basisSize];
  
  {
    int basisSizeLim = getLimBasisSize(basisSize);
    idx2 idx(basisSize, basisSizeLim);
    uh2lim = new double[basisSize*basisSizeLim](); //initialized at 0, Fortran ref: uh2lim(nSubLim,N+1)
    
    const double dxi = 1. / basisSizeLim;
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
  
  {
    idx2 idx(basisSize, basisSize);
    uh2lob = new double[basisSize*basisSize];
    
    for(i=0; i<basisSize; i++) {
      BaseFunc1D(phi, gaussLobattoNodes[basisSize-1][i], basisSize);
      for(j=0; j<basisSize; j++) {
        uh2lob[idx(j,i)] = phi[j]; //Fortran: uh2lob(ii,:) = phi(:) 
      }
    }
  }
  
  {
    int basisSizeLim = getLimBasisSize(basisSize);
    idx2 idx(basisSizeLim, basisSize);
    lim2uh = new double[basisSizeLim*basisSize];
    
    double* lsqm = new double[(basisSize+1)*(basisSize+1)];
    double* lsqrhs = new double[basisSizeLim*(basisSize+1)];
    idx2 idxLSQM((basisSize+1),(basisSize+1));
    idx2 idxLSQrhs(basisSizeLim,(basisSize+1));
    idx2 idxUh2Lim(basisSize, basisSizeLim);
    const double dxi = 1.0 / basisSizeLim;
    for(i=0; i<basisSize; i++) {
      for(j=0; j<basisSize; j++) {
        lsqm[idxLSQM(i,j)] = 0.;
        for(k=0; k<basisSizeLim; k++) {
          lsqm[idxLSQM(i,j)] += 2* uh2lim[idxUh2Lim(i,k)] * uh2lim[idxUh2Lim(j,k)];
        }
      }
      lsqm[idxLSQM(i,basisSize)] = gaussLegendreWeights[basisSize-1][i];
    }
    for(i=0; i<basisSize; i++) {
      lsqm[idxLSQM(basisSize,i)] = -gaussLegendreWeights[basisSize-1][i];
    }
    lsqm[idxLSQM(basisSize,basisSize)] = 0.;
    
    double* ilsqm = matrixInverse(basisSize+1, lsqm);   
    
    for(i=0; i<basisSizeLim; i++) {
      for(j=0; j<basisSize; j++) {
        lsqrhs[idxLSQrhs(i,j)] = 2*uh2lim[idxUh2Lim(j,i)];
      }
      lsqrhs[idxLSQrhs(i,basisSize)] = dxi;
    }
    
    for(i=0; i<basisSizeLim; i++) {
      for(j=0; j<basisSize; j++) {
        lim2uh[idx(i,j)] = 0.;
        for(k=0; k<basisSize+1; k++) {
          lim2uh[idx(i,j)] += ilsqm[idxLSQM(k,j)] * lsqrhs[idxLSQrhs(i,k)];
        }
      }
    }
    
    delete[] lsqm;
    delete[] ilsqm;
    delete[] lsqrhs;
  }
  
  delete[] phi;
}

void freeProjectionMatrices(const int basisSize) {
  delete[] uh2lim;
  delete[] uh2lob;
  delete[] lim2uh;
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
