// g++ -c -DDIMENSIONS=3 limiterUtils.cpp && g++ -c tmpTestLimiterUtils_Main.cpp && g++ -o main tmpTestLimiterUtils_Main.o ../../../../GaussLegendreQuadrature.o limiterUtils.o 

// 
#include <iostream>
#include <iomanip>
#include <set>

#include "../../Limiter.h"

void getTestData(double* &luh, int& numberOfVariable, int& basisSize) { 
  numberOfVariable = 5;
  basisSize = 4;
  int luhSize = numberOfVariable*basisSize*basisSize*basisSize;
  luh = new double[luhSize];
  for(int i=0; i<luhSize; i++) {
    luh[i] = i+ 0.5;
  }
}


/*
int main() {
  double* luh;
  int numberOfVariable, basisSize;
  getTestData(luh, numberOfVariable, basisSize);
  
  double* localMin = new double[numberOfVariable];
  double* localMax = new double[numberOfVariable];
  
  std::cout << "Hello World: " << luh[10] << std::endl;
  
  std::set<int> orders;
  kernels::initGaussLegendreNodesAndWeights(orders);
  std::cout << "gaussLegendreNodes: ";
  for(int i=0; i< basisSize; i++)
    std::cout << kernels::gaussLegendreNodes[basisSize-1][i] << " ";
  std::cout <<std::endl;
  
  //test BaseFunc1D
  double xi = 0.37;
  double* phi = new double[basisSize];
  kernels::limiter::generic::c::BaseFunc1D(phi, xi, basisSize);
  std::cout << "phi for xi = 0.37: ";
  for(int i=0; i< basisSize; i++)
    std::cout << phi[i] << " ";
  std::cout <<std::endl;
  
  kernels::limiter::generic::c::initProjectionMatrices(basisSize);
  int bLim = 2*(basisSize-1)+1;
  for(int i=0; i< basisSize; i++) {
    for(int j=0; j<bLim; j++)
      std::cout << kernels::limiter::generic::c::uh2lim[i*bLim+j] << " ";
    std::cout <<std::endl;
  }
  std::cout <<std::endl;
  std::cout << "uh2lim[1][1]: " <<kernels::limiter::generic::c::uh2lim[1*bLim+1] << std::endl;
  
  kernels::limiter::generic::c::findCellLocallocalMinlocalMax(luh, numberOfVariable, basisSize, localMin, localMax);
  
  std::cout << localMin[4] << std::endl;
  std::cout << localMax[4] << std::endl;
  
  double* troubledMin = new double[numberOfVariable];
  double* troubledMax = new double[numberOfVariable];
  
  double tMin = 0;
  double tMax = 319.499;
  std::cout << "troubled cell limit: " << tMin << " - " << tMax <<std::endl;
  
  for(int i=0;i<numberOfVariable;i++) {
    troubledMin[i] = tMin;
    troubledMax[i] = tMax;
  }
  
  if(kernels::limiter::generic::c::isTroubledCell(luh, numberOfVariable, basisSize, troubledMin, troubledMax))
  {
    std::cout << "troubled cell" << std::endl;
  } 
  else {
    std::cout << "not troubled cell" << std::endl;
  }
  
  int basisSizeLim = 0;
  double* lim = kernels::limiter::generic::c::getFVMData(luh, numberOfVariable, basisSize, basisSizeLim);
 
  std::cout << "basisSizeLim: " << basisSizeLim << std::endl;  
  kernels::idx4 idxLim(basisSizeLim,basisSizeLim,basisSizeLim,numberOfVariable);
  std::cout << "lim[iVar=0][x=0][y=1][z=3]: " << std::setprecision (15) << lim[idxLim(3,1,0,0)] << std::endl;
  
  return 0;
}
*/
