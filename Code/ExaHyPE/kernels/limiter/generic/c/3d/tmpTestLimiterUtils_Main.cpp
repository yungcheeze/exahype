#include <iostream>

#define DIMENSIONS 3 
#include "limiterUtils.cpph"

void getTestData(double* &luh, int& numberOfVariable, int& basisSize) { 
  numberOfVariable = 5;
  basisSize = 4;
  int luhSize = numberOfVariable*basisSize*basisSize*basisSize;
  luh = new double[luhSize];
  for(int i=0; i<luhSize; i++) {
    luh[i] = i+ 0.5;
  }
}


/**
int main() {
  double* luh;
  int numberOfVariable, basisSize;
  getTestData(luh, numberOfVariable, basisSize);
  
  double* localMin = new double[numberOfVariable];
  double* localMax = new double[numberOfVariable];
  
  std::cout << "Hello World: " << luh[10] << std::endl;
  
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
  
  return 0;
}
*/
