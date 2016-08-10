// ./make.sh
#include "sources/Kernels.h"
#include "sources/DGMatrices.h"
#include "sources/GaussLegendreQuadrature.h"
#include <iostream>
#include <fstream>
#include <set>
#include <string.h>
#include <iomanip>

#define order 3

using namespace std;
// order (see define)
// 2D
// nVar 5

int main(int argc, char* argv[]) {
  if(argc != 2) {
    cout << "path to input data missing" << endl;
    return 0;
  } 
  
  // set input stream
  ifstream source;
  source.open(argv[1], ios_base::in);  // open data
  if (!source)  {                     // if it does not work
    cerr << "Can't open Data!\n";
    return 0;
  }
  
  constexpr int pad = 3;
  constexpr int nVar = 5;
  constexpr int nVarPadded = nVar+pad;
  constexpr int nDOF =  order+1;
  constexpr int dim = 2;
  std::set<int> orders;

  int chunkSize = nDOF;

  double* lduh = new double[nVar * nDOF * nDOF];
  double lFbnd[4*5*6]; // chunkSize * nVar * 6
  
  const int startAddr_face1 = 0*4*5;
  const int startAddr_face2 = 1*4*5;
  const int startAddr_face3 = 2*4*5;
  const int startAddr_face4 = 3*4*5;
  const int startAddr_face5 = 4*4*5;
  const int startAddr_face6 = 5*4*5;
  
  
  // assume isotropic cell widths
  double* dx = new double[dim];
  dx[0] = 0.1; dx[1] = 0.1;
  
  // orders unused, for compatiblity with generic kernels
  kernels::aderdg::optimised::initDGMatrices(orders);
  kernels::aderdg::optimised::initGaussLegendreNodesAndWeights(orders);
  
  for(int i=0;i<nVar*nDOF*nDOF;i++) {
    source >> lduh[i];
  }
    
  for(int k=0;k<1;k++) {
    for(int j=0;j<nDOF;j++) {
      for(int iVar=0;iVar<nVar;iVar++) {
        source >> lFbnd[startAddr_face1+iVar*chunkSize+k*nDOF+j];
        source >> lFbnd[startAddr_face2+iVar*chunkSize+k*nDOF+j];
        source >> lFbnd[startAddr_face3+iVar*chunkSize+k*nDOF+j];
        source >> lFbnd[startAddr_face4+iVar*chunkSize+k*nDOF+j];
        source >> lFbnd[startAddr_face5+iVar*chunkSize+k*nDOF+j];
        source >> lFbnd[startAddr_face6+iVar*chunkSize+k*nDOF+j];
      }
    }
  }
 
  
  kernels::aderdg::optimised::surfaceIntegral(&lduh[0], &lFbnd[0], &dx[0]);

  for(int i=0;i<nVar*nDOF*nDOF;i++) {
    cout << setprecision (15) << lduh[i] << endl;
  }

  kernels::aderdg::optimised::freeGaussLegendreNodesAndWeights(orders);
  return 0;
}
