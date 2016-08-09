// ./make.sh
#include "Kernels.h"
#include "DGMatrices.h"
#include "GaussLegendreQuadrature.h"
#include <iostream>
#include <fstream>
#include <set>
#include <string.h>
#include <iomanip>

using namespace std;
// order 3 (nDOF 4)
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
  constexpr int nDOF = 4;
  constexpr int dim = 2;
  std::set<int> orders;
  
  cout.precision(6);

  double* lduh = new double[nVar * nDOF * nDOF]; // not initialised
  double *lFhi = new double[nVarPadded * nDOF* nDOF * dim]; 
  double *lFhi_x = &lFhi[0];
  double *lFhi_y = &lFhi[nVarPadded*nDOF*nDOF]; // (nVar + Pad) * nDOF**2
  
  memset(lFhi, 0.0, 256*sizeof(double)); 
  
  // assume isotropic cell widths
  double* dx = new double[dim];
  dx[0] = 0.1; dx[1] = 0.1;
  
  // orders unused, for compatiblity with generic kernels
  kernels::aderdg::optimised::initDGMatrices(orders);
  kernels::aderdg::optimised::initGaussLegendreNodesAndWeights(orders);
  
  // Get input data
  for(int iVar=0;iVar<nVar;iVar++) {
    for(int i=0;i<nDOF;i++) {
      for(int j=0;j<nDOF;j++) {
        source >> lFhi_x[iVar + i*nVarPadded + j*nVarPadded*nDOF];      
      }
    }
  }
  for(int iVar=0;iVar<nVar;iVar++) {
    for(int j=0;j<nDOF;j++) {
      for(int i=0;i<nDOF;i++) {
        source >> lFhi_y[iVar + i*nVarPadded + j*nVarPadded*nDOF];        
      }
    }
  }
  
/*  
  double cnt = 1.0;
  
  for(int j=0;j<nDOF;j++) {
    for(int i=0;i<nDOF;i++) {
      for(int iVar=0;iVar<nVar;iVar++) {
        lFhi_x[iVar + i*nVarPadded + j*nVarPadded*nDOF] = cnt;
        cnt++;         
      }
    }
  }
  for(int i=0;i<nDOF;i++) {
    for(int j=0;j<nDOF;j++) {
      for(int iVar=0;iVar<nVar;iVar++) {
        lFhi_y[iVar + j*nVarPadded + i*nVarPadded*nDOF] = cnt;
        cnt++;         
      }
    }
  }
 */ 
 /*
  for(int i=0;i<256;i++) {
    cout << lFhi[i] << endl;
  }
  */
  //cout << endl << lFhi_y[2*8*4+3*8+2] << endl;
  
  kernels::aderdg::optimised::volumeIntegral( 
    lduh, 
    lFhi, 
    dx
  );

  for(int i=0;i<80;i++) {
    cout << setprecision (15) << lduh[i] << endl;
  }

  kernels::aderdg::optimised::freeGaussLegendreNodesAndWeights(orders);
  return 0;
}
