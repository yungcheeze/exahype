#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

#include <fstream> 

using std::endl;
using std::cout;

extern "C" 
{
  void adervolumeintegrallinear_(double *lduh, double *lFhi, double *dx);
}

void kernels::aderdg::generic::fortran::volumeIntegralLinear(
    double * lduh,
    const double * const lFhi,
    const tarch::la::Vector<DIMENSIONS,double>& dx,
    const int numberOfVariables,
    const int basisSize
) {
  // todo Angelika
  // Please remove the typedefs in generic kernels again since numberOf(...)Dof is not
  // a compile time variable anymore

 
  double* lFhiFortran = new double[numberOfVariables*DIMENSIONS*basisSize*basisSize*basisSize];
  for(int i=0; i < numberOfVariables*DIMENSIONS*basisSize*basisSize*basisSize; i++)
    lFhiFortran[i] = -123.45;

  for (int ii=0; ii<basisSize; ii++) {  // loop over dof
    for (int jj=0; jj<basisSize; jj++) {
      for (int kk=0; kk<basisSize; kk++) {
        for(int ivar=0; ivar < numberOfVariables; ivar++) {
          for(int dim=0; dim < DIMENSIONS; dim++) {
            lFhiFortran[f2p5(ivar, dim, ii, jj, kk)] = lFhi[p2f5(ivar, dim, ii, jj, kk)];
          }
        }
      }
    }
  }

  
  // std::ofstream ofs;
  // ofs.open ("boutput_lFhi.txt", std::ofstream::out);
  // for (int ii=0; ii<numberOfVariables*DIMENSIONS*basisSize*basisSize*basisSize; ii++) {
    // ofs << lFhi[ii] << "\n";
  // }
  // ofs.close();

  // ofs.open ("boutput_lFhiFortran.txt", std::ofstream::out);
  // for (int ii=0; ii<numberOfVariables*DIMENSIONS*basisSize*basisSize*basisSize; ii++) {
    // ofs << lFhiFortran[ii] << "\n";
  // }
  // ofs.close();

  // cout << "-------------lFhi in volumeIntegral.cpph------------------" << "\n";
  // cout << lFhi[0] << "\n";
  // cout << lFhi[1]<< "\n";
  // cout << lFhi[2] << "\n";
  // cout << lFhi[3] << "\n";
  // cout << lFhi[4] << "\n";
  // cout << lFhi[5] << "\n";
  // cout << "-------------lFhi in volumeIntegral.cpph------------------" << "\n";
  
  
  double* dxTemp = new double[3];
  dxTemp[0]= dx[0];
  dxTemp[1]= dx[1];
  dxTemp[2]= dx[2];
  
  adervolumeintegrallinear_(lduh, lFhiFortran, dxTemp);

  delete[] lFhiFortran;
  delete[] dxTemp;
}

