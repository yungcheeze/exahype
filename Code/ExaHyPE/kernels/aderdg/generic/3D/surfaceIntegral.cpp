#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

using std::endl;
using std::cout;

extern "C" {
void adersurfaceintegral_(double* lduh, double* lFhi, double* dx);
}

void kernels::aderdg::generic::fortran::surfaceIntegral(
    double* lduh, const double* const lFbnd,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    const int numberOfVariables, const int basisSize) {
  double* lFbndFortran =
      new double[numberOfVariables * 6 * basisSize * basisSize];
  for (int i = 0; i < numberOfVariables * 6 * basisSize * basisSize; i++) {
    lFbndFortran[i] = -123.45;
  }

  // Permutation of lFbnd and lFbnd

  // for(int i=0; i < numberOfVariables*6*basisSize*basisSize; i++){
  // lFbndFortran[i] = i;
  // }

  for (int bb = 0; bb < basisSize; bb++) {  // loop over dof
    for (int aa = 0; aa < basisSize; aa++) {
      for (int ivar = 0; ivar < numberOfVariables; ivar++) {
        for (int face = 0; face < 6; face++) {
          lFbndFortran[f2p4(ivar, face, aa, bb)] =
              lFbnd[p2f4(ivar, face, aa, bb)];
        }
      }
    }
  }

  // std::ofstream ofs;
  // ofs.open ("boutput_lFbnd.txt", std::ofstream::out);
  // for (int ii=0; ii<numberOfVariables*6*basisSize*basisSize; ii++) {
  // ofs << lFbnd[ii] << "\n";
  // }
  // ofs.close();

  // ofs.open ("boutput_lFbndFortran.txt", std::ofstream::out);
  // for (int ii=0; ii<numberOfVariables*6*basisSize*basisSize; ii++) {
  // ofs << lFbndFortran[ii] << "\n";
  // }
  // ofs.close();

  double* dxTemp = new double[3];
  dxTemp[0] = dx[0];
  dxTemp[1] = dx[1];
  dxTemp[2] = dx[2];

  adersurfaceintegral_(lduh, lFbndFortran, dxTemp);

  // for(int i=0; i < numberOfVariables*basisSize*basisSize*basisSize; i++){
  // cout << lduh[i] << endl;
  // }
  // exit(0);

  delete[] lFbndFortran;
  delete[] dxTemp;
}
