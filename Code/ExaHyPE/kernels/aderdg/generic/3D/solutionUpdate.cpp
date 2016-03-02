#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

using std::endl;
using std::cout;

extern "C" 
{
  void elementupdate_(double *luh, double *lduh, double *dt);
}
void kernels::aderdg::generic::solutionUpdate(double * luh,
                                              const double * const lduh,
                                              const double dt,
                                              const int numberOfVariables,
                                              const int basisSize) {
const int order = basisSize-1;

double* lduhFortran = new double[numberOfVariables*basisSize*basisSize*basisSize];
for(int i=0; i < numberOfVariables*basisSize*basisSize*basisSize; i++){
  lduhFortran[i] = lduh[i];
}

double* dtTemp = new double[1];
dtTemp[0] = dt;

elementupdate_(luh, lduhFortran, dtTemp);

delete[] lduhFortran;
delete dtTemp;

}


