#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"


void kernels::aderdg::generic::solutionUpdate(double * luh,
                                              const double * const lduh,
                                              const double dt,
                                              const int numberOfVariables,
                                              const int basisSize) {
  const int order = basisSize-1;

  for (int ii=0; ii<basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {
      const int nodeIndex     = jj + basisSize * ii;
      const int dofStartIndex = nodeIndex * numberOfVariables;

      const double weight     =  kernels::gaussLegendreWeights[order][ii] * kernels::gaussLegendreWeights[order][jj];
      const double updateSize = dt/weight;

      for(int ivar=0; ivar < numberOfVariables; ivar++) {
        luh[dofStartIndex+ivar] +=  lduh[dofStartIndex+ivar]*updateSize;
      }
    }
  }
}


