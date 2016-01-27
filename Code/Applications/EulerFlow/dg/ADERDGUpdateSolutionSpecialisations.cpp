#include "EulerFlow/dg/ADERDG.h"
#include "EulerFlow/quad/GaussLegendre.h"


// 3D specialisation
template <>
void exahype::dg::updateSolution<3>(
    double * restrict /*inout*/ luh,
    const double * restrict const /*in*/ lduh,
    const double /*in*/ dt
) {
  // todo insert your code here
}

// 2D specialisation
template <>
void exahype::dg::updateSolution<2>(
    double * restrict /*inout*/ luh,
    const double * restrict const /*in*/ lduh,
    const double /*in*/ dt
) {
  constexpr int nvar        = EXAHYPE_NVARS;
  constexpr int basisSize   = EXAHYPE_ORDER+1;

  for (int ii=0; ii<basisSize; ii++) {
    for (int jj=0; jj<basisSize; jj++) {
      const int nodeIndex     = jj + basisSize * ii;
      const int dofStartIndex = nodeIndex * nvar;

      double weight =  exahype::quad::gaussLegendreWeights[ii] * exahype::quad::gaussLegendreWeights[jj];
      const double updateSize = dt/weight;

      for(int ivar=0; ivar < nvar; ivar++) {
        luh[dofStartIndex+ivar] +=  lduh[dofStartIndex+ivar]*updateSize;
      }
    }
  }
}
