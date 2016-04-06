#include "kernels/aderdg/generic/Kernels.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

#if DIMENSIONS == 2
void kernels::aderdg::generic::c::faceUnknownsProlongation(double* lQhbndFine,
                              double* lFhbndFine,
                              const double* lQhbndCoarse,
                              const double* lFhbndCoarse,
                              const int coarseGridLevel,
                              const int fineGridLevel,
                              const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex,
                              const int numberOfVariables, const int basisSize){
  // @todo Please implement!
}

void kernels::aderdg::generic::c::faceUnknownsRestriction(double* lQhbndCoarse,
                             double* lFhbndCoarse,
                             const double* lQhbndFine,
                             const double* lFhbndFine,
                             const int coarseGridLevel,
                             const int fineGridLevel,
                             const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex,
                             const int numberOfVariables, const int basisSize){
  // @todo Please implement!
}

void kernels::aderdg::generic::c::volumeUnknownsProlongation(double* luhFine,
                                const double* luhCoarse,
                                const int coarseGridLevel,
                                const int fineGridLevel,
                                const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
                                const int numberOfVariables, const int basisSize){
  // @todo Please implement!
}

void kernels::aderdg::generic::c::volumeUnknownsRestriction(double* luhCoarse,
                               const double* luhFine,
                               const int coarseGridLevel,
                               const int fineGridLevel,
                               const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
                               const int numberOfVariables, const int basisSize){
  // @todo Please implement!
}
#endif
