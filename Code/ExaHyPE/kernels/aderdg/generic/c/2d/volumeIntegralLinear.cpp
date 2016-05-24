#include "../../Kernels.h"

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

#if DIMENSIONS == 2

void volumeIntegralLinear(double* lduh, const double* const lFhi,
                          const tarch::la::Vector<DIMENSIONS, double>& dx,
                          const int numberOfVariables, const int basisSize) {
  // TODO(guera): Implement
}

#endif

}  // namespace c
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels
