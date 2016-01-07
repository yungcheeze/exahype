#include "EulerFlow/Constants.h"

namespace exahype {
  namespace dg {
    /**
     * \note Unused function. The flux matrices are precomputed and hardcoded in DGMatrices.
     */
    void BaseFunc(double phi[EXAHYPE_ORDER+1], double phi_xi[EXAHYPE_ORDER+1], double xi, double xin[EXAHYPE_ORDER+1]);
  }  // namespace dg
}  // namespace exahype
