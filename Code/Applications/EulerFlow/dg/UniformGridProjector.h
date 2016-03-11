#ifndef EXAHYPE_DG_BASISPROJECTOR_H_
#define EXAHYPE_DG_BASISPROJECTOR_H_

namespace exahype {
namespace dg {
template <int dim, int order, int basisSizePowerDim>
struct UniformGridProjector;
}  // namespace dg
}  // namespace exahype

template <int dim, int order, int basisSizePowerDim>
class exahype::dg::UniformGridProjector {
 public:
  const static double matrix[basisSizePowerDim][basisSizePowerDim];
};

// explicit specialisations
#include "EulerFlow/dg/UniformGridProjector23Specialisation.h"

#endif /* EXAHYPE_DG_BASISPROJECTOR_H_ */
