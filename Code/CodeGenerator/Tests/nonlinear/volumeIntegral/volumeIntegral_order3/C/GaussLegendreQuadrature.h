#ifndef _EXAHYPE_KERNELS_ADERDG_OPTIMISED_WEIGHTS_H_
#define _EXAHYPE_KERNELS_ADERDG_OPTIMISED_WEIGHTS_H_

#include <set>

namespace kernels { 
namespace aderdg {
namespace optimised {

void initGaussLegendreNodesAndWeights(const std::set<int>& orders);
void freeGaussLegendreNodesAndWeights(const std::set<int>& orders);

extern double **gaussLegendreNodes;
extern double **gaussLegendreWeights;
extern double *weights1;
extern double *weights2;
extern double *weights3;
}
}
}
#endif /* _EXAHYPE_KERNELS_ADERDG_OPTIMISED_WEIGHTS_H_ */