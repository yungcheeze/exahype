#ifndef _EXAHYPE_KERNELS_ADERDG_OPTIMISED_DGMATRICES_H_
#define _EXAHYPE_KERNELS_ADERDG_OPTIMISED_DGMATRICES_H_

#include <set>

namespace kernels { 
namespace aderdg {
namespace optimised {

void initDGMatrices(const std::set<int>& orders);
void freeDGMatrices(const std::set<int>& orders);

extern double *Kxi;
extern double *Kxi_T;
extern double *iK1;
extern double *dudx;
extern double *s_m;
extern double *s_v;
extern double *tmp_bnd;
extern double *F0; 
extern double *FLCoeff;
extern double *FRCoeff;
extern double ***equidistantGridProjector1d;
extern double **** fineGridProjector1d;

}
}
}
#endif /* _EXAHYPE_KERNELS_ADERDG_OPTIMISED_DGMATRICES_H_ */