#ifndef _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_
#define _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_

//#define GAMMA 1.4


namespace kernels {
  namespace aderdg {
    namespace generic {
      template <void (*PDEFlux2d)(const double * const Q,const int nvar,double * f,double * g)>
      void spaceTimePredictor(
        double * lQi,
        double * lFi,
        const double * const luh,
        double * lQhi,
        double * lFhi,
        double * lQhbnd,
        double * lFhbnd,
        const double * const dx,
        const double dt,
        int          order,
        int          numberOfVariables
      );
    }
  }
}


#include "kernels/aderdg/generic/Kernels.cpph"

#endif
