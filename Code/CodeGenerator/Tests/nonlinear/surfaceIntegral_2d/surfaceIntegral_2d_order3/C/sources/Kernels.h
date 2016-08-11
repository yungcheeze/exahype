#ifndef EXAHYPE_KERNELS_OPTIMISED_KERNELS_H_
#define EXAHYPE_KERNELS_OPTIMISED_KERNELS_H_


#define NDEBUG
namespace kernels {
  namespace aderdg {
    namespace optimised {

     
      void surfaceIntegral( 
        double* restrict lduh, 
        const double* restrict const lFbnd, 
        const double* dx
      );

    }
  }
}


#endif // EXAHYPE_KERNELS_OPTIMISED_KERNELS_H_