#ifndef EXAHYPE_KERNELS_OPTIMISED_KERNELS_H_
#define EXAHYPE_KERNELS_OPTIMISED_KERNELS_H_


#define NDEBUG
namespace kernels {
  namespace aderdg {
    namespace optimised {


      void volumeIntegral( 
        double* restrict lduh, 
        const double* restrict const lFhi, 
        const double* dx
      );

    

    }
  }
}


#endif // EXAHYPE_KERNELS_OPTIMISED_KERNELS_H_
