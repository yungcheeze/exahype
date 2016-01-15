#ifndef EXAHYPE_KERNELS_DEFAULT_PDEDESCRIPTION_H_
#define EXAHYPE_KERNELS_DEFAULT_PDEDESCRIPTION_H_

#define GAMMA 1.4

namespace exahype {
  namespace pde {
    /**
     * @brief Computes the physical (flux) tensor values.
     *
     * @note: DEC: Function does probably not exist if optimised kernel for single PDE runs is used.
     *             Function must be provided by user for  uncertainty quantification/parameter sweep runs.
     * @note: DEC: Needs to be dynamically modified in the uncertainty quantification/parameter sweep case.
     * @note: DEC: Needs to accept user data (void* userData) for
     *             uncertainty quantification/parameter sweep runs.
     *
     * @param[in] Q the conserved variables
     * @param[in] nvar the number of consered quantities
     * @param[out] f x component of the physical flux
     * @param[out] g x component of the physical flux
     */
    void PDEFlux(const double * const Q,const int nvar,double * f,double * g);


    /**
     * @brief Computes the physical (flux) tensor values.
     *
     * @note: DEC: Function does possibly not exist if optimised kernel is used.
     *
     * @param[in] Q the conserved variables
     * @param[in] nvar the number of consered quantities
     * @param[n]  normal vector
     * @param[n]  d the length of the normal vector, i.e., the space dimension
     * @param[in] nvar the number of consered quantities
     * @param[out] f lambda
     */
    void PDEEigenvalues(const double * const Q,const int nvar,const double * const n,const int d,double * lambda);
  }  // namespace kernels
}  // namespace exahype

#endif /* EXAHYPE_KERNELS_DEFAULT_PDEDESCRIPTION_H_ */
