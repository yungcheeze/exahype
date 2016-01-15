#ifndef EXAHYPE_KERNELS_DEFAULT_PDEDESCRIPTION_H_
#define EXAHYPE_KERNELS_DEFAULT_PDEDESCRIPTION_H_

namespace exahype {
  namespace kernels {
    /**
     * @brief Computes the initial value of the kernels at point (x,y)
     *
     * @note: DEC: Needs to be defined by the user.
     * @note: DEC: Needs to provide an enumerator and state information for users for
     *        uncertainty quantification/parameter sweep runs.
     * @note: DEC: Needs to accept user data (void* userData) for
     *        uncertainty quantification/parameter sweep runs.
     *
     * @param[in] x physical x coordinate
     * @param[in] y physical y coordinate
     * @param[in] nvar the number of consered quantities
     * @param[out] values array of size nvar
     */
    void PDEInitialValue2d(const double x,const double y,const int nvar,double * Q);
  }  // namespace kernels
}  // namespace exahype

#endif /* EXAHYPE_KERNELS_DEFAULT_PDEDESCRIPTION_H_ */
