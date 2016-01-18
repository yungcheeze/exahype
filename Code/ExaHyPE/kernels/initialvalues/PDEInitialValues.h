#ifndef EXAHYPE_PDE_PDEINITIALVALUES_H_
#define EXAHYPE_PDE_PDEINITIALVALUES_H_

namespace exahype {
  namespace pde {
    /**
     * @brief Provides the initial value of the kernels at point (x,y)
     *
     * @todo: DEC: Needs to be defined by the user.
     * @todo: DEC: Needs to provide an enumerator and state information for
     *        uncertainty quantification/parameter sweep runs.
     * @todo: DEC: Needs to accept user data (void* userData) for
     *        uncertainty quantification/parameter sweep runs.
     *
     * @param[in]  x    physical x coordinate.
     * @param[in]  y    physical y coordinate.
     * @param[in]  nvar the number of conserved quantities.
     * @param[out] Q    values array of size nvar.
     */
    void PDEInitialValue2d(const double x,const double y,const int nvar,double * Q);

    /**
     * @brief Provides the initial value of the kernels at point (x,y,z).
     *
     * @todo: DEC: Needs to be defined by the user.
     * @todo: DEC: Needs to provide an enumerator and state information for
     *        uncertainty quantification/parameter sweep runs.
     * @todo: DEC: Needs to accept user data (void* userData) for
     *        uncertainty quantification/parameter sweep runs.
     *
     * @param[in]  x    physical x coordinate.
     * @param[in]  y    physical y coordinate.
     * @param[in]  nvar the number of conserved quantities.
     * @param[out] Q    values array of size nvar.
     */
    void PDEInitialValue3d(const double x,const double y,const double z,const int nvar,double * Q);
  }  // namespace pde
}  // namespace exahype

#endif /* EXAHYPE_PDE_PDEINITIALVALUES_H_ */
