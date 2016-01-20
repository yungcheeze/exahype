/**
 * Declares flux functions for a elementwise constant DG method.
 */

#ifndef EXAHYPE_EULERFLOW3D_PROBLEM_H_
#define EXAHYPE_EULERFLOW3D_PROBLEM_H_

#define GAMMA 1.4

namespace exahype {
  namespace problem {
    /**
     * @brief Computes the initial value of the problem at point (x,y)
     *
     * @param[in] x physical x coordinate
     * @param[in] y physical y coordinate
     * @param[out] values array of size nvar
     */
    void PDEInitialValue2d(const double x,
                           const double y,
                           double * restrict Q);

    /**
     * @brief Computes the physical (flux) tensor values for the 2D Euler equations.
     *
     * @param[in]  Q the conserved variables
     * @param[out] f x component of the physical flux
     * @param[out] g x component of the physical flux
     */
    void PDEFlux(const double * restrict const Q,
                 double * restrict f,
                 double * restrict g);

    /**
     * @brief Computes the physical (flux) tensor values for the 3D Euler equations.
     *
     * @param[in]  Q the conserved variables
     * @param[out] f x component of the physical flux
     * @param[out] g y component of the physical flux
     * @param[out] h z component of the physical flux
     */
    void PDEFlux(const double * restrict const Q,
                 double * restrict f,
                 double * restrict g,
                 double * restrict h);

    /**
     * @brief Computes the physical (flux) tensor values.
     *
     * @param[in]  Q      the conserved variables
     * @param[in]  n      normal vector
     * @param[out] lambda vector of eigenvalues
     */
    void PDEEigenvalues(const double * restrict const Q,
                        const double * restrict const n,
                        double * restrict lambda);

    /**
     * @note Unused. Leftover from debugging.
     */
    double PDENormalFlux(const double x,const double y,const double nx,const double ny);

    /**
     * @note Unused. Leftover from debugging.
     */
    double DGNormalFlux(const double x,const double y,const double nx,const double ny);

    /**
     * @note Unused. Leftover from debugging.
     */
    void   DGRiemannSolver(const double x,const double y,const double nx,const double ny,double *selfFlux,double *neighbourFlux);
  }  // namespace problem
}  // namespace exahype

#endif /* EXAHYPE_EULERFLOW3D_PROBLEM_H_ */
