#ifndef EXAHYPE_DG_ADERDG_H_
#define EXAHYPE_DG_ADERDG_H_

namespace exahype {
  namespace dg {
    /**
     * Order depending PNPM factor.
     */
    constexpr double PNPM[10] = {
        1.0, 0.33, 0.17, 0.1, 0.069, 0.045,  0.038, 0.03, 0.02, 0.015
    };

    /**
     * @brief Returns a stable time step size.
     *
     * @param[in] luh       an array of size basisSize**2 containing the solution dof
     * @param[in] dx        an array of size dim containing the extent of the cell
     * @param[in] lambda    a temporary work vector to store the eigenvalues
     * @param[in] nvar      the number of conserved variables
     * @param[in] basisSize the size of the 1D basis
     */
    template <int dim>
    double stableTimeStepSize(
        const double * const luh,
        const double * const dx,
        double * lambda,
        const int nvar,
        const int basisSize
    );

    /**
     * @todo docu
     */
    template <int dim>
    void updateSolution(
        double * luh,
        const double * const lduh,
        const double * const dx,
        const double dt,
        const int nvar,
        const int basisSize
    );

    /**
     * @todo docu
     */
    template <int dim>
    void volumeIntegral(
        double * lduh,
        const double * const lFhi,
        const double * const dx
    );

    /**
     * @todo docu
     * @brief Computes the 2d surface integral contributions
     * to the cell update.
     */
    void surfaceIntegral(
        double * lduh,
        const double * const dx,
        const int nvar,
        const int basisSize,
        const double * const FLeft,
        const double * const FRight,
        const double * const FFront,
        const double * const FBack
    );

    /**
     * @todo docu
     * @brief Computes the 3d surface integral contributions
     * to the cell update.
     */
    void surfaceIntegral(
        double * lduh,
        const double * const dx,
        const int nvar,
        const int basisSize,
        const double * const FLeft,
        const double * const FRight,
        const double * const FFront,
        const double * const FBack,
        const double * const FBottom,
        const double * const FTop
    );

    /**
     * \todo docu
     * Computes the normal fluxes/fluctuations of the two cells adjacent to the current face/
     * Needs the extrapolated predictor values at the interface.
     *
     * QavL,QavR,lambdaL, and lambdaR are work vectors of size nvar.
     */
    template <int dim>
    void solveRiemannProblem(
        double * FL,
        double * FR,
        const double * const QL,
        const double * const QR,
        double * QavL,
        double * QavR,
        double * lambdaL,
        double * lambdaR,
        const double dt,
        const double hFace,
        const double * const n,
        const int nvar,
        const int basisSize
    );

    /**
     * @todo docu
     * Computes the space-time predictor lQi, the space-time volume flux lFi,
     * the predictor lQhi, the volume flux lFhi, the boundary
     * extrapolated predictor lQhbnd and normal flux lFhbnd.
     *
     * rhs0, rhs, and tmp are work vectors.
     * luh will not be modified.
     */
    template <int dim>
    void spaceTimePredictor(
        double * lQi,
        double * lFi,
        const double * const luh, // const
        double * lQhi,
        double * lFhi,
        double * lQhbnd,
        double * lFhbnd,
        double * rhs0,
        double * rhs,
        double * tmp,
        const double * const dx,
        const double dt
    );

    // 2D specialisations
    template <>
    double stableTimeStepSize<2>(
        const double * const luh,
        const double * const dx,
        double * lambda,
        const int nvar,
        const int basisSize
    );

    template <>
    void updateSolution<2>(
        double * luh,
        const double * const lduh,
        const double * const dx,
        const double dt,
        const int nvar,
        const int basisSize
    );

    template <>
    void volumeIntegral<2>(
        double * lduh,
        const double * const lFhi,
        const double * const dx
    );

    template <>
    void solveRiemannProblem<2>(
        double * FL,
        double * FR,
        const double * const QL,
        const double * const QR,
        double * QavL,
        double * QavR,
        double * lambdaL,
        double * lambdaR,
        const double dt,
        const double hFace,
        const double * const n,
        const int nvar,
        const int basisSize
    );

    template <>
    void spaceTimePredictor<2>(
        double * lQi,
        double * lFi,
        const double * const luh,
        double * lQhi,
        double * lFhi,
        double * lQhbnd,
        double * lFhbnd,
        double * rhs0,
        double * rhs,
        double * tmp,
        const double * const dx,
        const double dt
    );

    // 3D specialisations
    template <>
    double stableTimeStepSize<3>(
        const double * const luh,
        const double * const dx,
        double * lambda,
        const int nvar,
        const int basisSize
    );

    template <>
    void updateSolution<3>(
        double * luh,
        const double * const lduh,
        const double * const dx,
        const double dt,
        const int nvar,
        const int basisSize
    );

    template <>
    void volumeIntegral<3>(
        double * lduh,
        const double * const lFhi,
        const double * const dx
    );

    template <>
    void solveRiemannProblem<3>(
        double * FL,
        double * FR,
        const double * const QL,
        const double * const QR,
        double * QavL,
        double * QavR,
        double * lambdaL,
        double * lambdaR,
        const double dt,
        const double hFace,
        const double * const n,
        const int nvar,
        const int basisSize
    );

    template <>
    void spaceTimePredictor<3>(
        double * lQi,
        double * lFi,
        const double * const luh,
        double * lQhi,
        double * lFhi,
        double * lQhbnd,
        double * lFhbnd,
        double * rhs0,
        double * rhs,
        double * tmp,
        const double * const dx,
        const double dt
    );

  }  // namespace dg
}  // namespace exahype


// No includes since there are no partial template function specialisations allowed by the standard
// and full specialisations belong in a source file.

#endif /* EXAHYPE_DG_ADERDG_H_ */
