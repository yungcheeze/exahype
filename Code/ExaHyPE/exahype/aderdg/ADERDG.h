#ifndef EXAHYPE_DG_ADERDG_H_
#define EXAHYPE_DG_ADERDG_H_

/**
 * Indexing of the faces of a cell.
 */
///@{
#define EXAHYPE_FACE_LEFT   0
#define EXAHYPE_FACE_RIGHT  1
#define EXAHYPE_FACE_FRONT  2
#define EXAHYPE_FACE_BACK   3
#define EXAHYPE_FACE_BOTTOM 4
#define EXAHYPE_FACE_TOP    5
///@}

namespace exahype {
  namespace aderdg {
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
    void exahype::aderdg::solutionUpdate(
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
    void exahype::aderdg::volumeIntegral(
        double * lduh,
        const double * const lFhi,
        const double * const dx
    );

    /**
     * @todo docu
     * @brief Computes the 2d surface integral contributions
     * to the cell update.
     */
    void exahype::aderdg::surfaceIntegral(
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
    void exahype::aderdg::surfaceIntegral(
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
    void exahype::aderdg::riemannSolver(
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
        const double * const n
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
    void exahype::aderdg::spaceTimePredictor(
        double * lQi,
        double * lFi,
        const double * const luh, // const
        double * lQhi,
        double * lFhi,
        double * lQhbnd[],
        double * lFhbnd[],
        double * rhs0,
        double * rhs,
        double * tmp,
        const double * const dx,
        const double dt
    );

    /**
     * Sets the initial values for an element.
     *
     * @param[in/out] luh       The local solution DoF.
     * @param[in]     center    The element center.
     * @param[in]     nvar      Number of physical coordinates.
     * @param[in]     basisSize number of DoF in one coordinate direction.
     *
     * @todo: DEC: Will need more input parameters for parameter sweep runs (void* userData).
     */
    template <int dim>
    void exahype::aderdg::initialValues(
        double * luh,
        const double * const center,
        const double * const dx,
        const int nvar,
        const int basisSize
    );
  }  // namespace aderdg
}  // namespace exahype

// 2D specialisations
template <>
double exahype::aderdg::stableTimeStepSize<2>(
    const double * const luh,
    const double * const dx,
    double * lambda,
    const int nvar,
    const int basisSize
);

template <>
void exahype::aderdg::solutionUpdate<2>(
    double * luh,
    const double * const lduh,
    const double * const dx,
    const double dt,
    const int nvar,
    const int basisSize
);

template <>
void exahype::aderdg::volumeIntegral<2>(
    double * lduh,
    const double * const lFhi,
    const double * const dx
);

template <>
void exahype::aderdg::riemannSolver<2>(
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
    const double * const n
);

template <>
void exahype::aderdg::spaceTimePredictor<2>(
    double * lQi,
    double * lFi,
    const double * const luh,
    double * lQhi,
    double * lFhi,
    double * lQhbnd[],
    double * lFhbnd[],
    double * rhs0,
    double * rhs,
    double * tmp,
    const double * const dx,
    const double dt
);

template <>
void exahype::aderdg::initialValues<2>(
    double * luh,
    const double * const center,
    const double * const dx,
    const int nvar,
    const int basisSize
);

// 3D specialisations
template <>
double exahype::aderdg::stableTimeStepSize<3>(
    const double * const luh,
    const double * const dx,
    double * lambda,
    const int nvar,
    const int basisSize
);

template <>
void exahype::aderdg::solutionUpdate<3>(
    double * luh,
    const double * const lduh,
    const double * const dx,
    const double dt,
    const int nvar,
    const int basisSize
);

template <>
void exahype::aderdg::volumeIntegral<3>(
    double * lduh,
    const double * const lFhi,
    const double * const dx
);

template <>
void exahype::aderdg::riemannSolver<3>(
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
    const double * const n
);

template <>
void exahype::aderdg::spaceTimePredictor<3>(
    double * lQi,
    double * lFi,
    const double * const luh,
    double * lQhi,
    double * lFhi,
    double * lQhbnd[],
    double * lFhbnd[],
    double * rhs0,
    double * rhs,
    double * tmp,
    const double * const dx,
    const double dt
);

template <>
void exahype::aderdg::initialValues<3>(
    double * luh,
    const double * const center,
    const double * const dx,
    const int nvar,
    const int basisSize
);

// No includes since there are no partial template function specialisations allowed by the standard
// and full specialisations belong in a source file.

#endif /* EXAHYPE_DG_ADERDG_H_ */
