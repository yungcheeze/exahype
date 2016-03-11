#ifndef EXAHYPE_DG_ADERDG_H_
#define EXAHYPE_DG_ADERDG_H_

namespace exahype {
namespace dg {
/**
 * Order depending PNPM factor.
 */
constexpr double PNPM[10] = {1.0,   0.33,  0.17, 0.1,  0.069,
                             0.045, 0.038, 0.03, 0.02, 0.015};

/**
 * @brief Returns a stable time step size.
 *
 * @param[in] luh       an array of size basisSize**2 containing the solution
 *dof
 * @param[in] dx        an array of size dim containing the extent of the cell
 * @param[in] lambda    a temporary work vector to store the eigenvalues
 * @param[in] nvar      the number of conserved variables
 * @param[in] basisSize the size of the 1D basis
 */
template <int dim>
double stableTimeStepSize(const double* restrict const luh,
                          const double* restrict const dx,
                          double* restrict lambda);

/**
 * @brief Adds the delta to the solution and thereby finalises the updates of
 * the cell
 * @param[in,out] luh    the current solution that is to be updated
 * @param[in]     lduh   the update quantities that are incorporated into \p luh
 * @param[in]     dt     the time step
 */
template <int dim>
void updateSolution(double* restrict luh, const double* restrict const lduh,
                    const double dt);

/**
 * @brief Computes the volume integral
 *
 * @param[out] lduh   spatial degrees of freedom, array of size
 *nVar*basisSize**dim
 * @param[in]  lFhi   flux tensor, array of size nVar*dim*basisSize**dim
 * @param[in]  dx     array of size dim containing the spacing of the cell
 */
template <int dim>
void volumeIntegral(double* restrict lduh, const double* restrict const lFhi,
                    const double* restrict const dx);

/**
 * @todo docu
 * @brief Computes the 2d surface integral contributions
 * to the cell update.
 */
void surfaceIntegral(double* restrict lduh, const double* restrict const dx,
                     const double* restrict const FLeft,
                     const double* restrict const FRight,
                     const double* restrict const FFront,
                     const double* restrict const FBack);

/**
 * @todo docu
 * @brief Computes the 3d surface integral contributions
 * to the cell update.
 */
void surfaceIntegral(double* restrict lduh, const double* restrict const dx,
                     const double* restrict const FLeft,
                     const double* restrict const FRight,
                     const double* restrict const FFront,
                     const double* restrict const FBack,
                     const double* restrict const FBottom,
                     const double* restrict const FTop);

/**
 * \todo docu
 * Computes the normal fluxes/fluctuations of the two cells adjacent to the
 *current face/
 * Needs the extrapolated predictor values at the interface.
 *
 * QavL,QavR,lambdaL, and lambdaR are work vectors of size nvar.
 */
template <int dim>
void solveRiemannProblem(double* restrict FL, double* restrict FR,
                         const double* restrict const QL,
                         const double* restrict const QR, double* restrict QavL,
                         double* restrict QavR, double* restrict lambdaL,
                         double* restrict lambdaR, const double dt,
                         const double hFace, const double* restrict const n);

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
void spaceTimePredictor(double* restrict lQi, double* restrict lFi,
                        const double* restrict const luh,  // const
                        double* restrict lQhi, double* restrict lFhi,
                        double* restrict lQhbnd, double* restrict lFhbnd,
                        double* restrict rhs0, double* restrict rhs,
                        double* restrict tmp, const double* restrict const dx,
                        const double dt);

// 2D specialisations
template <>
double stableTimeStepSize<2>(const double* restrict const luh,
                             const double* restrict const dx,
                             double* restrict lambda);

template <>
void updateSolution<2>(double* restrict luh, const double* restrict const lduh,
                       const double dt);

template <>
void volumeIntegral<2>(double* restrict lduh, const double* restrict const lFhi,
                       const double* restrict const dx);

template <>
void solveRiemannProblem<2>(double* restrict FL, double* restrict FR,
                            const double* restrict const QL,
                            const double* restrict const QR,
                            double* restrict QavL, double* restrict QavR,
                            double* restrict lambdaL, double* restrict lambdaR,
                            const double dt, const double hFace,
                            const double* restrict const n);

template <>
void spaceTimePredictor<2>(double* restrict lQi, double* restrict lFi,
                           const double* restrict const luh,
                           double* restrict lQhi, double* restrict lFhi,
                           double* restrict lQhbnd, double* restrict lFhbnd,
                           double* restrict rhs0, double* restrict rhs,
                           double* restrict tmp,
                           const double* restrict const dx, const double dt);

// 3D specialisations
template <>
double stableTimeStepSize<3>(const double* restrict const luh,
                             const double* restrict const dx,
                             double* restrict lambda);

template <>
void updateSolution<3>(double* restrict luh, const double* restrict const lduh,
                       const double dt);

template <>
void volumeIntegral<3>(double* restrict lduh, const double* restrict const lFhi,
                       const double* restrict const dx);

template <>
void solveRiemannProblem<3>(double* restrict FL, double* restrict FR,
                            const double* restrict const QL,
                            const double* restrict const QR,
                            double* restrict QavL, double* restrict QavR,
                            double* restrict lambdaL, double* restrict lambdaR,
                            const double dt, const double hFace,
                            const double* restrict const n);

template <>
void spaceTimePredictor<3>(double* restrict lQi, double* restrict lFi,
                           const double* restrict const luh,
                           double* restrict lQhi, double* restrict lFhi,
                           double* restrict lQhbnd, double* restrict lFhbnd,
                           double* restrict rhs0, double* restrict rhs,
                           double* restrict tmp,
                           const double* restrict const dx, const double dt);

}  // namespace dg
}  // namespace exahype

// No includes since there are no partial template function specialisations
// allowed by the standard
// and full specialisations belong in a source file.

#endif /* EXAHYPE_DG_ADERDG_H_ */
