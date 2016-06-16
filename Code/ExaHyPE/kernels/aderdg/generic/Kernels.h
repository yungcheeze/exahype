#ifndef _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_
#define _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_

#include "string.h"

#include "peano/utils/Globals.h"
#include "tarch/la/Vector.h"

#include "kernels/GaussLegendreQuadrature.h"

#include "kernels/DGMatrices.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#define MbasisSize 4
#define Mvar 5
#define Mdim 3
#define f2p5(var, dim, i, j, k)                                        \
  (var + Mvar * dim + Mvar * Mdim * i + Mvar * Mdim * MbasisSize * j + \
   Mvar * Mdim * MbasisSize * MbasisSize * k)
#define p2f5(var, dim, i, j, k)                                   \
  (dim * MbasisSize * MbasisSize * MbasisSize * Mvar + Mvar * i + \
   Mvar * MbasisSize * j + Mvar * MbasisSize * MbasisSize * k + var)

#define Mface 6
#define f2p4(var, face, a, b) \
  (var + Mvar * face + Mvar * Mface * a + Mvar * Mface * MbasisSize * b)
#define p2f4(var, face, a, b)                                                 \
  (face * MbasisSize * MbasisSize * Mvar + Mvar * a + Mvar * MbasisSize * b + \
   var)

// todo Dominic Etienne Charrier
// Possibly redundant definition of face indices
// see exahype/solvers/Solver.h
// On the other hand, the kernels should be
// more or less independent of ExaHyPE/exahype.
#define EXAHYPE_FACE_LEFT 0
#define EXAHYPE_FACE_RIGHT 1
#define EXAHYPE_FACE_FRONT 2
#define EXAHYPE_FACE_BACK 3
#define EXAHYPE_FACE_BOTTOM 4
#define EXAHYPE_FACE_TOP 5

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <void PDEFlux(const double* const Q, double** F)>
void spaceTimePredictorNonlinear(
    double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd,
    double* lFhbnd, const double* const luh,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    const double predictorTimeStepSize, const int numberOfVariables,
    const int basisSize);

#if DIMENSIONS == 2

template <void PDEFlux(const double* const Q, double** F)>
void spaceTimePredictorNonlinear(
    double* lQi, double* lFi, const double* const luh,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    const double predictorTimeStepSize, const int numberOfVariables,
    const int basisSize);

#endif  // DIMENSIONS == 2

/**
 * (At the moment, we always evaluate the time averaged space-time
 * predictor unknowns.)
 * todo docu
 */
void predictor(double* lQhi, double* lFhi, const double* const lQi,
               const double* const lFi, const double predictorTimeStepSize,
               const int numberOfVariables, const int basisSize);

/**
 * @todo Dominic Etienne Charrier
 * This is just a "parent" function that
 * invokes the function going by the same
 * name 2*dim times.
 */
void extrapolatedPredictor(double* lQhbnd, double* lFhbnd,
                           const double* const lQhi, const double* const lFhi,
                           const double predictorTimeStepSize,
                           const int numberOfVariables, const int basisSize);

#if DIMENSIONS == 2
/**
 * @todo Dominic Etienne Charrier
 * docu
 * Note that we need to replace lQhi and LFhi by
 * the space-time predictor unknowns if we want to employ
 * local/anarchic time stepping. Since we will have
 * to perform evaluations of the extrapolated boundary fluxes
 * at various appropriate times in this case.
 * The evaluation of the extrapolated predictor requires a time integration
 * of the space-time predictor unknowns. Clearly,
 * \p evaluationTimeStepSize must be smaller than or equal to
 * \p predictorTimeStepSize.
 * At the moment, we always evaluate the time averaged space-time
 * predictor unknowns. Thus it is not necessary to pass these values.
 */
void extrapolatedPredictorXDirection(
    double* lQhbnd, double* lFhbnd, const double* const lQhi,
    const double* const lFhi,
    const int facePosition,  // 0 for "left" face, 1 far "right" face
    const double evaluationTimeStepSize, const double predictorTimeStepSize,
    const int numberOfVariables, const int basisSize);

void extrapolatedPredictorYDirection(
    double* lQhbnd, double* lFhbnd, const double* const lQhi,
    const double* const lFhi,
    const int facePosition,  // 0 for "left" face, 1 far "right" face
    const double evaluationTimeStepSize, const double predictorTimeStepSize,
    const int numberOfVariables, const int basisSize);
#endif

// todo Dominic Etienne Charrier:
// The DIMENSIONS depending mesh size vector enables overloading at the moment.
// If we replace it by scalar mesh size, we have to add a template argument "int
// dim".

void solutionUpdate(double* luh, const double* const lduh, const double dt,
                    const int numberOfVariables, const int basisSize);

void volumeIntegralLinear(double* lduh, const double* const lFhi,
                          const tarch::la::Vector<DIMENSIONS, double>& dx,
                          const int numberOfVariables, const int basisSize);

void volumeIntegralNonlinear(double* lduh, const double* const lFhi,
                             const tarch::la::Vector<DIMENSIONS, double>& dx,
                             const int numberOfVariables, const int basisSize);

// todo 10/02/16: Dominic
// Keep only one surfaceIntegral.
void surfaceIntegralNonlinear(double* lduh, const double* const lFbnd,
                              const tarch::la::Vector<DIMENSIONS, double>& dx,
                              const int numberOfVariables, const int basisSize);

void surfaceIntegralLinear(double* lduh, const double* const lFbnd,
                           const tarch::la::Vector<DIMENSIONS, double>& dx,
                           const int numberOfVariables, const int basisSize);

/*void surfaceIntegral2(
    double* lduh,
    const double* const lFhbnd,
    const tarch::la::Vector<DIMENSIONS,double>&  dx,
    const int numberOfVariables,
    const int basisSize
);*/

#if DIMENSIONS == 2
void surfaceIntegralXDirection(
    double* lduh, const double* const lFhbnd, const double dx,
    const int facePosition,   // 0 for "left" face, 1 for "right" face.
    const double updateSign,  // -1 for "left" face, 1 for "right" face.
    const int numberOfVariables, const int basisSize);

void surfaceIntegralYDirection(
    double* lduh, const double* const lFhbnd, const double dy,
    const int facePosition,   // 0 for "left" face, 1 for "right" face.
    const double updateSign,  // -1 for "left" face, 1 for "right" face.
    const int numberOfVariables, const int basisSize);

#endif

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <void PDESolutionAdjustment(const double* const x, const double J_w,
                                     const double t, const double dt,
                                     double* Q)>
void solutionAdjustment(double* luh,
                        const tarch::la::Vector<DIMENSIONS, double>& center,
                        const tarch::la::Vector<DIMENSIONS, double>& dx,
                        const double t, const double dt,
                        const int numberOfVariables, const int basisSize);

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments
// template argument functions and non-template argument function.
template <void PDEEigenvalues(const double* const Q, const int normalNonZero,
                              double* lambda)>
void riemannSolverNonlinear(double* FL, double* FR, const double* const QL,
                            const double* const QR, const double dt,
                            const int normalNonZero,
                            const int numberOfVariables, const int basisSize);

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <void PDEEigenvalues(const double* const Q, const int normalNonZero,
                              double* lambda)>
double stableTimeStepSize(const double* const luh,
                          const tarch::la::Vector<DIMENSIONS, double>& dx,
                          const int numberOfVariables, const int basisSize);

void faceUnknownsProlongation(
    double* lQhbndFine, double* lFhbndFine, const double* lQhbndCoarse,
    const double* lFhbndCoarse, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize);

void faceUnknownsRestriction(
    double* lQhbndCoarse, double* lFhbndCoarse, const double* lQhbndFine,
    const double* lFhbndFine, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize);

void volumeUnknownsProlongation(
    double* luhFine, const double* luhCoarse, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize);

void volumeUnknownsRestriction(
    double* luhCoarse, const double* luhFine, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize);
}  // namespace c
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels

#if DIMENSIONS == 2
#include "kernels/aderdg/generic/c/2d/riemannSolverLinear.cpph"
#include "kernels/aderdg/generic/c/2d/riemannSolverNonlinear.cpph"
#include "kernels/aderdg/generic/c/2d/solutionAdjustment.cpph"
#include "kernels/aderdg/generic/c/2d/spaceTimePredictorLinear.cpph"
#include "kernels/aderdg/generic/c/2d/spaceTimePredictorNonlinear.cpph"
#include "kernels/aderdg/generic/c/2d/stableTimeStepSize.cpph"
#elif DIMENSIONS == 3
#include "kernels/aderdg/generic/c/3d/riemannSolverLinear.cpph"
#include "kernels/aderdg/generic/c/3d/riemannSolverNonlinear.cpph"
#include "kernels/aderdg/generic/c/3d/solutionAdjustment.cpph"
#include "kernels/aderdg/generic/c/3d/spaceTimePredictorLinear.cpph"
#include "kernels/aderdg/generic/c/3d/spaceTimePredictorNonlinear.cpph"
#include "kernels/aderdg/generic/c/3d/stableTimeStepSize.cpph"
#endif

namespace kernels {
namespace aderdg {
namespace generic {
namespace fortran {

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <void PDEFlux(const double* const Q, double** F)>
void spaceTimePredictorNonlinear(
    double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd,
    double* lFhbnd, const double* const luh,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    const double predictorTimeStepSize, const int numberOfVariables,
    const int basisSize);

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <void PDENCP(const double* const Q, const double* const gradQ,
                      double* BgradQ)>
void spaceTimePredictorLinear(double* lQi, double* lFi, double* lQhi,
                              double* lFhi, double* lQhbnd, double* lFhbnd,
                              const double* const luh,
                              const tarch::la::Vector<DIMENSIONS, double>& dx,
                              const double predictorTimeStepSize,
                              const int numberOfVariables, const int basisSize);

/**
 * (At the moment, we always evaluate the time averaged space-time
 * predictor unknowns.)
 * todo docu
 */
void predictor(double* lQhi, double* lFhi, const double* const lQi,
               const double* const lFi, const double predictorTimeStepSize,
               const int numberOfVariables, const int basisSize);

/**
 * @todo Dominic Etienne Charrier
 * This is just a "parent" function that
 * invokes the function going by the same
 * name 2*dim times.
 */
void extrapolatedPredictor(double* lQhbnd, double* lFhbnd,
                           const double* const lQhi, const double* const lFhi,
                           const double predictorTimeStepSize,
                           const int numberOfVariables, const int basisSize);

// todo Dominic Etienne Charrier:
// The DIMENSIONS depending mesh size vector enables overloading at the moment.
// If we replace it by scalar mesh size, we have to add a template argument "int
// dim".

void solutionUpdate(double* luh, const double* const lduh, const double dt,
                    const int numberOfVariables, const int basisSize);

void volumeIntegralNonlinear(double* lduh, const double* const lFhi,
                             const tarch::la::Vector<DIMENSIONS, double>& dx,
                             const int numberOfVariables, const int basisSize);

void volumeIntegralLinear(double* lduh, const double* const lFhi,
                          const tarch::la::Vector<DIMENSIONS, double>& dx,
                          const int numberOfVariables, const int basisSize);

// todo 10/02/16: Dominic
// Keep only one surfaceIntegral.
void surfaceIntegralNonlinear(double* lduh, const double* const lFbnd,
                              const tarch::la::Vector<DIMENSIONS, double>& dx,
                              const int numberOfVariables, const int basisSize);

void surfaceIntegralLinear(double* lduh, const double* const lFbnd,
                           const tarch::la::Vector<DIMENSIONS, double>& dx,
                           const int numberOfVariables, const int basisSize);

/*void surfaceIntegral2(
    double* lduh,
    const double* const lFhbnd,
    const tarch::la::Vector<DIMENSIONS,double>&  dx,
    const int numberOfVariables,
    const int basisSize
);*/

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <void PDESolutionAdjustment(const double* const x, const double w,
                                     const double t, const double dt,
                                     double* Q)>
void solutionAdjustment(double* luh,
                        const tarch::la::Vector<DIMENSIONS, double>& center,
                        const tarch::la::Vector<DIMENSIONS, double>& dx,
                        const double t, const double dt,
                        const int numberOfVariables, const int basisSize);

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments
// template argument functions and non-template argument function.
template <void PDEEigenvalues(const double* const Q, const int normalNonZero,
                              double* lambda)>
void riemannSolverNonlinear(double* FL, double* FR, const double* const QL,
                            const double* const QR, const double dt,
                            const int normalNonZero,
                            const int numberOfVariables, const int basisSize);

template <void PDEEigenvalues(const double* const Q, const int normalNonZero,
                              double* lambda),
          void PDEMatrixB(const double* const Q, const int normalNonZero,
                          double* Bn)>
void riemannSolverLinear(double* FL, double* FR, const double* const QL,
                         const double* const QR, const double dt,
                         const int normalNonZero, const int numberOfVariables,
                         const int basisSize);

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <void PDEEigenvalues(const double* const Q, const int normalNonZero,
                              double* lambda)>
double stableTimeStepSize(const double* const luh,
                          const tarch::la::Vector<DIMENSIONS, double>& dx,
                          const int numberOfVariables, const int basisSize);

void faceUnknownsProlongation(
    double* lQhbndFine, double* lFhbndFine, const double* lQhbndCoarse,
    const double* lFhbndCoarse, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize);

void faceUnknownsRestriction(
    double* lQhbndCoarse, double* lFhbndCoarse, const double* lQhbndFine,
    const double* lFhbndFine, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize);

void volumeUnknownsProlongation(
    double* luhFine, const double* luhCoarse, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize);

void volumeUnknownsRestriction(
    double* luhCoarse, const double* luhFine, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize);

}  // namespace fortran
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels

#if DIMENSIONS == 3
#include "kernels/aderdg/generic/fortran/3d/riemannSolverLinear.cpph"
#include "kernels/aderdg/generic/fortran/3d/riemannSolverNonlinear.cpph"
#include "kernels/aderdg/generic/fortran/3d/solutionAdjustment.cpph"
#include "kernels/aderdg/generic/fortran/3d/spaceTimePredictorLinear.cpph"
#include "kernels/aderdg/generic/fortran/3d/spaceTimePredictorNonlinear.cpph"
#include "kernels/aderdg/generic/fortran/3d/stableTimeStepSize.cpph"
// #elif DIMENSIONS == 2
// //@todo
// #include "kernels/aderdg/generic/fortran/2d/solutionAdjustment.cpph"
// #include "kernels/aderdg/generic/fortran/2d/stableTimeStepSize.cpph"
// #include "kernels/aderdg/generic/fortran/2d/spaceTimePredictorNonlinear.cpph"
// #include "kernels/aderdg/generic/fortran/2d/spaceTimePredictorLinear.cpph"
// #include "kernels/aderdg/generic/fortran/2d/riemannSolverNonlinear.cpph"
// #include "kernels/aderdg/generic/fortran/2d/riemannSolverLinear.cpph"
#endif

#endif
