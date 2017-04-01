/*
 * rusanov.h
 *
 *  Created on: Apr 1, 2017
 *      Author: dominic
 */

#ifndef KERNELS_FINITEVOLUMES_RIEMANNSOLVERS_C_RIEMANNSOLVERS_H_
#define KERNELS_FINITEVOLUMES_RIEMANNSOLVERS_C_RIEMANNSOLVERS_H_

namespace kernels {
namespace finitevolumes {
namespace riemannsolvers {
namespace c {

/**
 * A simple Rusanov flux considering pointwise
 * left and right values.
 *
 * \note This does not result in a well-balanced scheme.
 * It is not a good Riemann solver for the Shallow Water Equations (SWE) e.g.
 */
template <typename SolverType>
double rusanov(SolverType& solver, double* fL, double *fR, const double* qL, const double* qR,
               int normalNonZero,double** tempStateSizedVectors);
} // namespace c
} // namespace riemansolvers
} // namespace finitevolumes
} // namespace kernels

#include "rusanov.cpph"

#endif /* KERNELS_FINITEVOLUMES_RIEMANNSOLVERS_C_RIEMANNSOLVERS_H_ */
