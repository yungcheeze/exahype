///////////////////////////////////////////////////////////////////////////////
//                    HEADER FILE FOR ExaHyPE CODE                           //
///////////////////////////////////////////////////////////////////////////////
// Authors:   D.E. Charrier                                                  //
// Code:      ExaHyPE                                                        //
// File:      gausslegendre.h                                                //
///////////////////////////////////////////////////////////////////////////////

/** \file gausslegendre.h
 *  \brief Header including Gauss-Legendre quadrature weights and abscissas
 *
 *  A Gauss-Legendre quadrature with \f$n\f$ nodes integrates a polynomial
 *  of order \f$2\,n-1\f$ exactly.
 */

///////////////////////////////////////////////////////////////////////////////
#ifndef GAUSSLEGENDRE_H_
#define GAUSSLEGENDRE_H_

namespace kernels {
/**
 * Initialises the lookup tables \p gaussLegendreWeights
 * and \p gaussLegendreNodes.
 *
 * \see freeGaussLegendreNodesAndWeights
 */
void initGaussLegendreNodesAndWeights(const int maxOrder);

/**
 * Frees the memory that was allocated for the lookup tables \p gaussLegendreWeights
 * and \p gaussLegendreNodes.
 *
 * \see initGaussLegendreNodesAndWeights
 */
void freeGaussLegendreNodesAndWeights(const int maxOrder);

/**
 * The Gauss-Legendre weights mapped onto [0,1]. Array of arrays. The first
 *entry is the order, the second entry the Legendre point.
 **/
extern double** gaussLegendreWeights;

/**
 * The Gauss-Legendre nodes mapped onto [0,1]. Array of arrays. The first entry
 *is the order, the second entry the Legendre point.
 **/
extern double** gaussLegendreNodes;
}

#endif
