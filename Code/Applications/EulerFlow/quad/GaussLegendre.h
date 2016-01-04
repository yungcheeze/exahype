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
 *  A Gauss-Legendre quadrature with @f$n@f$ nodes integrates a polynomial
 *  of order @f$2\,n-1@f$ exactly.
 *
 *  \note This is one way to manage the quadrature points with the disadvantage that
 *  both large arrays are inserted in each object file including this header file (unless we use the
 *  extern keyword). Other strategies are to use a lookup table or to compute the quadrature points.
 */

///////////////////////////////////////////////////////////////////////////////
#ifndef GAUSSLEGENDRE_H_
#define GAUSSLEGENDRE_H_

namespace exahype {
  namespace quad {
    /** \brief The maximum number of supported nodes.**/
    extern const double gaussLegendreMaxNodes;

    /**
     * \brief The Gauss-Legendre weights w[n][k] in [-1,1] (n=number of nodes, k=node index).
     **/
    extern const double gaussLegendreWeights[10][10];

    /**
     * \brief The Gauss-Legendre nodes x[n][k] in [-1,1] (n=number of nodes, k=node index).
     **/
    extern const double gaussLegendreNodes[10][10];
  }
}

#endif /* GAUSSLEGENDRE_H_ */
