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
 */

///////////////////////////////////////////////////////////////////////////////
#ifndef GAUSSLEGENDRE_H_
#define GAUSSLEGENDRE_H_

#include "EulerFlow/Constants.h"

namespace exahype {
  namespace quad {
    /** \brief The maximum number of supported nodes.**/
    extern const double gaussLegendreMaxNodes;

    /**
     * \brief The Gauss-Legendre weights mapped onto [0,1]
     **/
    extern const double gaussLegendreWeights[EXAHYPE_ORDER+1];

    /**
     * \brief The Gauss-Legendre nodes mapped onto [0,1]
     **/
    extern const double gaussLegendreNodes[EXAHYPE_ORDER+1];
  }
}

#endif /* GAUSSLEGENDRE_H_ */
