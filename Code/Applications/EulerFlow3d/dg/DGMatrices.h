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
#ifndef EXAHYPE_DG_DGMATRICES_H
#define EXAHYPE_DG_DGMATRICES_H

namespace exahype {
  namespace dg {
    /** \brief todo */
    extern const double Kxi[4][4];
    /**
     * \brief todo
     **/
    extern const double F0[4];

    /**
     * \brief todo
     **/
    extern const double F0[4];

    /**
     * \brief todo
     **/
    extern const double iK1[4][4];

    /**
     * \brief todo
     **/
    extern const double FLCoeff[4];

    /**
     * \brief todo
     **/
    extern const double FRCoeff[4];
  }
}

#endif /* EXAHYPE_DG_DGMATRICES_H */
