///////////////////////////////////////////////////////////////////////////////
//                    HEADER FILE FOR ExaHyPE CODE                           //
///////////////////////////////////////////////////////////////////////////////
// Authors:   D.E. Charrier                                                  //
// Code:      ExaHyPE                                                        //
// File:      DGMatrices.h                                                   //
///////////////////////////////////////////////////////////////////////////////
#ifndef EXAHYPE_ADERDG_DGMATRICES_H
#define EXAHYPE_ADERDG_DGMATRICES_H

#include "exahype/Constants.h"

namespace exahype {
  namespace aderdg {
    /**
      * \brief Element stiffness matrix
      */
    extern const double Kxi[EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];

    /**
     * \brief Time flux matrix (left)
     */
    extern const double F0[EXAHYPE_ORDER+1];

    /**
     * \brief Time flux matrix (right)
     * \note Unused.
     */
    extern const double F1[EXAHYPE_ORDER+1];

    /**
     * \brief Inverse stiffness matrix
     */
    extern const double iK1[EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];

    /**
     * \brief Left extrapolation coefficients
     */
    extern const double FLCoeff[EXAHYPE_ORDER+1];

    /**
     * \brief Right extrapolation coefficients
     */
    extern const double FRCoeff[EXAHYPE_ORDER+1];

    /**
     * \brief Joint extrapolation coefficients
     *
     * FCoeff = [[FLCoeff];[FRCoeff]]
     */
    extern const double FCoeff[2][EXAHYPE_ORDER+1];

    /**
     * \brief derivative operator
     * \note for debugging purposes
     */
    extern const double dudx[EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];

    // is evaluated at compile time
    constexpr int power(int basis, int exp) {
      return (exp==0)? 1 : basis * power(basis, exp-1);
    }
    #define OUTPUTGRID_SIZE power(EXAHYPE_ORDER+1,DIMENSIONS)

    /**
     * \brief Projects nodal values located at Gauss-Legendre nodes to a equidistant grid.
     * \warning Only DIMENSIONS=2 supported!
     */
    extern const double subOutputMatrix[OUTPUTGRID_SIZE][OUTPUTGRID_SIZE];
    #undef OUTPUTGRID_SIZE

  }
}

#endif /* EXAHYPE_ADERDG_DGMATRICES_H */
