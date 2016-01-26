///////////////////////////////////////////////////////////////////////////////
//                    HEADER FILE FOR ExaHyPE CODE                           //
///////////////////////////////////////////////////////////////////////////////
// Authors:   D.E. Charrier                                                  //
// Code:      ExaHyPE                                                        //
// File:      DGMatrices.h                                                   //
///////////////////////////////////////////////////////////////////////////////
#ifndef EXAHYPE_KERNELS_DGMATRICES_H
#define EXAHYPE_KERNELS_DGMATRICES_H



#define EXAHYPE_ORDER 9

namespace kernels {
  void initDGMatrices();

  /**
      * \brief Element stiffness matrix
      */
    extern const double Kxi[EXAHYPE_ORDER+1][EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];

    /**
     * \brief Time flux matrix (left)
     */
    extern const double F0[EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];

    /**
     * \brief Time flux matrix (right)
     * \note Unused.
     */
    extern const double F1[EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];

    /**
     * \brief Inverse stiffness matrix
     */
    extern const double iK1[EXAHYPE_ORDER+1][EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];

    /**
     * \brief Left extrapolation coefficients
     */
    extern const double FLCoeff[EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];

    /**
     * \brief Right extrapolation coefficients
     */
    extern const double FRCoeff[EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];

    /**
     * \brief Joint extrapolation coefficients
     *
     * FCoeff = [[FLCoeff];[FRCoeff]]
     */
    extern const double FCoeff[EXAHYPE_ORDER+1][2][EXAHYPE_ORDER+1];

    /**
     * \brief derivative operator
     * \note for debugging purposes
     */
    extern const double dudx[EXAHYPE_ORDER+1][EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];
}

#endif
