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
    // order,row,column
    // [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1][EXAHYPE_ORDER+1]
    extern const double*** Kxi;

    /**
     * \brief Time flux matrix (left)
     */
    // order, row
    // [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1]
    extern const double** F0;

    /**
     * \brief Time flux matrix (right)
     * \note Unused.
     */
    // order, row
    // [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1]
    extern const double** F1;

    /**
     * \brief Inverse stiffness matrix
     */
    // order, row, column
    // [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1][EXAHYPE_ORDER+1]
    extern const double*** iK1;

    /**
     * \brief Left extrapolation coefficients
     */
    // order, row,
    // [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1]
    extern const double** FLCoeff;

    /**
     * \brief Right extrapolation coefficients
     */
    // order, row
    // [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];
    extern const double** FRCoeff;

    /**
     * \brief Joint extrapolation coefficients
     *
     * FCoeff = [[FLCoeff];[FRCoeff]]
     */
    // order, left/right, row
    // [EXAHYPE_ORDER+1][2][EXAHYPE_ORDER+1];
    extern const double*** FCoeff;

    /**
     * \brief derivative operator
     * \note for debugging purposes
     */
    // order, row, column
    // [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];
    extern const double*** dudx;
}

#endif
