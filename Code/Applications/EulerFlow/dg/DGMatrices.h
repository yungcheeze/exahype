///////////////////////////////////////////////////////////////////////////////
//                    HEADER FILE FOR ExaHyPE CODE                           //
///////////////////////////////////////////////////////////////////////////////
// Authors:   D.E. Charrier                                                  //
// Code:      ExaHyPE                                                        //
// File:      DGMatrices.h                                                   //
///////////////////////////////////////////////////////////////////////////////
#ifndef EXAHYPE_DG_DGMATRICES_H
#define EXAHYPE_DG_DGMATRICES_H

namespace exahype {
  namespace dg {
    /** @brief todo */
    extern const double Kxi[4][4];
    /**
     * @brief todo
     **/
    extern const double F0[4];

    /**
     * @brief todo
     **/
    extern const double F0[4];

    /**
     * @brief todo
     **/
    extern const double iK1[4][4];

    /**
     * @brief todo
     **/
    extern const double FLCoeff[4];

    /**
     * @brief todo
     **/
    extern const double FRCoeff[4];


    /*
     * @brief todo
     */
    extern const double dudx[4][4];


    /*
     * Projects nodal values located at Gauss-Legendre nodes to a equidistant grid.
     */
    extern const double subOutputMatrix[16][16];
  }
}

#endif /* EXAHYPE_DG_DGMATRICES_H */
