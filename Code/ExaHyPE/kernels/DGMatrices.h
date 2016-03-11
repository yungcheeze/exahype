///////////////////////////////////////////////////////////////////////////////
//                    HEADER FILE FOR ExaHyPE CODE                           //
///////////////////////////////////////////////////////////////////////////////
// Authors:   D.E. Charrier                                                  //
// Code:      ExaHyPE                                                        //
// File:      DGMatrices.h                                                   //
///////////////////////////////////////////////////////////////////////////////
#ifndef EXAHYPE_KERNELS_DGMATRICES_H
#define EXAHYPE_KERNELS_DGMATRICES_H

namespace kernels {
void initDGMatrices();

/**
 * \brief Element stiffness matrix
 */
// todo Dominic Etienne Charrier
// order,row,column
// [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1][EXAHYPE_ORDER+1]
extern double*** Kxi;

/**
 * \brief Time flux matrix (left)
 */
// todo Dominic Etienne Charrier
// order, row
// [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1]
extern double** F0;

/**
 * \brief Time flux matrix (right)
 * \note Unused.
 */
// todo Dominic Etienne Charrier
// order, row
// [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1]
extern double** F1;

/**
 * \brief Inverse stiffness matrix
 */
// todo Dominic Etienne Charrier
// order, row, column
// [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1][EXAHYPE_ORDER+1]
extern double*** iK1;

/**
 * \brief Left extrapolation coefficients
 */
// todo Dominic Etienne Charrier
// order, row,
// [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1]
extern double** FLCoeff;

/**
 * \brief Right extrapolation coefficients
 */
// todo Dominic Etienne Charrier
// order, row
// [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];
extern double** FRCoeff;

/**
 * \brief Joint extrapolation coefficients
 *
 * FCoeff = [[FLCoeff];[FRCoeff]]
 */
// todo Dominic Etienne Charrier
// order, left/right, row
// [EXAHYPE_ORDER+1][2][EXAHYPE_ORDER+1];
extern double*** FCoeff;

/**
 * \brief Projects the nodal DoF located at the Gauss-Legendre nodes
 * onto an uniform grid.
 */
// todo Dominic Etienne Charrier
// order, row, column
// [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];
extern double*** subOutputMatrix;
}

#endif
