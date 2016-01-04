///////////////////////////////////////////////////////////////////////////////
//                    HEADER FILE FOR ExaHyPE CODE                              //
///////////////////////////////////////////////////////////////////////////////
//            D.E. Charrier                                                  //
// Code:      ExaHyPe                                                           //
// File:      dghelpers.h                                                    //
///////////////////////////////////////////////////////////////////////////////
/** \file dghelpers.h
 *  \brief Defines mappings from the reference element to physical mesh elements.
 *
 *  -DEC: The functions should be modified such that they take vector arguments and return vectors.
 */
///////////////////////////////////////////////////////////////////////////////

#ifndef EXAHYPE_GEOMETRY_MAPPING_H_
#define EXAHYPE_GEOMETRY_MAPPING_H_

namespace exahype {
  namespace geometry {
    /**
     * @brief: Maps reference coordinates to physical coordinates.
     *
     * @param[in] centerX x coordinate of the center of the cell
     * @param[in] centerY y coordinate of the center of the cell
     * @param[in] r reference coordinate in [0,1]
     * @param[in] s reference coordinate in [0,1]
     * @param[out] x physical coordinate
     * @param[out] y physical coordinate
     */
    inline void
    mapping2d(
        const double centerX,const double centerY,
        const double dx,
        const double dy,
        const double r,const double s,
        double* x,double* y) {
      *x=centerX + 0.5*dx*(r-1.0);
      *y=centerY + 0.5*dy*(s-1.0);
    }

    /**
     * @brief: Maps reference coordinates to physical coordinates.
     *
     * @param[in] centerX x coordinate of the center of the cell
     * @param[in] centerY y coordinate of the center of the cell
     * @param[in] dx extent of the cell in x direction
     * @param[in] dy extent of the cell in y direction
     * @param[in] dxPatch extent of the patch cell in x direction
     * @param[in] dyPatch extent of the patch cell in y directio
     * @param[in] i index of the patch cell
     * @param[in] j index of the patch cell
     * @param[in] r reference coordinate in [0,1]
     * @param[in] s reference coordinate in [0,1]
     * @param[out] x physical coordinate
     * @param[out] y physical coordinate
     */
    inline void
    mapping2d(
        const double centerX,const double centerY,
        const double dx,const double dy,
        const double dxPatch,const double dyPatch,
        const int i,const int j,
        const double r,const double s,
        double* x,double* y) {
      *x=centerX - 0.5*dx + (i-1) * dxPatch + dxPatch * r;
      *y=centerY - 0.5*dy + (j-1) * dyPatch + dyPatch * s;
    }


    /**
     * @brief: Maps reference coordinates to physical coordinates.
     *
     * @param[in] centerX x coordinate of the center of the cell
     * @param[in] centerY y coordinate of the center of the cell
     * @param[in] r reference coordinate in [0,1]
     * @param[in] s reference coordinate in [0,1]
     * @param[out] x physical coordinate
     * @param[out] y physical coordinate
     */
    inline void
    mapping3d(
        const double centerX,const double centerY,
        const double dx,
        const double dy,
        const double r,const double s,
        double* x,double* y) {
      *x=centerX + 0.5*dx*(r-1.0);
      *y=centerY + 0.5*dy*(s-1.0);
      // todo
    }
  }
}

#endif /* EXAHYPE_GEOMETRY_MAPPING_H_ */
