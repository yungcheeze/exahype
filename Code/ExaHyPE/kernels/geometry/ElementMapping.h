///////////////////////////////////////////////////////////////////////////////
//                    HEADER FILE FOR ExaHyPE CODE                              //
///////////////////////////////////////////////////////////////////////////////
//            D.E. Charrier                                                  //
// Code:      ExaHyPe                                                           //
// File:      dghelpers.h                                                    //
///////////////////////////////////////////////////////////////////////////////
/** \file ElementGeometry.h
 *  \brief Defines mappings from the reference element to physical mesh elements.
 *
 *  @todo DEC: The functions should be modified such that they take vector arguments and return vectors.
 */
///////////////////////////////////////////////////////////////////////////////

#ifndef EXAHYPE_GEOMETRY_ELEMENT_MAPPING_H_
#define EXAHYPE_GEOMETRY_ELEMENT_MAPPING_H_

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
        const double dx,const double dy,
        const double r,const double s,
        double* x,double* y) {
      *x=centerX + dx*(r-0.5);
      *y=centerY + dy*(s-0.5);
    }

    /**
     * @brief: Maps reference coordinates to physical coordinates.
     *
     * @param[in] centerX x coordinate of the center of the cell
     * @param[in] centerY y coordinate of the center of the cell
     * @param[in] centerZ z coordinate of the center of the cell
     * @param[in] r reference coordinate in [0,1]
     * @param[in] s reference coordinate in [0,1]
     * @param[in] t reference coordinate in [0,1]
     * @param[out] x physical coordinate
     * @param[out] y physical coordinate
     * @param[out] z physical coordinate
     */
    inline void
    mapping3d(
        const double centerX,
        const double centerY,
        const double centerZ,
        const double dx,
        const double dy,
        const double dz,
        const double r,
        const double s,
        const double t,
        double* x,
        double* y,
        double* z
) {
      *x=centerX + dx*(r-0.5);
      *y=centerY + dy*(s-0.5);
      *z=centerZ + dz*(t-0.5);
    }
  }
}

#endif /* EXAHYPE_GEOMETRY_ELEMENT_MAPPING_H_ */
