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
     * @param x0 x coordinate of the center of the cell
     * @param y0 y coordinate of the center of the cell
     * @param r reference coordinate in [-1,1]
     * @param s reference coordinate in [-1,1]
     * @param x physical coordinate
     * @param y physical coordinate
     */
    inline void
    mapping2d(
        const double x0,const double y0,
        const double h,
        const double r,const double s,
        double* x,double* y) {
      *x=x0+.5*h*r;
      *y=y0+.5*h*s;
    }

    /**
     * @brief: Maps reference coordinates to physical coordinates.
     *
     * @param x0 x coordinate of the center of the cell
     * @param y0 y coordinate of the center of the cell
     * @param r reference coordinate
     * @param s reference coordinate
     * @param x physical coordinate
     * @param y physical coordinate
     */
    inline void
    mapping3d(
        const double x0,const double y0,const double z0,
        const double h,
        const double r,const double s,const double t,
        double* x,double* y,double *z) {
      *x=x0+.5*h*r;
      *y=y0+.5*h*s;
      *z=z0+.5*h*t;
    }
  }
}

#endif /* EXAHYPE_GEOMETRY_MAPPING_H_ */
