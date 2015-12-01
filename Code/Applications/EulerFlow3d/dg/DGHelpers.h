///////////////////////////////////////////////////////////////////////////////
//                    HEADER FILE FOR ExaHyPE CODE                              //
///////////////////////////////////////////////////////////////////////////////
//            D.E. Charrier                                                  //
// Code:      ExaHyPe                                                           //
// File:      dghelpers.h                                                    //
///////////////////////////////////////////////////////////////////////////////
/** \file dghelpers.h
 *  \brief Defines helper functions for DG matrix assembling.
 *
 *  -DEC: The quadrature lookup functions should be modified such that they return vectors.
 */
///////////////////////////////////////////////////////////////////////////////
#ifndef EXAHYPE_DG_DGHELPERS_H_
#define EXAHYPE_DG_DGHELPERS_H_

#include "EulerFlow3d/dg/Constants.h"

#include "EulerFlow3d/math/quad/Gausslegendre.h"

#define DG_ONE 1

namespace exahype {
  namespace dg {
    /**
     * Outward directed normal vectors on each face of the current cell.
     */
    ///@{
    extern const double normal[3][6];
    ///@}


    ///////////////////////////////////////////////////////////////////////////////
    /**
     * \brief Returns the neighbor face id to a given.
     */
    inline void GetNeighborFace(const int face,int *face_neighbor) {
      // pure integer arithmetic: *face_neighbor = face + ((face+1) % 2) - ((face) % 2);
      switch (face) {
      case EXAHYPE_FACE_LEFT:
      case EXAHYPE_FACE_FRONT:
      case EXAHYPE_FACE_BOTTOM:
        *face_neighbor = face+1;
        break;
      case EXAHYPE_FACE_RIGHT:
      case EXAHYPE_FACE_BACK:
      case EXAHYPE_FACE_TOP:
        *face_neighbor = face-1;
        break;
      default:
        break;
      }
    }

    /**
     * \brief Returns a quadrature node value.
     */
    inline void GetFaceQr(const int n_q_points,const int iq,const int face,double *qr) {
      switch (face) {
      case EXAHYPE_FACE_LEFT:
        *qr=-DG_ONE;
        break;
      case EXAHYPE_FACE_RIGHT:
        *qr=+DG_ONE;
        break;
      default:
        *qr=exahype::quad::gaussLegendreNodes[n_q_points-1][iq];
        break;
      }
    }

    /**
     * \brief Returns a quadrature node value.
     */
    inline void GetFaceQs(const int n_q_points,const int iq,const int face,double *qs) {
      switch (face) {
      case EXAHYPE_FACE_FRONT:
        *qs=-DG_ONE;
        break;
      case EXAHYPE_FACE_BACK:
        *qs=+DG_ONE;
        break;
      default:
        *qs=exahype::quad::gaussLegendreNodes[n_q_points-1][iq];
        break;
      }
    }

    /**
     * \brief Returns a quadrature node value.
     */
    inline void GetFaceQt(const int n_q_points,const int iq,const int face,double *qt) {
      switch (face) {
      case EXAHYPE_FACE_BOTTOM:
        *qt=-DG_ONE;
        break;
      case EXAHYPE_FACE_TOP:
        *qt=+DG_ONE;
        break;
      default:
        *qt=exahype::quad::gaussLegendreNodes[n_q_points-1][iq];
        break;
      }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
#endif /* EXAHYPE_DG_DGHELPERS_H_ */
