#ifndef _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_
#define _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_

//#define GAMMA 1.4

#include "tarch/la/Vector.h"
#include "peano/utils/Globals.h"

namespace kernels {
  namespace aderdg {
    namespace generic {
      template <void PDEFlux2d(const double * const Q, double * f,double * g)>
      void spaceTimePredictor(
        double * lQi,
        double * lFi,
        const double * const luh,
        double * lQhi,
        double * lFhi,
        double * lQhbnd,
        double * lFhbnd,
        const tarch::la::Vector<DIMENSIONS,double>&  dx,
        const double dt,
        int          order,
        int          numberOfVariables
      );


      void solutionUpdate(double * luh, const double * const lduh, const tarch::la::Vector<DIMENSIONS,double>&  dx, const double dt);
      void volumeIntegral(double * lduh, const double * const lFhi, const tarch::la::Vector<DIMENSIONS,double>&  dx);
      void surfaceIntegral(double * lduh, const double * const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>&  dx);

      template <void PDEFlux2d(const double * const Q, double * f,double * g)>
      void riemannSolver(double * FL, double * FR, const double * const QL, const double * const QR, const double dt, const double hFace, const tarch::la::Vector<DIMENSIONS,double>&  n);

      template <void PDEEigenvalues2d(const double * const Q, const tarch::la::Vector<DIMENSIONS,double>& n, const int d, double * lambda)>
      double stableTimeStepSize(const double * const luh, const tarch::la::Vector<DIMENSIONS,double>&  dx);
    }
  }
}


#include "kernels/aderdg/generic/Kernels.cpph"

#endif
