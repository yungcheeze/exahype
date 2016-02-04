#ifndef _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_
#define _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_

#include "tarch/la/Vector.h"
#include "peano/utils/Globals.h"

namespace kernels {
  namespace aderdg {
    namespace generic {

      // @todo Dominic Etienne Charrier
      // Inconsistent ordering of inout and in arguments for
      // template argument functions and non-template argument function.
      template <void PDEFlux(const double* const Q,double* f,double* g)>
      void spaceTimePredictor(
          double* lQi,
          double* lFi,
          double* lQhi,
          double* lFhi,
          double* lQhbnd,
          double* lFhbnd,
          const double* const luh,
          const tarch::la::Vector<DIMENSIONS,double>&  dx,
          const double dt,
          int          numberOfVariables,
          int          basisSize
      );


      // todo Dominic Etienne Charrier:
      // The DIMENSIONS depending mesh size vector enables overloading at the moment.
      // If we replace it by scalar mesh size, we have to add a template argument "int dim".

      void solutionUpdate(
          double* luh,
          const double* const lduh,
          const double dt,
          const int numberOfVariables,
          const int basisSize
      );

      void volumeIntegral(
          double* lduh,
          const double* const lFhi,
          const tarch::la::Vector<DIMENSIONS,double>&  dx,
          const int numberOfVariables,
          const int basisSize

      );
      void surfaceIntegral(
          double* lduh,
          const double* const lFhbnd,
          const tarch::la::Vector<DIMENSIONS,double>&  dx,
          const int numberOfVariables,
          const int basisSize
      );

      // @todo Dominic Etienne Charrier
      // Inconsistent ordering of inout and in arguments for
      // template argument functions and non-template argument function.
      template <void PDEInitialValues(const double* const x,double* Q)>
      void initialCondition(
          double* luh,
          const tarch::la::Vector<DIMENSIONS,double>& center,
          const tarch::la::Vector<DIMENSIONS,double>& dx,
          const int numberOfVariables,
          const int basisSize
      );

      // @todo Dominic Etienne Charrier
      // Inconsistent ordering of inout and in arguments
      // template argument functions and non-template argument function.
      template <void PDEEigenvalues(const double* const Q,const int normalNonZero,double* lambda)>
      void riemannSolver(
          double* FL,
          double* FR,
          const double* const QL,
          const double* const QR,
          const double dt,
          const int normalNonZero,
          const int numberOfVariables,
          const int basisSize
      );

      // @todo Dominic Etienne Charrier
      // Inconsistent ordering of inout and in arguments for
      // template argument functions and non-template argument function.
      template <void PDEEigenvalues(const double* const Q,const int normalNonZero,double* lambda)>
      double stableTimeStepSize(
          const double* const luh,
          const tarch::la::Vector<DIMENSIONS,double>& dx,
          const int numberOfVariables,
          const int basisSize
      );
    }
  }
}


#include "kernels/aderdg/generic/Kernels.cpph"

#endif
