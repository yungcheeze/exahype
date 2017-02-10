#include "MHDSolver.h"
//#include "fortran.h" _ltob

#include "InitialDataAdapter.h"
#include "PDE.h"

#include <memory>
#include <cstring>
#include <stdio.h>
#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h" // matrix indexing

/* This is the MHDSolver.cpp binding to Fortran functions, as done in SRHD. */


void SRMHD::MHDSolver::init(std::vector<std::string>& cmdargs){ //, exahype::Parser::ParserView _constants) {
  // just pass the pointer to the crazy Fortran glue code. Should be improved.
  // constants = &_constants;
}

void SRMHD::MHDSolver::flux(const double* const Q, double** F) {
  // Caveats: Fortran accepts a uniform array of size (nVar*nDim), however C passes an array of pointers.
  // This Fortran interface works only if F is a continous array and F[1]==F[nDim+1] etc!
  pdeflux_(F[0], Q);
}



void SRMHD::MHDSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}



bool SRMHD::MHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t, double dt) {
  return (t < 1e-10);
}



void SRMHD::MHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Fortran call:
  adjustedsolutionvalues_(x, &w, &t, &dt, Q);
}

void SRMHD::MHDSolver::source(const double* const Q, double* S) {
  pdesource_(S, Q);
}



exahype::solvers::Solver::RefinementControl SRMHD::MHDSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void SRMHD::MHDSolver::boundaryValues(const double* const x,const double t, const double dt, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut) {
/*  
  // These are the no-boundary conditions:
  for(int i=0; i < SRMHD::MHDSolver::numberOfVariables; i++) {
      fluxOut[i] = fluxIn[i];
      stateOut[i] = stateIn[i];
  }
  return;
  */

  /*
  // previous solution without integration
  double f[nVar], g[nVar], h[nVar];
  double* F[3];
  F[0] = f; F[1] = g; F[2] = h;
  alfenwave_(x, stateOut, &t);
  flux(stateOut, F);
  for(int i=0; i < SRMHD::MHDSolver::numberOfVariables; i++) {
    fluxOut[i] = F[normalNonZero][i];
  }
  */

  const int nVar = SRMHD::AbstractMHDSolver::NumberOfVariables;
  const int order = SRMHD::AbstractMHDSolver::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;

  double Qgp[nVar];
  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut, 0, nVar * sizeof(double));

  double F[3 * nVar]; // Fortran needs continous storage!
                      // Use always 3 dimensions here since the kernels works with those internally; see nDim in PDE.f90;


  kernels::idx2 F_idx(nDim, nVar);

  // Integrate solution in gauss points (Qgp) in time
  if (faceIndex==2 || faceIndex==3) {     
     for(int m=0; m < nVar; m++) {
        stateOut[m] = stateIn[m];
        fluxOut[m]  = fluxIn[m];
     }
     return;
  }

  for(int i=0; i < basisSize; i++)  { // i == time
     const double weight = kernels::gaussLegendreWeights[order][i];
     const double xi = kernels::gaussLegendreNodes[order][i];
     double ti = t + xi * dt;

     alfenwave_(x, Qgp, &ti);
     pdeflux_(F, Qgp);
     for(int m=0; m < nVar; m++) {
  //if(m==checkm) printf("fluxOut[%d] += %.20e\n", m, weight * F[normalNonZero][m]);
        stateOut[m] += weight * Qgp[m];
        fluxOut[m] += weight * F[F_idx(normalNonZero, m)];
     }
  }

  // debugging stuff:
  //printf("stateOut[%d]=%.5e == stateIn[%d]=%.5e at t=%f, dt=%f x=%.5e, y=%.5e, Equals=%e\n", statem, stateOut[statem], statem, stateIn[statem], t, dt, x[0], x[1], stateOut[statem] - stateIn[statem]);
  //printf("FOut[%d]=%e == Fin[%d]=%e\n", statem, fluxOut[statem], statem, fluxIn[statem]);
}


void SRMHD::MHDSolver::ncp(const double* const Q, const double* const gradQ, double* BgradQ) {
  const int nVar = SRMHD::AbstractMHDSolver::NumberOfVariables;
  std::memset(BgradQ, 0, nVar * sizeof(double));
}

void SRMHD::MHDSolver::matrixb(const double* const Q, const int normalNonZero, double* Bn) {
  const int nVar = SRMHD::AbstractMHDSolver::NumberOfVariables;
  std::memset(Bn, 0, nVar * nVar * sizeof(double));
}

bool SRMHD::MHDSolver::isDummyKRequired() const {
  return false;
}
