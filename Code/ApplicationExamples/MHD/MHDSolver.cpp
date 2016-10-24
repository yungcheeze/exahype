#include "MHDSolver.h"
//#include "fortran.h" _ltob

// exact BC
#include "InitialDataAdapter.h"
#include "GeneratedConstants.h"

#include <memory>
#include <cstring>
#include <stdio.h>
#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h" // matrix indexing

/* This is the MHDSolver.cpp binding to Fortran functions, as done in SRHD. */

// Fortran functions:
extern "C" {
void minimumtreedepth_(int* depth);
void hastoadjustsolution_(double* t, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
void pdeflux_(double* F, const double* const Q);
void pdeeigenvalues_(double* lambda, const double* const Q, double* nv);
void registerinitialdata_(const char* const id_name, int* id_name_len);
}/* extern "C" */

// storage for static class members
int MHDSolver::MHDSolver::numberOfVariables;
int MHDSolver::MHDSolver::numberOfParameters;
exahype::Parser::ParserView* MHDSolver::MHDSolver::constants;

void MHDSolver::MHDSolver::init() {
  // This function is called inside the constructur.
  // @todo Please implement/augment if required

  // cf issue #61
  numberOfVariables = getNumberOfVariables();
  numberOfParameters = getNumberOfParameters();

  printf("MHD static: numberOfVariables: %d, numberOfParameters: %d, DIMENSIONS=%d\n", numberOfVariables, numberOfParameters, DIMENSIONS);
}

void MHDSolver::MHDSolver::flux(const double* const Q, double** F) {
  // Caveats: Fortran accepts a uniform array of size (nVar*nDim), however C passes an array of pointers.
  // This Fortran interface works only if F is a continous array and F[1]==F[nDim+1] etc!
  pdeflux_(F[0], Q);
}



void MHDSolver::MHDSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}



bool MHDSolver::MHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  return (t < 1e-10);
}



void MHDSolver::MHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Fortran call:
  adjustedsolutionvalues_(x, &w, &t, &dt, Q);
}

void MHDSolver::MHDSolver::source(const double* const Q, double* S) {
  // TODO: pass this to Fortran.
  for(int i=0; i < MHDSolver::MHDSolver::numberOfVariables; i++) {
    S[i] = 0.0;
  }
}



exahype::solvers::Solver::RefinementControl MHDSolver::MHDSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void MHDSolver::MHDSolver::boundaryValues(const double* const x,const double t, const double dt, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut) {
/*	
  // These are the no-boundary conditions:
  for(int i=0; i < MHDSolver::MHDSolver::numberOfVariables; i++) {
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
  for(int i=0; i < MHDSolver::MHDSolver::numberOfVariables; i++) {
    fluxOut[i] = F[normalNonZero][i];
  }
  */

  // too lazy to access MHDSolver::MHDSolver::numberOfVariables etc.
  const int nVar = MY_NUMBER_OF_VARIABLES;
  const int nDim = MY_DIMENSIONS;
  const int order = MY_POLYNOMIAL_DEGREE;
  const int basisSize = order + 1;

  double Qgp[nVar];
  std::memset(stateOut, 0.0, nVar * sizeof(double));
  std::memset(fluxOut, 0.0, nVar * sizeof(double));

  double F[nDim * nVar]; // Fortran needs continous storage!
  kernels::idx2 F_idx(nDim, nVar);

  // Integrate solution in gauss points (Qgp) in time
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
  const int statem=6; // the interesting component

  // debugging stuff:
  //printf("stateOut[%d]=%.5e == stateIn[%d]=%.5e at t=%f, dt=%f x=%.5e, y=%.5e, Equals=%e\n", statem, stateOut[statem], statem, stateIn[statem], t, dt, x[0], x[1], stateOut[statem] - stateIn[statem]);
  //printf("FOut[%d]=%e == Fin[%d]=%e\n", statem, fluxOut[statem], statem, fluxIn[statem]);
}

