#include "MHDSolver.h"
//#include "fortran.h" _ltob

#include <memory>

/* This is the MHDSolver.cpp binding to Fortran functions, as done in SRHD. */

// Fortran functions:
extern "C" {
void minimumtreedepth_(int* depth);
void hastoadjustsolution_(double* t, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
void pdeflux_(double* F, const double* const Q);
void pdeeigenvalues_(double* lambda, const double* const Q, double* nv);
}

int MHDSolver::MHDSolver::numberOfVariables = 0;
int MHDSolver::MHDSolver::numberOfParameters = 0;

void MHDSolver::MHDSolver::init() {
  // This function is called inside the constructur.
  // @todo Please implement/augment if required

  // cf issue #61
  numberOfVariables = getNumberOfVariables();
  numberOfParameters = getNumberOfParameters();
}

void MHDSolver::MHDSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)
  // Fortran call
  pdeflux_(F[0], Q);
}



void MHDSolver::MHDSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}



bool MHDSolver::MHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  return (t < 1e-10);

  // Fixed the following.
  bool refine;
  hastoadjustsolution_(&t, &refine);
  //return _ltob(refine); // something like this is needed
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

void MHDSolver::MHDSolver::boundaryValues(const double* const x,const double t, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut) {
  // TODO: Pass this to Fortran
  for(int i=0; i < MHDSolver::MHDSolver::numberOfVariables; i++) {
      fluxOut[i] = fluxIn[i];
      stateOut[i] = stateIn[i];
  }
}

