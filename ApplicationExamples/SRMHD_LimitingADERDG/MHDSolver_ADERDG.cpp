#include "MHDSolver_ADERDG.h"
//#include "fortran.h" _ltob

#include "PDE.h"

#include <memory>
#include <cstring>
#include <stdio.h>
#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h" // matrix indexing

/* This is the MHDSolver_ADERDG.cpp binding to Fortran functions, as done in SRHD. */


void MHD::MHDSolver_ADERDG::init(std::vector<std::string>& cmdargs, exahype::Parser::ParserView constants) {
  // do nothing
}

void MHD::MHDSolver_ADERDG::flux(const double* const Q, double** F) {
  // Caveats: Fortran accepts a uniform array of size (nVar*nDim), however C passes an array of pointers.
  // This Fortran interface works only if F is a continous array and F[1]==F[nDim+1] etc!
  pdeflux_(F[0], Q);
}

void MHD::MHDSolver_ADERDG::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

bool MHD::MHDSolver_ADERDG::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t, double dt) {
  return (t < 1e-10);
}

void MHD::MHDSolver_ADERDG::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Fortran call:
  adjustedsolutionvalues_(x, &w, &t, &dt, Q);
}

void MHD::MHDSolver_ADERDG::source(const double* const Q, double* S) {
  pdesource_(S, Q);
}

exahype::solvers::Solver::RefinementControl MHD::MHDSolver_ADERDG::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}

/* This is work for the Alfven wave
void MHD::MHDSolver_ADERDG::boundaryValues(const double* const x,const double t, const double dt, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut) {
  // These are the no-boundary conditions:
  constexpr int nVar = 9;
  
  for(int i=0; i < nVar; i++) {
      fluxOut[i]  = fluxIn[i];
      stateOut[i] = stateIn[i];
  }
}
*/

void MHD::MHDSolver_ADERDG::boundaryValues(const double* const x,const double t, const double dt, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut) {
  // These are the no-boundary conditions:
  constexpr int nVar = 9;

  if (faceIndex==1){
    injectjet_(x, stateOut);
                             // Here we need commpute the fluxIn in the circular region
  }
  else{
  for(int i=0; i < nVar; i++) {
      fluxOut[i]  = fluxIn[i];
      stateOut[i] = stateIn[i];
  }
  }
  
}





void MHD::MHDSolver_ADERDG::ncp(const double* const Q, const double* const gradQ, double* BgradQ) {
  constexpr int nVar = 9;
  std::memset(BgradQ, 0, nVar * sizeof(double));
}

void MHD::MHDSolver_ADERDG::matrixb(const double* const Q, const int normalNonZero, double* Bn) {
  constexpr int nVar = 9;
  std::memset(Bn, 0, nVar * nVar * sizeof(double));
}

bool MHD::MHDSolver_ADERDG::physicalAdmissibilityDetection(const double* const QMin,const double* const QMax) {
  if (QMin[0] < 0.0) return false;
  if (QMin[4] < 0.0) return false;

  for (int i=0; i<5; ++i) {
    if (!std::isfinite(QMin[i])) return false;
    if (!std::isfinite(QMax[i])) return false;
  }

  return true;
}

bool MHD::MHDSolver_ADERDG::isDummyKRequired() const {
  return false;
}

void MHD::MHDSolver_ADERDG::dummyK_Value(const double* const x,const double t,const double dt, double* forceVector, double* x0) {
  // do nothing
}
