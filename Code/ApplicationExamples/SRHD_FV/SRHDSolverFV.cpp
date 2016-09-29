#include "SRHDSolverFV.h"

#include <memory>

using std::endl;
using std::cout;

extern "C" {
void hastoadjustsolution_(double* time, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
void pdeflux_(double* F, const double* const Q);
void pdeeigenvalues_(double* lambda, const double* const Q, const double* nv);
}

SRHD::SRHDSolverFV::SRHDSolverFV(int cellsPerCoordinateAxis, double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler):
  exahype::solvers::FiniteVolumesSolver("SRHDSolverFV", 5, 0, cellsPerCoordinateAxis, maximumMeshSize, timeStepping, std::move(profiler)) {
  // @todo Please implement/augment if required
}



void SRHD::SRHDSolverFV::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  if ( tarch::la::equals(t,0.0) ) {
  // @todo Please implement and set initial conditions
    adjustedsolutionvalues_(x, &w, &t, &dt, Q);
  }
  // @todo Feel free to add further conditions
}


bool SRHD::SRHDSolverFV::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t) {
  // @todo Please implement
  if ( tarch::la::equals(t,0.0) ) {
    // Tell kernel that you want to set initial conditions 
    return true;
  }
  else {
    // @todo Please implement
    return false; 
  }
}


exahype::solvers::Solver::RefinementControl SRHD::SRHDSolverFV::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void SRHD::SRHDSolverFV::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // @todo Please implement
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void SRHD::SRHDSolverFV::flux(const double* const Q, double** F) {
  // @todo Please implement
  pdeflux_(F[0], Q);
}


void SRHD::SRHDSolverFV::source(const double* const Q, double* S) {
  S[0] = 0.;
  S[1] = 0.;
  S[2] = 0.;
  S[3] = 0.;
  S[4] = 0.;
}


