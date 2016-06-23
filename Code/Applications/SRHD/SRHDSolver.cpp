#include "SRHDSolver.h"

#include <memory>

using std::endl;
using std::cout;

extern "C" {
void minimumtreedepth_(int* depth);
void hastoadjustsolution_(double* time, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
void flux_(const double* const Q, double* F);
void eigenvalues_(const double* const Q, const int* normalNonZeroIndex, double* lambda);
}

SRHD::SRHDSolver::SRHDSolver( int kernelNumber, std::unique_ptr<exahype::profilers::Profiler> profiler):
  exahype::solvers::Solver("SRHDSolver",exahype::solvers::Solver::ADER_DG,kernelNumber,5,1+1,exahype::solvers::Solver::GlobalTimeStepping, std::move(profiler)) {
}



int SRHD::SRHDSolver::getMinimumTreeDepth() const {
  int depth;

  minimumtreedepth_(&depth);
  
  return depth;
}



void SRHD::SRHDSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)

  flux_(Q, F[0]);
  
}



void SRHD::SRHDSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  
  const int nnzi = normalNonZeroIndex+1;
  eigenvalues_(Q, &nnzi, lambda);
  
}



bool SRHD::SRHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  
  bool refine;
  
  hastoadjustsolution_(&t, &refine);
  
  return refine;
  
}



void SRHD::SRHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  
  adjustedsolutionvalues_(x, &w, &t, &dt, Q);
}



exahype::solvers::Solver::RefinementControl SRHD::SRHDSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {

  return exahype::solvers::Solver::Keep;
}



