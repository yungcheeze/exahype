#include "SRHDSolver.h"

#include "fortran-interface.h"

#include <memory>

srhdhybrid::SRHDSolver::SRHDSolver( int kernelNumber, std::unique_ptr<exahype::profilers::Profiler> profiler):
  exahype::solvers::Solver("SRHDSolver",exahype::solvers::Solver::ADER_DG,kernelNumber,5,1+1,exahype::solvers::Solver::GlobalTimeStepping, std::move(profiler)) {
  // @todo Please implement/augment if required
}



int srhdhybrid::SRHDSolver::getMinimumTreeDepth() const {
  // @todo Please implement
  return 3;
}



void srhdhybrid::SRHDSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)

  // I would love to see these numbers as variables provided by the Toolkit:
  int dimensions = 2, nVar = 5;

  int flength = nVar * dimensions;
  int qlength = nVar;

  // F is a multimendional array, however the semantics of
  // f = F[0], g = F[1], h = F[2]
  // are the same as in Fortran's F(-,1), F(-,2), F(-,3)
  // That is, no conversion is neccessary at this point

  PDEFlux_(F, (fdvec) Q, flength, qlength);
}



void srhdhybrid::SRHDSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)

  // I would love to see these numbers as variables provided by the Toolkit:
  int dimensions = 2, nVar = 5;

  int lambdalength = nVar;
  int qlength = nVar;
  int nvlength = dimensions;

  // construct the nv vector as given in Fortran kernels
  double nv[nvlength] = { 0. };
  nv[normalNonZeroIndex-1] = 1.

  PDEEigenvalues_(lambda, (fdvec) Q, nv, lambdalength, qlength, nvlength);
}



bool srhdhybrid::SRHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  // @todo Please implement
  return false;
}



void srhdhybrid::SRHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  // @todo Please implement
  Q[0] = 0.0;
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
}



exahype::solvers::Solver::RefinementControl srhdhybrid::SRHDSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::Keep;
}



