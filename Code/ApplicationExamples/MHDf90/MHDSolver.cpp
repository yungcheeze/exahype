#include "MHDSolver.h"

#include <memory>

// Fortran functions:
extern "C" {
void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
}/* extern "C" */


void MHDSolver::MHDSolver::init(std::vector<std::string>& cmdlineargs) {
  // This function is called by the constructor.
  // You can access spec file parameters as well as command line arguments (argv as std::vector).
  // @todo Please implement/augment if required.
}

//************************************************* 
//for FORTRAN kernels the fluxes and eigenvalues 
//have to be implemented in the file ./PDE.f90 and ./typesDef.f90 
//
//You have to create these files yourself
//and follow the sample applications in the wiki
//************************************************* 






void MHDSolver::MHDSolver::boundaryValues(const double* const x,const double t, const double dt, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)


  // @todo Please implement
  // fluxOut
  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];
  fluxOut[3] = fluxIn[3];
  fluxOut[4] = fluxIn[4];
  fluxOut[5] = fluxIn[5];
  fluxOut[6] = fluxIn[6];
  fluxOut[7] = fluxIn[7];
  fluxOut[8] = fluxIn[8];
  // stateOut
  // @todo Please implement
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];
  stateOut[5] = stateIn[5];
  stateOut[6] = stateIn[6];
  stateOut[7] = stateIn[7];
  stateOut[8] = stateIn[8];
}






bool MHDSolver::MHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  // @todo Please implement
  return (t < 1e-10);
}



void MHDSolver::MHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)
  // @todo Please implement

  adjustedsolutionvalues_(x, &w, &t, &dt, Q);
}



exahype::solvers::Solver::RefinementControl MHDSolver::MHDSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}









