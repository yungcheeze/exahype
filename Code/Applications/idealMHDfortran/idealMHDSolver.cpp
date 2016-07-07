#include "idealMHDSolver.h"

#include <memory>

srhd3dfortran::idealMHDSolver::idealMHDSolver( int kernelNumber, std::unique_ptr<exahype::profilers::Profiler> profiler):
  exahype::solvers::Solver("idealMHDSolver",exahype::solvers::Solver::ADER_DG,kernelNumber,9,3+1,exahype::solvers::Solver::GlobalTimeStepping, std::move(profiler)) {
  // @todo Please implement/augment if required
}



int srhd3dfortran::idealMHDSolver::getMinimumTreeDepth() const {
  // @todo Please implement
  return 3;
}



bool srhd3dfortran::idealMHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  // @todo Please implement
  return false;
}



void srhd3dfortran::idealMHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
	// Dimensions             = 3
	// Number of variables    = 9 (#unknowns + #parameters)

	if (tarch::la::equals(t, 0.0)) {
		double V[9];
		
		// specifiy initial values for V[0]...V[8] here.
		
		// prim2cons(V, Q);
	}
	
	// why does the toolkit propose this?:

	Q[0] = 0.0;
	Q[1] = 0.0;
	Q[2] = 0.0;
	Q[3] = 0.0;
	Q[4] = 0.0;
	Q[5] = 0.0;
	Q[6] = 0.0;
	Q[7] = 0.0;
	Q[8] = 0.0;
}



exahype::solvers::Solver::RefinementControl srhd3dfortran::idealMHDSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}






//************************************************* 
//for FORTRAN kernels the fluxes and eigenvalues 
//have to be implemented in the file ./PDE.f90 
//************************************************* 



