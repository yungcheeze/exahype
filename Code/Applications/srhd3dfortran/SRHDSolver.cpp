#include "SRHDSolver.h"



srhd3dfortran::SRHDSolver::SRHDSolver(const std::string& identifier, exahype::solvers::Solver::Type type, int kernelNumber, int numberOfVariables, int numberOfParameters, int nodesPerCoordinateAxis, tarch::la::Vector<DIMENSIONS,double> maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler)
  : exahype::solvers::Solver(
            identifier, type, kernelNumber, numberOfVariables, numberOfParameters, nodesPerCoordinateAxis, maximumMeshSize, timeStepping, std::move(profiler)) {
  // @todo Please implement/augment if required
}



int srhd3dfortran::SRHDSolver::getMinimumTreeDepth() const {
  // @todo Please implement
  return 3;
}



bool srhd3dfortran::SRHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  // @todo Please implement
  if (tarch::la::equals(t, 0.0)) {
    return true;
  }
  return false;

}

// This is a port of PDE_prim2con from PDE.f90 to C.
// It would be much nicer to be able to call the Fortran routine directly.
bool PDEC_prim2con(double* Q, const double* const V) {
	// Make sure this gamma is the same as in typesDef.f90!
	double const gamma = 5./3.;
	
	double rho=V[0], vx=V[1], vy=V[2], vz=V[3], p=V[4];
	double v2 = vx*vx + vy*vy + vz*vz;
	
	if(v2 > 1.0) {
		printf("Superluminal velocity in PDEC_prim2con!");
		return false;
	}
	
	double lf     =  1.0 / sqrt(1.0 - v2);
	double gamma1 = gamma/(gamma-1.0);
	double w      = rho + gamma1*p;
	double ww     = w*lf * w*lf;

	Q[0]   = rho*lf;
	Q[1]   = ww*vx;
	Q[2]   = ww*vy;
	Q[3]   = ww*vz;
	Q[4]   = ww - p - Q[1];

	return true;
}



void srhd3dfortran::SRHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
	// Dimensions             = 3
	// Number of variables    = 5 (#unknowns + #parameters)
	// @todo Please implement

	/* Problem 1 of Table 1 of Zanotti, Dumbser j.cpc.2014.11.015 */
	const double rho_L = 1, v_L = -0.6, p_L = 10, rho_R = 10, v_R = 0.5, p_R = 20, gamma = 5./3., t_f = 0.4, xsep=0.5;

	/* Riemann problem, 1D. */
	if (tarch::la::equals(t, 0.0, 1e-15)) {
		// How can I access nVar from typesDef.f90?
		const int nVar = 5;
		double V[nVar];
		
		if(x[0] < xsep) {
			V[0] = rho_L;
			V[1] = v_L;
			V[2] = 0.;
			V[3] = 0.;
			V[4] = p_L;
		} else {
			V[0] = rho_R;
			V[1] = v_R;
			V[2] = 0.;
			V[3] = 0.;
			V[4] = p_R;
		}
		
		PDEC_prim2con(Q, V);
	}
}

exahype::solvers::Solver::RefinementControl srhd3dfortran::SRHDSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
  // Problem here: Determine how to access conserved quantities Q
  // which give the criterion whether to refine or not.
}






//************************************************* 
//for FORTRAN kernels the fluxes and eigenvalues 
//have to be implemented in the file ./PDE.f90 
//************************************************* 



