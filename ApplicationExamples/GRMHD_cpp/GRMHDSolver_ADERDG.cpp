#include "GRMHDSolver_ADERDG.h"
#include <algorithm> // fill_n
#include <cstring> // memset
#include <fenv.h> // enable nan tracker

#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

#include "peano/utils/Dimensions.h" // Defines DIMENSIONS
//namespace GRMHD { constexpr int nVar = GRMHDSolver_FV::NumberOfVariables; } // ensure this is 19 or so
#include "PDE/PDE.h"
#include "InitialData.h"

#include "GRMHDSolver_ADERDG_Variables.h"


constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
constexpr int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
constexpr int basisSize = order + 1;
constexpr int nDim = DIMENSIONS;

tarch::logging::Log GRMHD::GRMHDSolver_ADERDG::_log( "GRMHD::GRMHDSolver_ADERDG" );

constexpr double magicCheck = 123456789;

void GRMHD::GRMHDSolver_ADERDG::init(std::vector<std::string>& cmdlineargs) { // ,  exahype::Parser::ParserView constants) {

	// feenableexcept(FE_INVALID | FE_OVERFLOW);  // Enable all floating point exceptions but FE_INEXACT
	
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue GRMHD::GRMHDSolver_ADERDG::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
	// @todo Please implement/augment if required
	return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PointWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;
}

void zeroHelpers(double* Q) {
	// debugging variables
	GRMHD::AbstractGRMHDSolver_ADERDG::Variables var(Q);
	DFOR(i) var.pos(i) = 0;
	var.check() = 0;
}

void GRMHD::GRMHDSolver_ADERDG::adjustPointSolution(const double* const x,const double w,const double t,const double dt,double* Q) {
	// Number of variables    = 23 + #parameters
	
	// currently, the C++ AlfenWave spills out primitive data
	double V[nVar];
	AlfenWave id(x,t,V);
	GRMHD::Prim2Cons(Q, V).copyFullStateVector();
	
	// also store the positions for debugging
	GRMHD::AbstractGRMHDSolver_ADERDG::Variables var(Q);
	DFOR(i) var.pos(i) = x[i];
	var.check() = magicCheck;
	

}

void GRMHD::GRMHDSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
	// Provide NVARS eigenvalues
	PDE::eigenvalues(Q, d, lambda);
}


void GRMHD::GRMHDSolver_ADERDG::flux(const double* const Q,double** F) {
	GRMHD::Fluxes(F, Q).zeroMaterialFluxes();
	DFOR(d) zeroHelpers(F[d]);
}


void GRMHD::GRMHDSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int d,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
	// for debugging, to make sure BC are set correctly
	/*
	double snan = std::numeric_limits<double>::signaling_NaN();
	double weird_number = -1.234567;
	std::memset(stateOut, weird_number, nVar * sizeof(double));
	std::memset(fluxOut,  weird_number, nVar * sizeof(double));
	*/
	
	// employ time-integrated exact BC for AlfenWave.

	double Qgp[nVar], Fs[nDim][nVar], *F[nDim];
	for(int dd=0; dd<nDim; dd++) F[dd] = Fs[dd];
	// zeroise stateOut, fluxOut
	for(int m=0; m<nVar; m++) {
		stateOut[m] = 0;
		fluxOut[m] = 0;
	}
	for(int i=0; i < basisSize; i++)  { // i == time
		const double weight = kernels::gaussLegendreWeights[order][i];
		const double xi = kernels::gaussLegendreNodes[order][i];
		double ti = t + xi * dt;

		adjustPointSolution(x, weight/*not sure, not used anyway*/, ti, dt, Qgp);
		flux(Qgp, F);

		for(int m=0; m < nVar; m++) {
			stateOut[m] += weight * Qgp[m];
			fluxOut[m] += weight * Fs[d][m];
		}
	}
}


exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void GRMHD::GRMHDSolver_ADERDG::fusedSource(const double* const Q, const double* const gradQ, double* S_) {
	/*
	// ExaHyPE workaround:
	// if the input are zeros everywhere, complain
	if(isUnphysical(Q)) {
		//printf("WRONG FusedSource input, all zero!\n");
		//NVARS(m) printf("Q[%d]=%e\n", m, Q[m]);
		//std::abort();
		// Set everything to NAN
		NVARS(m) S_[m] = NAN;
		return;
	}
	*/
	PDE::Source S(S_);
	PDE pde(Q);
	GRMHD::AbstractGRMHDSolver_ADERDG::ReadOnlyVariables var(Q);
/*
	constexpr double eps = 1e-8;
	// ExaHyPE workaround:
	// if the input has zeros at weird places, don't do anything
	if(std::abs(pde.gam.det) < eps || var.check()!=magicCheck) {
		printf("Weird FusedSource input (det=%e), skipping. ",pde.gam.det);
		//SYMFOR(i,j) printf("gam(%d,%d)=%e\n", i,j, pde.gam.lo(i,j));
		std::cout << "pos: " << var.pos() << " check: " << var.check() <<  std::endl;
		// set everything to NANs
		NVARS(m) S_[m] = NAN;
		return;
	} else {
		//printf("Good FusedSource input (det=%e)\n",pde.gam.det);
		//NVARS(m) printf("FusedSource input: Q[%d]=%e\n", m, Q[m]);
		
		
		//std::cout << "Good FusedSource, x= " << var.pos() <<  std::endl;
	}
*/
	pde.RightHandSide(gradQ, S);
	S.zero_adm();
	zeroHelpers(S_);
}


void GRMHD::GRMHDSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
	PDE::NCP ncp(BgradQ);
	const Gradients g(gradQ);
	PDE(Q).nonConservativeProduct(g, ncp);
	ncp.zero_adm();
	zeroHelpers(BgradQ);
}
