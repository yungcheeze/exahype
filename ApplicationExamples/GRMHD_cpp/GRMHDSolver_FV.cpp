#include "GRMHDSolver_FV.h"
#include <algorithm> // fill_n

#include "peano/utils/Dimensions.h" // Defines DIMENSIONS
//namespace GRMHD { constexpr int nVar = GRMHDSolver_FV::NumberOfVariables; } // ensure this is 19 or so
#include "PDE/PDE.h"

#include "InitialData.h"
#include "GRMHDSolver_FV_Variables.h"

constexpr int nDim = DIMENSIONS;

using namespace GRMHD;

tarch::logging::Log GRMHD::GRMHDSolver_FV::_log( "GRMHD::GRMHDSolver_FV" );

void GRMHD::GRMHDSolver_FV::init(std::vector<std::string>& cmdlineargs, exahype::Parser::ParserView& constants) {
}

bool GRMHD::GRMHDSolver_FV::useAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) const {
  return tarch::la::equals(t,0.0);
}

constexpr double magicCheck = 123456789;

void GRMHD::GRMHDSolver_FV::adjustSolution(const double* const x,const double w,const double t,const double dt, double* Q) {
	// Set the 9 SRMHD variables (D,S_j,tau,B^j) and the 10 [11] ADM material parameters (N^i,g_ij,[detg])
	
	// currently, the C++ AlfenWave spills out primitive data
	double V[nVar];
	AlfenWave id(x,t,V);
	GRMHD::Prim2Cons(Q, V).copyFullStateVector();
	
	// also store the positions for debugging
	GRMHD::AbstractGRMHDSolver_FV::Variables var(Q);
	DFOR(i) var.pos(i) = x[i];
	var.check() = magicCheck;
	
	// NVARS(m) printf("Qid[%d]=%e\n", m, Q[m]);
	// std::cout << "GRMHDSolver_FV::adjustSolution(" << var.pos() << ")\n";

}

exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver_FV::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
	return exahype::solvers::Solver::RefinementControl::Keep;
}


void GRMHD::GRMHDSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
	// Provide eigenvalues for the 9 SRMHD variables (D,S_j,tau,B^j),
	// we split off the 11 ADM material parameters (N^i,g_ij,detg)
	PDE::eigenvalues(Q, dIndex, lambda);
	//NVARS(m) printf("EV[%d]=%f\n", m, lambda[m]);
}

// Detection of unphysical states. In these cases, the user PDE functions shall never be called.
// We workaround by returning some kind of "neutral" values which go well with the scheme.
bool isUnphysical(const double* const Q) {
	bool allzero=true; NVARS(i) { if(Q[i]!=0) allzero=false; }
	return allzero;
}

void GRMHD::GRMHDSolver_FV::flux(const double* const Q, double** F) {
	// Provide fluxes for the 9 SRMHD variables (D,S_j,tau,B^j),
	// we split off the 11 ADM material parameters (N^i,g_ij,detg)
	//PDE(Q).flux(F);
	
	// zero detection
	if(isUnphysical(Q)) {
		//printf("WRONG Flux input, all zero!\n");
		//NVARS(m) printf("Q[%d]=%e\n", m, Q[m]);
		//std::abort();
		// Set everything to some neutral value
		DFOR(i) NVARS(m) F[i][m] = 0;
		return;
	}
	
	GRMHD::Fluxes(F, Q).zeroMaterialFluxes();
}



void GRMHD::GRMHDSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateIn,
    double* stateOut) {
	// SET BV for all 9 variables + 11 parameters.
	// EXACT FV AlfenWave BC
	adjustSolution(x, 0/*not sure, not used anyway*/, t, dt, stateOut);
}

void GRMHD::GRMHDSolver_FV::fusedSource(const double* const Q, const double* const gradQ, double* S_) {
	
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
	
	PDE::Source S(S_);
	PDE pde(Q);
	GRMHD::AbstractGRMHDSolver_FV::ReadOnlyVariables var(Q);

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
	
	pde.RightHandSide(gradQ, S);
	S.zero_adm();
}

void GRMHD::GRMHDSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
	PDE::NCP ncp(BgradQ);
	const Gradients g(gradQ);
	PDE(Q).nonConservativeProduct(g, ncp);
	ncp.zero_adm();
}

