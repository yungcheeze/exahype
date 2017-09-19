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

void GRMHD::GRMHDSolver_FV::adjustSolution(const double* const x,const double w,const double t,const double dt, double* Q) {
	// TODO: CALL ID CODE HERE.
	// Set the 9 SRMHD variables (D,S_j,tau,B^j) and the 11 ADM material parameters (N^i,g_ij,detg)
	
	// currently, the C++ AlfenWave spills out primitive data
	double V[nVar];
	AlfenWave id(x,t,V);
	GRMHD::Prim2Cons(Q, V).copyFullStateVector();
	
	// also store the positions for debugging
	GRMHD::AbstractGRMHDSolver_FV::Variables var(Q);
	DFOR(i) var.pos(i) = x[i];
	var.check() = 15;
	
	//NVARS(m) printf("Qid[%d]=%e\n", m, Q[m]);

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

void GRMHD::GRMHDSolver_FV::flux(const double* const Q, double** F) {
	// Provide fluxes for the 9 SRMHD variables (D,S_j,tau,B^j),
	// we split off the 11 ADM material parameters (N^i,g_ij,detg)
	//PDE(Q).flux(F);
	
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
	NVARS(m) printf("FusedSource input: Q[%d]=%e\n", m, Q[m]);
	GRMHD::AbstractGRMHDSolver_FV::ReadOnlyVariables var(Q);
	std::cout << "pos: " << var.pos() << " check: " << var.check() <<  std::endl;
	
	PDE::Source S(S_);
	PDE(Q).RightHandSide(gradQ, S);
	S.zero_adm();
}

void GRMHD::GRMHDSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
	PDE::NCP ncp(BgradQ);
	const Gradients g(gradQ);
	PDE(Q).nonConservativeProduct(g, ncp);
	ncp.zero_adm();
}

