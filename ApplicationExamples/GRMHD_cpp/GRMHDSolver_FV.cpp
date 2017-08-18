#include "GRMHDSolver_FV.h"
#include <algorithm> // fill_n

#include "peano/utils/Dimensions.h" // Defines DIMENSIONS
namespace GRMHD { constexpr int nVar = GRMHDSolver_FV::NumberOfVariables; }
#include "PDE/PDE.h"

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
}

exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver_FV::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
	return exahype::solvers::Solver::RefinementControl::Keep;
}


void GRMHD::GRMHDSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
	// Provide eigenvalues for the 9 SRMHD variables (D,S_j,tau,B^j),
	// we split off the 11 ADM material parameters (N^i,g_ij,detg)
	PDE::eigenvalues(Q, dIndex, lambda);
}

void GRMHD::GRMHDSolver_FV::flux(const double* const Q, double** F) {
	// Provide fluxes for the 9 SRMHD variables (D,S_j,tau,B^j),
	// we split off the 11 ADM material parameters (N^i,g_ij,detg)
	PDE(Q).flux(F);
}



void GRMHD::GRMHDSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {

	// TODO: SET BV for all 9 variables + 11 parameters.
}

void GRMHD::GRMHDSolver_FV::fusedSource(const double* const Q, const double* const gradQ, double* S) {
	PDE(Q).RightHandSide(gradQ, S);
}

