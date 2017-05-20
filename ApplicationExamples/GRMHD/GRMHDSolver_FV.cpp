#include "GRMHDSolver_FV.h"

#include "GRMHDSolver_FV_Variables.h"

#include "AbstractGRMHDSolver_ADERDG.h"
#include "InitialData/InitialData.h"
#include "Fortran/PDE.h"

#include <stdio.h>
#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

#include "exahype/disableOptimization.h" // bugs when limiting is on. whatevers

const double excision_radius = 1.0;

tarch::logging::Log GRMHD::GRMHDSolver_FV::_log("GRMHDSolver_FV");
constexpr int nVar = GRMHD::AbstractGRMHDSolver_FV::NumberOfVariables;


// enable nan tracker
#include <fenv.h>

void GRMHD::GRMHDSolver_FV::init(std::vector<std::string>& cmdlineargs) {
  // try NaN catcher
  feenableexcept(FE_INVALID | FE_OVERFLOW);  // Enable all floating point exceptions but FE_INEXACT
	
  prepare_id();
}

bool GRMHD::GRMHDSolver_FV::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  using namespace tarch::la;
  // Do excision only in 3D.
  bool insideExcisionBall = norm2(center) < excision_radius;
  insideExcisionBall = false;
  return tarch::la::equals(t,0.0) || insideExcisionBall;
}

void GRMHD::GRMHDSolver_FV::adjustSolution(const double* const x,const double w,const double t,const double dt, double* Q) {
  id->Interpolate(x, t, Q);
}

exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver_FV::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void GRMHD::GRMHDSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void GRMHD::GRMHDSolver_FV::flux(const double* const Q, double** F) {
  pdeflux_(F[0], F[1], (DIMENSIONS==3)?F[2]:nullptr, Q);
}


void GRMHD::GRMHDSolver_FV::algebraicSource(const double* const Q, double* S) {
  pdesource_(S, Q);
}

void GRMHD::GRMHDSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateIn,
    double* stateOut) {

	// Reflection BC at lower faces
	// face indices: 0 x=xmin 1 x=xmax, 2 y=ymin 3 y=ymax 4 z=zmin 5 z=zmax
	// corresponding 0-left, 1-right, 2-front, 3-back, 4-bottom, 5-top
	constexpr int EXAHYPE_FACE_LEFT = 0;
	constexpr int EXAHYPE_FACE_RIGHT = 1;
	constexpr int EXAHYPE_FACE_FRONT = 2;
	constexpr int EXAHYPE_FACE_BACK = 3;
	constexpr int EXAHYPE_FACE_BOTTOM = 4;
	constexpr int EXAHYPE_FACE_TOP = 5;
	
	switch(faceIndex) {
		case EXAHYPE_FACE_LEFT:
		case EXAHYPE_FACE_FRONT:
		case EXAHYPE_FACE_BOTTOM:

		// Reflection BC
		for(int m=0; m<nVar; m++) {
			stateOut[m] = stateIn[m];
			//fluxOut[m] = -fluxIn[m];
		}
		
		break;
		
		case EXAHYPE_FACE_RIGHT:
		case EXAHYPE_FACE_BACK:
		case EXAHYPE_FACE_TOP:

		// Should probably use:
		//    stateOut = vacuum
		//    fluxOut = fluxIn
		// Use for the time being: Exact BC
		id->Interpolate(x, t, stateOut);
		break;
		
		default:
			logError("boundaryValues", "faceIndex not supported");
			std::abort();
	}
	
	
  
}


void GRMHD::GRMHDSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
}

void GRMHD::GRMHDSolver_FV::coefficientMatrix(const double* const Q,const int d,double* Bn) {
  // new FV scheme has no coefficient matrix
  static tarch::logging::Log _log("GRMHDSolver");
  logError("coefficientMatrix()", "Coefficient Matrix invoked");
  exit(-2);

  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
}

