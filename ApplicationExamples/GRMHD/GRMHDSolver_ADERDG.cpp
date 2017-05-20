#include "GRMHDSolver_ADERDG.h"

#include "GRMHDSolver_ADERDG_Variables.h"
#include "InitialData/InitialData.h"
#include "Fortran/PDE.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

#include "exahype/disableOptimization.h" // we experience compiler bugs sometimes.


const double excision_radius = 1.0;
constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
constexpr int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
constexpr int basisSize = order + 1;
constexpr int nDim = DIMENSIONS;

tarch::logging::Log GRMHD::GRMHDSolver_ADERDG::_log("GRMHDSolver_ADERDG");

void GRMHD::GRMHDSolver_ADERDG::init(std::vector<std::string>& cmdlineargs) {
  prepare_id();
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue  __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  bool insideExcisionBall = std::sqrt(center[0]*center[0] + center[1]*center[1] + center[2]*center[2]) < excision_radius;
  insideExcisionBall = false;
  bool hastoadjust = tarch::la::equals(t,0.0) || insideExcisionBall;
  return hastoadjust ? AdjustSolutionValue::PointWisely : AdjustSolutionValue::No;
}

void __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::adjustPointSolution(const double* const x,const double w,const double t,const double dt,double* Q) {
  id->Interpolate(x, t, Q);
  //printf("Interpoalted at x=[%f,%f,%f], t=%f, Q0=%f\n", x[0],x[1],x[2], t, Q[0]);
  for(int i=0; i<NumberOfVariables; i++) {
	if(!std::isfinite(Q[i])) {
		printf("NAN in i=%d at t=%f, x=[%f,%f,%f], Q[%d]=%f\n", i, t, x[0],x[1],x[2], i, Q[i]);
	}
  }
}

void __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void GRMHD::GRMHDSolver_ADERDG::flux(const double* const Q,double** F) {
  pdeflux_(F[0], F[1], (DIMENSIONS==3)?F[2]:nullptr, Q);
}


void __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::algebraicSource(const double* const Q,double* S) {
  pdesource_(S, Q);
}


void GRMHD::GRMHDSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int d,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
	double Qgp[nVar], F[nDim][nVar];

	std::memset(stateOut, 0, nVar * sizeof(double));
	std::memset(fluxOut, 0, nVar * sizeof(double));
  
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
			fluxOut[m] = -fluxIn[m];
		}
		
		break;
		
		case EXAHYPE_FACE_RIGHT:
		case EXAHYPE_FACE_BACK:
		case EXAHYPE_FACE_TOP:

		// Should probably use:
		//    stateOut = vacuum
		//    fluxOut = fluxIn
		// Use for the time being: Exact time integrated BC
		for(int i=0; i < basisSize; i++)  { // i == time
			const double weight = kernels::gaussLegendreWeights[order][i];
			const double xi = kernels::gaussLegendreNodes[order][i];
			double ti = t + xi * dt;

			id->Interpolate(x, ti, Qgp);
			pdeflux_(F[0], F[1], (nDim==3)?F[2]:nullptr, Qgp);
			
			for(int m=0; m < nVar; m++) {
				stateOut[m] += weight * Qgp[m];
				fluxOut[m] += weight * F[d][m];
			}
		}
		break;
		
		default:
			logError("boundaryValues", "faceIndex not supported");
			std::abort();
	}
}


exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

// only evaluated in Limiting context
bool GRMHD::GRMHDSolver_ADERDG::isPhysicallyAdmissible(const double* const QMin, const double* const QMax, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx, const double t, const double dt) const {

  // Static criterium for startup: When the density makes a large jump,
  // ie. at the star crust
  //if ( QMin[0] != 0.0 && QMax[0]/QMin[0] > 1e3 ) return false;

  if (QMin[0] < 0.0) return false;
  if (QMin[4] < 0.0) return false;

  for (int i=0; i<nVar; ++i) {
    if (!std::isfinite(QMin[i])) return false;
    if (!std::isfinite(QMax[i])) return false;
  }

  return true;
}


void __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
}


void GRMHD::GRMHDSolver_ADERDG::coefficientMatrix(const double* const Q,const int d,double* Bn) {
  // new scheme has no coefficient matrix
  static tarch::logging::Log _log("GRMHDSolver");
  logError("coefficientMatrix()", "ADERDG Coefficient Matrix invoked");
  exit(-2);

  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
}
