#include "GRMHDSolver.h"

#include "GRMHDSolver_Variables.h"

#include "Fortran/InitialData.h"
#include "Fortran/PDE.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"



void GRMHD::GRMHDSolver::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue GRMHD::GRMHDSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  const double excision_radius = 1.0;
  bool insideExcisionBall = std::sqrt(center[0]*center[0] + center[1]*center[1] + center[2]*center[2]) < excision_radius;
  //bool insideExcisionBall = false;
  bool hastoadjust = tarch::la::equals(t,0.0) || insideExcisionBall;
  return hastoadjust ? AdjustSolutionValue::PointWisely : AdjustSolutionValue::No;
}

void GRMHD::GRMHDSolver::adjustPointSolution(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Fortran
  initialdata_(x, &t, Q);
}

void GRMHD::GRMHDSolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void GRMHD::GRMHDSolver::flux(const double* const Q,double** F) {
  pdeflux_(F[0], F[1], (DIMENSIONS==3) ? F[2] : nullptr, Q);
}


void GRMHD::GRMHDSolver::algebraicSource(const double* const Q,double* S) {
  pdesource_(S, Q);
}


void GRMHD::GRMHDSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int d,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {

  const int nVar = GRMHD::AbstractGRMHDSolver::NumberOfVariables;
  const int order = GRMHD::AbstractGRMHDSolver::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;

  double Qgp[nVar], F[nDim][nVar];

  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut, 0, nVar * sizeof(double));
  
  for(int i=0; i < basisSize; i++)  { // i == time
     const double weight = kernels::gaussLegendreWeights[order][i];
     const double xi = kernels::gaussLegendreNodes[order][i];
     double ti = t + xi * dt;

     initialdata_(x, &ti, Qgp);
     pdeflux_(F[0], F[1], (nDim==3)?F[2]:nullptr, Qgp);
     for(int m=0; m < nVar; m++) {
        stateOut[m] += weight * Qgp[m];
        fluxOut[m] += weight * F[d][m];
     }
  }
}


exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void GRMHD::GRMHDSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
  /*
  for(int i=0; i<NumberOfVariables; ++i) {
    if(!std::isfinite(BgradQ[i])) {
       printf("Nonfinite BgradQ[%d] = %f\n", i, BgradQ[i]);
    }
  }
  */
}


void GRMHD::GRMHDSolver::coefficientMatrix(const double* const Q,const int d,double* Bn) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
}

