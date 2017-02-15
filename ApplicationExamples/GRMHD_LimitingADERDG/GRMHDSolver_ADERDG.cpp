#include "GRMHDSolver_ADERDG.h"

#include "GRMHDSolver_ADERDG_Variables.h"
#include "Fortran/InitialData.h"
#include "Fortran/PDE.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

void GRMHD::GRMHDSolver_ADERDG::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

bool GRMHD::GRMHDSolver_ADERDG::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) {
  // excision is missing here, probably
  return tarch::la::equals(t,0.0);
}

void GRMHD::GRMHDSolver_ADERDG::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Fortran
  initialdata_(x, &t, Q);
}

void GRMHD::GRMHDSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void GRMHD::GRMHDSolver_ADERDG::flux(const double* const Q,double** F) {
  pdeflux_(F[0], Q);
}


void GRMHD::GRMHDSolver_ADERDG::source(const double* const Q,double* S) {
  pdesource_(S, Q);
}


void GRMHD::GRMHDSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int d,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {

  const int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
  const int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;

  double Qgp[nVar];
  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut, 0, nVar * sizeof(double));

  double F[3 * nVar]; // Fortran needs continous storage!
                      // Use always 3 dimensions here since the kernels works with those internally; see nDim in PDE.f90;
                      
  kernels::idx2 F_idx(nDim, nVar);

  for(int i=0; i < basisSize; i++)  { // i == time
     const double weight = kernels::gaussLegendreWeights[order][i];
     const double xi = kernels::gaussLegendreNodes[order][i];
     double ti = t + xi * dt;

     initialdata_(x, &ti, Qgp);
     pdeflux_(F, Qgp);
     for(int m=0; m < nVar; m++) {
        stateOut[m] += weight * Qgp[m];
        fluxOut[m] += weight * F[F_idx(d, m)];
     }
  }
}


exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


bool GRMHD::GRMHDSolver_ADERDG::physicalAdmissibilityDetection(const double* const QMin,const double* const QMax) {
  // Disable Limiter for AlfenWave:
  return false;
}


void GRMHD::GRMHDSolver_ADERDG::ncp(const double* const Q,const double* const gradQ,double* BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
}


void GRMHD::GRMHDSolver_ADERDG::matrixb(const double* const Q,const int d,double* Bn) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
}

bool GRMHD::GRMHDSolver_ADERDG::isDummyKRequired() const { 
  return false;
}


void GRMHD::GRMHDSolver_ADERDG::dummyK_Value(const double* const x,const double t,const double dt, double* forceVector, double* x0) {
  // ...
}