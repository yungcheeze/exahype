/*
 * This is a stripped-down version of the GRMHD-AlfenWave for the super simple
 * benchmarking.
 * 
 * The AlfenWave does not care too much about grid size or so, so play around with
 * it. You also might turn on limiting anywhere in order to test what happens then.
 * 
 * At least ensure you have a reasonable minimum resolution, otherwise the sine wave
 * will explode.
 *
 */

#include "GRMHDSolver_ADERDG.h"

#include "GRMHDSolver_ADERDG_Variables.h"
#include "Fortran/PDE.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

#include "exahype/disableOptimization.h" // we experience compiler bugs sometimes.


constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
constexpr int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
constexpr int basisSize = order + 1;
constexpr int nDim = DIMENSIONS;
tarch::logging::Log GRMHD::GRMHDSolver_ADERDG::_log("GRMHDSolver_ADERDG");


void GRMHD::GRMHDSolver_ADERDG::init(std::vector<std::string>& cmdlineargs,exahype::Parser::ParserView constants) {
  // nothing to do here.
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue GRMHD::GRMHDSolver_ADERDG::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  return tarch::la::equals(t,0.0) ? AdjustSolutionValue::PointWisely : AdjustSolutionValue::No;
}

void GRMHD::GRMHDSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  alfenwave_(x, &t, Q);
}

void GRMHD::GRMHDSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void GRMHD::GRMHDSolver_ADERDG::flux(const double* const Q,double** F) {
  pdeflux_(F[0], F[1], (DIMENSIONS==3)?F[2]:nullptr, Q);
}


void GRMHD::GRMHDSolver_ADERDG::algebraicSource(const double* const Q,double* S) {
  pdesource_(S, Q);
}


void GRMHD::GRMHDSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int d,
  const double * const fluxIn,const double* const stateIn, double *fluxOut,double* stateOut) {

  // employ time-integrated exact BC for AlfenWave.
  
  double Qgp[nVar], Fs[nDim][nVar], *F[nDim];
  for(int d=0; d<nDim; d++) F[d] = Fs[d];
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
  return exahype::solvers::Solver::RefinementControl::Keep;
}

// only evaluated in Limiting context
void GRMHD::GRMHDSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(
    double* observables,const int numberOfObservables,
    const double* const Q) const {
  assertion(numberOfObservables==2);

  observables[0] = Q[0]; // rho
  observables[1] = Q[4]; // dens
}


bool GRMHD::GRMHDSolver_ADERDG::isPhysicallyAdmissible(
  const double* const solution,
  const double* const observablesMin,const double* const observablesMax,const int numberOfObservables,
  const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,
  const double t, const double dt) const {

  // geometric criterion:
  //  if ((center[0]-0.5)*(center[0]-0.5)+(center[1]-0.5)*(center[1]-0.5)<0.25*dx[0]*dx[0]) return false;

  // Static criterium for startup: When the density makes a large jump,
  // ie. at the star crust
  //if ( QMin[0] != 0.0 && QMax[0]/QMin[0] > 1e3 ) return false;

  if (observablesMin[0] <= 0.0) return false;
  if (observablesMin[1] < 0.0) return false;
  
  // what about this kind of check?
  /*

  for (int i=0; i<nVar; ++i) {
    if (!std::isfinite(QMin[i])) return false;
    if (!std::isfinite(QMax[i])) return false;
  }
  
  */
  return true;
}

void GRMHD::GRMHDSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
}


void GRMHD::GRMHDSolver_ADERDG::coefficientMatrix(const double* const Q,const int d,double* Bn) {
  // new scheme has no coefficient matrix
  static tarch::logging::Log _log("GRMHDSolver");
  logError("coefficientMatrix()", "ADERDG Coefficient Matrix invoked");
  exit(-2);

  // this works, thought.
  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
}
