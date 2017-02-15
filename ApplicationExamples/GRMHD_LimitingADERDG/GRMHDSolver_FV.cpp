#include "GRMHDSolver_FV.h"

#include "GRMHDSolver_FV_Variables.h"

#include "AbstractGRMHDSolver_ADERDG.h"
#include "Fortran/InitialData.h"
#include "Fortran/PDE.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

void GRMHD::GRMHDSolver_FV::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

bool GRMHD::GRMHDSolver_FV::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) {
  // excsision is probably missing here
  return tarch::la::equals(t,0.0);
}

void GRMHD::GRMHDSolver_FV::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt, double* Q) {
  // Fortran
  initialdata_(x, &t, Q);
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
  pdeflux_(F[0], Q);
}


void GRMHD::GRMHDSolver_FV::source(const double* const Q, double* S) {
  pdesource_(S, Q);
}

void GRMHD::GRMHDSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateIn,
    double* stateOutside) {

  // Question: do we need time integration here?
  initialdata_(x, &t, stateOutside);
}

// Question: Are we missing NCP and BgradQ here?