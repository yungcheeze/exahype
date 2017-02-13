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
  // @todo Please implement/augment if required
  return tarch::la::equals(t,0.0);
}

void GRMHD::GRMHDSolver_FV::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt, double* Q) {
  // Fortran
  adjustedsolutionvalues_(x, &w, &t, &dt, Q);
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

  // TODO: Exact BC
  stateOutside[ 0] = 0.0;
  stateOutside[ 1] = 0.0;
  stateOutside[ 2] = 0.0;
  stateOutside[ 3] = 0.0;
  stateOutside[ 4] = 0.0;
  stateOutside[ 5] = 0.0;
  stateOutside[ 6] = 0.0;
  stateOutside[ 7] = 0.0;
  stateOutside[ 8] = 0.0;
  stateOutside[ 9] = 0.0;
  stateOutside[10] = 0.0;
  stateOutside[11] = 0.0;
  stateOutside[12] = 0.0;
  stateOutside[13] = 0.0;
  stateOutside[14] = 0.0;
  stateOutside[15] = 0.0;
  stateOutside[16] = 0.0;
  stateOutside[17] = 0.0;
  stateOutside[18] = 0.0;
}