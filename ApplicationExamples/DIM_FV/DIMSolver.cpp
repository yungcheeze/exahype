#include "DIMSolver.h"

#include "DIMSolver_Variables.h"

#include "InitialData.h"
#include "PDE.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"



void DIM::DIMSolver::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

bool DIM::DIMSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  return tarch::la::equals(t,0.0); 
}

void DIM::DIMSolver::adjustSolution(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Fortran
  initialdata_(x, &t, Q);
}

void DIM::DIMSolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void DIM::DIMSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int d,
  const double * const stateIn, double *stateOut) {

  initialdata_(x, &t, stateOut);
}  


exahype::solvers::Solver::RefinementControl DIM::DIMSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void  __attribute__((optimize("O0"))) DIM::DIMSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
  /*
  for(int i=0; i<NumberOfVariables; ++i) {
    if(!std::isfinite(BgradQ[i])) {
       printf("Nonfinite BgradQ[%d] = %f\n", i, BgradQ[i]);
    }
  }
  */
}


void DIM::DIMSolver::coefficientMatrix(const double* const Q,const int d,double* Bn) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
}

