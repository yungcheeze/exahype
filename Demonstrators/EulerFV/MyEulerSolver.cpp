#include "MyEulerSolver.h"

#include "MyEulerSolver_Variables.h"


tarch::logging::Log EulerFV::MyEulerSolver::_log( "EulerFV::MyEulerSolver" );


void EulerFV::MyEulerSolver::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

void EulerFV::MyEulerSolver::adjustSolution(const double* const x,const double w,const double t,const double dt, double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 + #parameters
  
  // @todo Please implement/augment if required
  // State variables:
  Q[0] = 0.0;
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
}

exahype::solvers::Solver::RefinementControl EulerFV::MyEulerSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void EulerFV::MyEulerSolver::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 5 + #parameters
  
  // @todo Please implement/augment if required
  lambda[0] = 1.0;
  lambda[1] = 1.0;
  lambda[2] = 1.0;
  lambda[3] = 1.0;
  lambda[4] = 1.0;
}

void EulerFV::MyEulerSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 2
  // Number of variables    = 5 + #parameters
  
  // @todo Please implement/augment if required
  F[0][0] = 0.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;

  F[1][0] = 0.0;
  F[1][1] = 0.0;
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;
}



void EulerFV::MyEulerSolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
  // Dimensions             = 2
  // Number of variables    = 5 + #parameters

  // @todo Please implement/augment if required
  stateOutside[0] = 0.0;
  stateOutside[1] = 0.0;
  stateOutside[2] = 0.0;
  stateOutside[3] = 0.0;
  stateOutside[4] = 0.0;
}
