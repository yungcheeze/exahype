#include "MyEulerSolver.h"



Euler2d::MyEulerSolver::MyEulerSolver( int kernelNumber):
  exahype::solvers::Solver("MyEulerSolver",exahype::solvers::Solver::ADER_DG,kernelNumber,5,3+1,exahype::solvers::Solver::GlobalTimeStepping) {
  // @todo Please implement/augment if required
}

int Euler2d::MyEulerSolver::getMinimumTreeDepth() const {
  // @todo Please implement
  return 4;
}

bool Euler2d::MyEulerSolver::hasToAdjustSolution( const tarch::la::Vector<DIMENSIONS,double>&   center, const tarch::la::Vector<DIMENSIONS,double>&   dx, double t) {
  if (tarch::la::equals( t, 0.0 ,1e-15 )) { // @todo precision
    return true;
  }
  return false;
}

void Euler2d::MyEulerSolver::flux(const double* const Q, double* f, double* g) {
  // @todo Please implement
  const double GAMMA = 1.4;
  
  const double irho = 1.0/Q[0];
  const double p = (GAMMA-1)*( Q[4] - 0.5* (Q[1]*Q[1] + Q[2]*Q[2]) * irho );

  f[0] = Q[1];
  f[1] = irho*Q[1]*Q[1] + p;
  f[2] = irho*Q[1]*Q[2];
  f[3] = irho*Q[1]*Q[3];
  f[4] = irho*Q[1]*(Q[4]+p);

  g[0] = Q[2];
  g[1] = irho*Q[2]*Q[1];
  g[2] = irho*Q[2]*Q[2] + p;
  g[3] = irho*Q[2]*Q[3];
  g[4] = irho*Q[2]*(Q[4]+p);
}



void Euler2d::MyEulerSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // @todo Please implement
  const double GAMMA = 1.4;
  
  double irho = 1.0/Q[0];
  double p    = (GAMMA-1)*( Q[4] - 0.5* (Q[1]*Q[1] + Q[2]*Q[2]) * irho );

  double u_n = Q[normalNonZeroIndex+1] * irho;
  double c = std::sqrt(GAMMA * p * irho);

  lambda[0] = u_n-c;
  lambda[1] = u_n;
  lambda[2] = u_n;
  lambda[3] = u_n;
  lambda[4] = u_n+c;
}



void Euler2d::MyEulerSolver::adjustedSolutionValues(const double* const x, const double J_w,const double t, const double dt, double* Q) {
  if (tarch::la::equals( t, 0.0, 1e-15 )) { // @todo precision
    const double GAMMA = 1.4;

    Q[0] = 1.;
    Q[1] = 0.;
    Q[2] = 0.;
    Q[3] = 0.;
    Q[4] = 1./(GAMMA-1) + std::exp(-((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5))/(0.05*0.05)) * 1.0e-3;
  }
}



