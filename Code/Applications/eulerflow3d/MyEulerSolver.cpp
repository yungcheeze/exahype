#include "MyEulerSolver.h"



Euler3d::MyEulerSolver::MyEulerSolver( int kernelNumber):
  exahype::solvers::Solver("MyEulerSolver",exahype::solvers::Solver::ADER_DG,kernelNumber,5,3+1) {
  // @todo Please implement/augment if required
}



int Euler3d::MyEulerSolver::getMinimumTreeDepth() const {
  // @todo Please implement
  return 3;
}
void Euler3d::MyEulerSolver::flux(const double* const Q, double* f, double* g, double* h) {

  // Dimensions             = 3
  // Number of variables    = 5

  const double GAMMA = 1.4;
  
  const double irho = 1.0/Q[0];
  const double p = (GAMMA-1)*( Q[4] - 0.5* (Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) * irho );

  // f
  // @todo Please implement
  f[0] = Q[1];
  f[1] = irho*Q[1]*Q[1] + p;
  f[2] = irho*Q[1]*Q[2];
  f[3] = irho*Q[1]*Q[3];
  f[4] = irho*Q[1]*(Q[4]+p);
  // g
  // @todo Please implement
  g[0] = Q[2];
  g[1] = irho*Q[2]*Q[1];
  g[2] = irho*Q[2]*Q[2] + p;
  g[3] = irho*Q[2]*Q[3];
  g[4] = irho*Q[2]*(Q[4]+p);
  // h
  // @todo Please implement
  h[0] = Q[3];
  h[1] = irho*Q[3]*Q[1];
  h[2] = irho*Q[3]*Q[2];
  h[3] = irho*Q[3]*Q[3] + p;
  h[4] = irho*Q[3]*(Q[4]+p);
}



void Euler3d::MyEulerSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // Dimensions             = 3
  // Number of variables    = 5
  const double GAMMA = 1.4;
  
  double irho = 1.0/Q[0];
  double p    = (GAMMA-1)*( Q[4] - 0.5* (Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) * irho );

  double u_n = Q[normalNonZeroIndex+1] * irho;
  double c = std::sqrt(GAMMA * p * irho);
  // @todo Please implement
  lambda[0] = u_n-c;
  lambda[1] = u_n;
  lambda[2] = u_n;
  lambda[3] = u_n;
  lambda[4] = u_n+c;
}



void Euler3d::MyEulerSolver::initialValues(const double* const x, double* Q) {
  // Dimensions             = 3
  // Number of variables    = 5
  const double GAMMA = 1.4;
  // @todo Please implement
  Q[0] = 1.;
  Q[1] = 0.;
  Q[2] = 0.;
  Q[3] = 0.;
//  Q[4] = 1./(GAMMA-1) + std::exp(-((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/(0.05*0.05)) * 1.0e-3;
  Q[4] = (x[0] + x[1] + x[2]);
}



