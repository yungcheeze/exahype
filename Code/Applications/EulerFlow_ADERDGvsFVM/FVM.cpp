#include "FVM.h"

void Euler::FVM::flux(const double* const Q, double** F) {
  
  const double GAMMA = 1.4;

  const double irho = 1.0 / Q[0];
  const double p =
      (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);

  double* f = F[0];
  double* g = F[1];

  // @todo Please implement
  // f
  f[0] = Q[1];
  f[1] = irho * Q[1] * Q[1] + p;
  f[2] = irho * Q[1] * Q[2];
  f[3] = irho * Q[1] * Q[3];
  f[4] = irho * Q[1] * (Q[4] + p);
  // g
  // @todo Please implement
  g[0] = Q[2];
  g[1] = irho * Q[2] * Q[1];
  g[2] = irho * Q[2] * Q[2] + p;
  g[3] = irho * Q[2] * Q[3];
  g[4] = irho * Q[2] * (Q[4] + p);
}


void Euler::FVM::source(const double* const Q, double* S) {
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
}


void Euler::FVM::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  
  const double GAMMA = 1.4;

  double irho = 1.0 / Q[0];
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);

  double u_n = Q[normalNonZeroIndex + 1] * irho;
  double c = std::sqrt(GAMMA * p * irho);

  lambda[0] = u_n - c;
  lambda[1] = u_n;
  lambda[2] = u_n;
  lambda[3] = u_n;
  lambda[4] = u_n + c;
}


bool Euler::FVM::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) {
  if ( tarch::la::equals(t,0.0) ) {
    return true;
  }
  return false; 
}


exahype::solvers::Solver::RefinementControl Euler::FVM::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void Euler::FVM::adjustedSolutionValues(const double* const x, const double w, const double t, const double dt, double* Q) {  
/*
  Q[0] = 1.0;
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  if((x[0] -0.5) *(x[0] -0.5) + (x[1] -0.5) *(x[1] -0.5) < 0.1) {
    Q[4] = 1.0;
  } else {
    Q[4] = 0.1;
  }
  //Q[4] = floor(10/(1+10*sqrt((x[0] -0.5) *(x[0] -0.5) + (x[1] -0.5) *(x[1] -0.5))));
*/
    const double GAMMA = 1.4;
    Q[0] = 1.;
    Q[1] = 0.;
    Q[2] = 0.;
    Q[3] = 0.;
    Q[4] =
        1. / (GAMMA -1) +
        std::exp(-((x[0] -0.5) *(x[0] -0.5) + (x[1] -0.5) *(x[1] -0.5) ) /
        (0.05 *0.05)) *
        1.0e-1;

}

void Euler::FVM::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int normalNonZero,
    const double* const stateIn,
    double* stateOut) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)

    // Compute boundary state.
//    InitialData(x, stateOut, t);

    // This copy is not neccessary as we do have one component of
    // F already pointing to fluxOut.
    /*
    for (int i=0; i<5; i++) {
      fluxOut[i] = F[normalNonZero][i];
    }
    */

  // The problem with these definitions is that in a simulation
  // with a global nonzero velocity (as in MovingGauss2D), errnous
  // values move into the simulation domain very quickly. So these
  // boundary conditions are not good at all. Instead, we should
  // have per default "vacuum" boundary conditions, so that vacuum
  // values enter the grid as soon as matter moves away.

  //  // stateOut
  //  // @todo Please implement
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];
}
