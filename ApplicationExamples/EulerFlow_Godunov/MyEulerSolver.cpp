#include "MyEulerSolver.h"

// @todo Has to be included by generator
#include "MyEulerSolver_Variables.h"


#include "InitialData.h"

// TODO(Dominic): Assess
//void EulerFVM::MyEulerSolver::init() {
//  // empty
//}

bool EulerFVM::MyEulerSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) {
  // @todo Please implement
  if ( tarch::la::equals(t,0.0) ) {
    // Tell kernel that you want to set initial conditions 
    return true;
  }
  else {
    // @todo Please implement
    return false; 
  }
}

void EulerFVM::MyEulerSolver::adjustedSolutionValues(const double* const x,
                                                  const double w,
                                                  const double t,
                                                  const double dt, double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    initialData(x,Q);
  } 
}

exahype::solvers::Solver::RefinementControl EulerFVM::MyEulerSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void EulerFVM::MyEulerSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // Dimensions             = 2/3
  // Number of variables    = 5 (#unknowns + #parameters)
  const double GAMMA = 1.4;

  double irho = 1.0 / Q[0];
  #ifdef Dim2
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);
  #else
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] *
  Q[3]) * irho);
  #endif

  double u_n = Q[normalNonZeroIndex + 1] * irho;
  double c = std::sqrt(GAMMA * p * irho);

  lambda[0] = u_n - c;
  lambda[1] = u_n;
  lambda[2] = u_n;
  lambda[3] = u_n;
  lambda[4] = u_n + c;
}

void EulerFVM::MyEulerSolver::flux(const double* const Q, double** F) {


  // Dimensions             = 2/3
  // Number of variables    = 5 (#unknowns + #parameters)
/*
  State state(Q);
  Flux flux(F);
*/


/*
  const double GAMMA = 1.4;
  const double irho = 1.0 / state.rho();

  double p = (GAMMA-1) * ( state.E() - 0.5 * irho*state.u()*state.u() );

  flux.writeRow(0, state.u());
  tarch::la::Vector<3,double> pContribution(0.0,0.0,0.0);
  for (int i=0; i<3; i++) {
    pContribution[i] = p;
    flux.writeRow(i+1, irho*state.u()[i]*state.u() + pContribution);
    pContribution[i] = 0.0;
  }
  flux.writeRow(1+3, irho*state.u()*(state.E()+p));
*/

  Fluxes flux(F);

  const double GAMMA = 1.4;
  const double irho = 1.0 / Q[0];
  #ifdef Dim2
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);
  #else
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * irho);
  #endif

//   @todo Please implement
//   f
  flux.rho( variables.u() );
  f[0] = Q[1];
  f[1] = irho * Q[1] * Q[1] + p;
  f[2] = irho * Q[1] * Q[2];
  f[3] = irho * Q[1] * Q[3];
  f[4] = irho * Q[1] * (Q[4] + p);
//   g
//   @todo Please implement
  g[0] = Q[2];
  g[1] = irho * Q[2] * Q[1];
  g[2] = irho * Q[2] * Q[2] + p;
  g[3] = irho * Q[2] * Q[3];
  g[4] = irho * Q[2] * (Q[4] + p);

#ifdef Dim3
  double* h = F[2];
//   h
//   @todo Please implement
  h[0] = Q[3];
  h[1] = irho * Q[3] * Q[1];
  h[2] = irho * Q[3] * Q[2];
  h[3] = irho * Q[3] * Q[3] + p;
  h[4] = irho * Q[3] * (Q[4] + p);
#endif


/*
  const double GAMMA = 1.4;
  const double irho = 1.0 / Q[0];
  #ifdef Dim2
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);
  #else
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * irho);
  #endif

  double* f = F[0];
  double* g = F[1];

//   @todo Please implement
//   f
  f[0] = Q[1];
  f[1] = irho * Q[1] * Q[1] + p;
  f[2] = irho * Q[1] * Q[2];
  f[3] = irho * Q[1] * Q[3];
  f[4] = irho * Q[1] * (Q[4] + p);
//   g
//   @todo Please implement
  g[0] = Q[2];
  g[1] = irho * Q[2] * Q[1];
  g[2] = irho * Q[2] * Q[2] + p;
  g[3] = irho * Q[2] * Q[3];
  g[4] = irho * Q[2] * (Q[4] + p);

#ifdef Dim3
  double* h = F[2];
//   h
//   @todo Please implement
  h[0] = Q[3];
  h[1] = irho * Q[3] * Q[1];
  h[2] = irho * Q[3] * Q[2];
  h[3] = irho * Q[3] * Q[3] + p;
  h[4] = irho * Q[3] * (Q[4] + p);
#endif
*/
}


void EulerFVM::MyEulerSolver::source(const double* const Q, double* S) {
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
}

void EulerFVM::MyEulerSolver::boundaryValues(
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



