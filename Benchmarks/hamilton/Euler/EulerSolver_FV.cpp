/**
 * This file is part of the ExaHyPE project. For copyright and information
 * please see www.exahype.eu.
 */
#include "EulerSolver_FV.h"
#include "EulerSolver_FV_Variables.h"

#include "EulerSolver_ADERDG.h"

void Euler::EulerSolver_FV::init(std::vector<std::string>& cmdlineargs, exahype::Parser::ParserView& constants) {
}


bool Euler::EulerSolver_FV::useAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) const {
  //  return true; // comment in for vel_y cleaning; do the same in the ADER-DG solver
  return tarch::la::equals(t,0.0);
}

void Euler::EulerSolver_FV::adjustSolution(
    const double* const x,
    const double w,
    const double t,
    const double dt, double* Q) {
  if (tarch::la::equals(t,0.0)) {
    EulerSolver_ADERDG::referenceSolution(x,t,Q);
  }
  // else {  // comment in for vel_y cleaning; do the same in the ADER-DG solver
  //   Q[2] = 0.0;
  // }
}


exahype::solvers::Solver::RefinementControl Euler::EulerSolver_FV::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void Euler::EulerSolver_FV::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  double u_n = Q[normalNonZeroIndex + 1] * irho;
  double c  = std::sqrt(GAMMA * p * irho);

  eigs.rho()=u_n - c;
  eigs.E()  =u_n + c;
  eigs.j(u_n,u_n,u_n);
}


void Euler::EulerSolver_FV::flux(const double* const Q, double** F) {
  ReadOnlyVariables vars(Q);
  Fluxes fluxes(F);

  tarch::la::Matrix<3,3,double> I;
  I = 1, 0, 0,
      0, 1, 0,
      0, 0, 1;

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  fluxes.rho ( vars.j()                                 );
  fluxes.j   ( irho * outerDot(vars.j(),vars.j()) + p*I );
  fluxes.E   ( irho * (vars.E() + p) * vars.j()         );
}


void Euler::EulerSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int direction,
    const double* const stateInside,
    double* stateOutside) {
  switch (EulerSolver_ADERDG::ReferenceChoice) {
  case EulerSolver_ADERDG::Reference::SodShockTube:
    std::copy_n(stateInside, NumberOfVariables, stateOutside);
//    stateOutside[1+direction] =  -stateOutside[1+direction];
    break;
  case EulerSolver_ADERDG::Reference::EntropyWave:
    EulerSolver_ADERDG::referenceSolution(x,t,stateOutside);
    break;
  }
}
