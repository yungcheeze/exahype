/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include "MyEulerSolver.h"
#include "InitialData.h"

#include "MyEulerSolver_Variables.h"
#include "tarch/la/MatrixVectorOperations.h"

#include <memory>


void Euler::MyEulerSolver::init(std::vector<std::string>& cmdlineargs) {
}

void Euler::MyEulerSolver::flux(const double* const Q, double** F) {
  ReadOnlyVariables vars(Q);
  Fluxes f(F);

  tarch::la::Matrix<3,3,double> I;
  I = 1, 0, 0,
      0, 1, 0,
      0, 0, 1;

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  f.rho ( vars.j()                                 );
  f.j   ( irho * outerDot(vars.j(),vars.j()) + p*I );
  f.E   ( irho * (vars.E() + p) * vars.j()         );
}

void Euler::MyEulerSolver::eigenvalues(const double* const Q,
                                       const int normalNonZeroIndex,
                                       double* lambda) {
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


exahype::solvers::ADERDGSolver::AdjustSolutionValue Euler::MyEulerSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  // @todo Please implement/augment if required
  return tarch::la::equals(t,0.0) ? AdjustSolutionValue::PointWisely : AdjustSolutionValue::No;
}


void Euler::MyEulerSolver::adjustPointSolution(const double* const x,
                                                  const double w,const double t,const double dt, double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    DiffusingGauss(x,Q,t);
  }
}

exahype::solvers::Solver::RefinementControl
Euler::MyEulerSolver::refinementCriterion(
    const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level) {
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void Euler::MyEulerSolver::boundaryValues(const double* const x, const double t,const double dt,
                                          const int faceIndex,const int normalNonZero,
                                          const double* const fluxIn,const double* const stateIn,
                                          double* fluxOut, double* stateOut) {
  //  fluxOut
  //  //@todo Please implement
  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];
  fluxOut[3] = fluxIn[3];
  fluxOut[4] = fluxIn[4];
  //  // stateOut
  //  // @todo Please implement
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];
}
