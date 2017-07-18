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

#include "EulerSolver_ADERDG.h"
#include "EulerSolver_ADERDG_Variables.h"
#include "tarch/la/MatrixVectorOperations.h"

#include <algorithm>

#include <math.h>

#include "kernels/GaussLegendreQuadrature.h"

void Euler::EulerSolver_ADERDG::init(std::vector<std::string>& cmdlineargs) {
}

void Euler::EulerSolver_ADERDG::flux(const double* const Q, double** F) {
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

void Euler::EulerSolver_ADERDG::eigenvalues(const double* const Q,
                                       const int direction,
                                       double* lambda) {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  double u_n = Q[direction + 1] * irho;
  double c  = std::sqrt(GAMMA * p * irho);

  eigs.rho()=u_n - c;
  eigs.E()  =u_n + c;
  eigs.j(u_n,u_n,u_n);
}


exahype::solvers::ADERDGSolver::AdjustSolutionValue Euler::EulerSolver_ADERDG::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  return tarch::la::equals(t,0.0) ? AdjustSolutionValue::PointWisely : AdjustSolutionValue::No;
}

void Euler::EulerSolver_ADERDG::entropyWave(const double* const x,double t, double* Q) {
  const double GAMMA     = 1.4;
  constexpr double width = 0.3;

  #if DIMENSIONS==2
  tarch::la::Vector<DIMENSIONS,double> xVec(x[0],x[1]);
  tarch::la::Vector<DIMENSIONS,double> v0(0.5,0.0);
  tarch::la::Vector<DIMENSIONS,double> x0(0.5,0.5);
  #else
  tarch::la::Vector<DIMENSIONS,double> xVec(x[0],x[1],x[2]);
  tarch::la::Vector<DIMENSIONS,double> v0(0.5,0.0,0.0);
  tarch::la::Vector<DIMENSIONS,double> x0(0.5,0.5,0.5);
  #endif
  const double distance  = tarch::la::norm2( xVec - x0 - v0 * t );

  Q[0] = 0.5 + 0.3 * std::exp(-distance / std::pow(width, DIMENSIONS));
  Q[1] = Q[0] * v0[0];
  Q[2] = Q[0] * v0[1];
  Q[3] = 0.0;
  // total energy = internal energy + kinetic energy
  const double p = 1.;
  Q[4] = p / (GAMMA-1)   +  0.5*Q[0] * (v0[0]*v0[0]+v0[1]*v0[1]); // v*v; assumes: v0[2]=0
}

void Euler::EulerSolver_ADERDG::adjustPointSolution(const double* const x,
                                                  const double w,const double t,const double dt, double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    entropyWave(x,0.0,Q);
  }
}

exahype::solvers::Solver::RefinementControl
Euler::EulerSolver_ADERDG::refinementCriterion(
    const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level) {
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void Euler::EulerSolver_ADERDG::boundaryValues(const double* const x, const double t,const double dt,
                                          const int faceIndex,const int direction,
                                          const double* const fluxIn,const double* const stateIn,
                                          double* fluxOut, double* stateOut) {
  double Q[NumberOfVariables]     = {0.0};
  double _F[3][NumberOfVariables] = {0.0};
  double* F[3] = {_F[0],_F[1],_F[2]};

  // initialise
  std::fill_n(stateOut, NumberOfVariables, 0.0);
  std::fill_n(fluxOut,  NumberOfVariables, 0.0);
  for (int i=0; i<Order+1; i++) {
    const double ti = t + dt * kernels::gaussLegendreNodes[Order][i];

    entropyWave(x,ti,Q);
    flux(Q,F);
    for (int v=0; v<NumberOfVariables; v++) {
      stateOut[v] += Q[v]            * kernels::gaussLegendreWeights[Order][i];
      fluxOut[v]  += F[direction][v] * kernels::gaussLegendreWeights[Order][i];
    }
  }
}
