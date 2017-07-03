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

#include "EulerSolver.h"
#include "EulerSolver_Variables.h"
#include "tarch/la/MatrixVectorOperations.h"

#include <algorithm>

#include <math.h>

#include "kernels/GaussLegendreQuadrature.h"

void EulerADERDG::EulerSolver::init(std::vector<std::string>& cmdlineargs) {
}

void EulerADERDG::EulerSolver::flux(const double* const Q, double** F) {
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

void EulerADERDG::EulerSolver::eigenvalues(const double* const Q,
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

//  std::cout << "p="<<c<<std::endl;
//  std::cout << "c="<<c<<std::endl;
//
//  for (int i=0; i<NumberOfVariables; i++) {
//    std::cout << "eigs["<<i<<"]="<<eigs.data()[i]<<",";
//  }
//  std::cout << std::endl;
}


exahype::solvers::ADERDGSolver::AdjustSolutionValue EulerADERDG::EulerSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  // @todo Please implement/augment if required
  return tarch::la::equals(t,0.0) ? AdjustSolutionValue::PointWisely : AdjustSolutionValue::No;
}

void EulerADERDG::EulerSolver::entropyWave(const double* const x,double t, double* Q) {
  constexpr double A      = 1.0;
  constexpr double rhoInf = 3.0; // rhoInf-A > 0
  constexpr double uInf   = 1.0;
  constexpr double vInf   = 1.0;
  constexpr double wInf   = 1.0;
  constexpr double EInf   = 10.0;  // p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() ) > 0
  constexpr double freq   = 4.0;

  Q[4] = EInf;
#if DIMENSIONS==2
  constexpr double qInf = uInf+vInf;
  Q[0] = rhoInf + A * std::sin( freq * M_PI*(x[0]+x[1] - qInf*t) );
  Q[1] = Q[0] * uInf;
  Q[2] = Q[0] * vInf;
  Q[3] = Q[0] * wInf;
#else
  constexpr double qInf   = uInf+vInf+wInf;
  Q[0] = rhoInf + A * std::sin( freq * M_PI*(x[0]+x[1]+x[2] - qInf*t) );
  Q[1] = Q[0] * uInf;
  Q[2] = Q[0] * vInf;
  Q[3] = Q[0] * wInf;
#endif
}

void EulerADERDG::EulerSolver::adjustPointSolution(const double* const x,
                                                  const double w,const double t,const double dt, double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    entropyWave(x,0.0,Q);
  }
}

exahype::solvers::Solver::RefinementControl
EulerADERDG::EulerSolver::refinementCriterion(
    const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level) {
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void EulerADERDG::EulerSolver::boundaryValues(const double* const x, const double t,const double dt,
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
//  std::copy_n(stateIn, NumberOfVariables, stateOut);
//  std::copy_n(fluxIn,  NumberOfVariables, fluxOut);
}
