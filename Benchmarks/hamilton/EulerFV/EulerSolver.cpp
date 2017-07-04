/**
 * This file is part of the ExaHyPE project. For copyright and information
 * please see www.exahype.eu.
 */
#include "EulerSolver.h"
#include "EulerSolver_Variables.h"

void EulerFV::EulerSolver::init(std::vector<std::string>& cmdlineargs) {
}


bool EulerFV::EulerSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) const {
  return tarch::la::equals(t,0.0);
}

/**
 * to solve Sod's Shock Tube problem
 * reference: "http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html"
 *   |             |         |       |         |
 *   |             |         |       |         |
 *   | rarefaction | contact | shock |         |
 *___|_____________|_________|_______|_________|_______________
 *   x1           x2         x0      x3       x4
 */
void EulerFV::EulerSolver::sodShockTube(const double* const x, const double t, double* Q) {
  // Initial conditions
  constexpr double gamma     =1.39999999999999991118;
  constexpr double mu        =0.40824829046386296172;
  constexpr double mu2       =mu*mu;
  constexpr double x0        =0.50000000000000000000;

  constexpr double rho_r     =0.12500000000000000000;
  constexpr double P_r       =0.10000000000000000555;
  constexpr double u_r       =0.00000000000000000000;
  constexpr double rho_l     =1.00000000000000000000;
  constexpr double P_l       =1.00000000000000000000;
  constexpr double u_l       =0.00000000000000000000;

  // Speed of sound
  constexpr double c_l       =1.18321595661992318149;
  // constexpr double c_r       =1.05830052442583610883;

  constexpr double P_post    =0.31900053530972960480;
  constexpr double v_post    =0.89095233349579605608;
  constexpr double rho_post  =0.27393934150152593476;
  constexpr double v_shock   =1.63869997736350780926;
  constexpr double rho_middle=0.44214558642791562670;

  // Key Positions
  const double x1 = x0 - c_l*t;
  const double x3 = x0 + v_post*t;
  const double x4 = x0 + v_shock*t;
  // Determining x2
  constexpr  double c_2 = c_l - ((gamma - 1)/2)*v_post;
  const      double x2 = x0 + (v_post - c_2)*t;

  double p = 0; // pressure
  Q[2] = 0; // y velocity
  Q[3] = 0; // z velocity
  if (t>0) {
    if (x[0] < x1) {
      Q[0] = rho_l;
      Q[1] = Q[0] * u_l;
      p    = P_l;
    } else if (x1 <= x[0] && x[0] < x2) {
      // rarefaction wave
      const double c = mu2*((x0 - x[0])/t) + (1 - mu2)*c_l;
      Q[0] = rho_l*std::pow(c/c_l,2/(gamma - 1));
      Q[1] = Q[0] * (1 - mu2)*( (-(x0-x[0])/t) + c_l);
      p    = P_l*power((Q[0]/rho_l),gamma);
    } else if (x2 <= x[0] && x[0] < x3) {
      Q[0] = rho_middle;
      Q[1] = Q[0] * v_post;
      p    = P_post;
    } else if (x3 <= x[0] && x[0] < x4) {
      Q[0] = rho_post;
      Q[1] = Q[0] * v_post;
      p    = P_post;
    } else if (x4 <= x[0]) {
      Q[0] = rho_r;
      Q[1] = Q[0] * u_r;
      p    = P_r;
    }
    // total energy = internal energy + kinetic energy
    Q[4] = p/(gamma-1) + 0.5*Q[0] * (Q[1]*Q[1]); // v*v; assumes: Q[1+i]=0, i=1,2.
  } else {
    if (x[0] < x1) {
      Q[0] = rho_l;
      Q[1] = Q[0] * u_l;
      p    = P_l;
    } else {
      Q[0] = rho_r;
      Q[1] = Q[0] * u_r;
      p    = P_r;
    }

    Q[4] = p/(gamma-1) + 0.5 / Q[0] * (Q[1]*Q[1]); // j*j, j=rho*v !!! ; assumes: Q[1+i]=0, i=1,2.
  }

//  std::cout << "x=" << x[0] << ": ";
//  std::cout << "Q=";
//  std::cout << Q[0] << ",";
//  std::cout << Q[1] << ",";
//  std::cout << Q[2] << ",";
//  std::cout << Q[3] << ",";
//  std::cout << Q[4] << std::endl;
}

void EulerFV::EulerSolver::adjustSolution(
    const double* const x,
    const double w,
    const double t,
    const double dt, double* Q) {
  sodShockTube(x,t,Q);
}


exahype::solvers::Solver::RefinementControl EulerFV::EulerSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void EulerFV::EulerSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
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


void EulerFV::EulerSolver::flux(const double* const Q, double** F) {
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


void EulerFV::EulerSolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int normalNonZero,
    const double* const stateInside,
    double* stateOutside) {
  ReadOnlyVariables varsInside(stateInside);
  Variables         varsOutside(stateOutside);

  varsOutside = varsInside;
}
