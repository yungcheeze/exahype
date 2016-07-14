#include "MyEulerSolver.h"

Euler2d::MyEulerSolver::MyEulerSolver(int kernelNumber, int nodesPerCoordinateAxis, double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler)
  : exahype::solvers::Solver(
            "MyEulerSolver", exahype::solvers::Solver::Type::ADER_DG, kernelNumber, 5, 0, nodesPerCoordinateAxis, maximumMeshSize, timeStepping, std::move(profiler)) {
}

bool Euler2d::MyEulerSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  if (tarch::la::equals(t, 0.0)) {
    return true;
  }
  return false;
}

void Euler2d::MyEulerSolver::flux(const double *const Q, double **F) {
  double* f = F[0];
  double* g = F[1];

  const double GAMMA = 1.4;

  const double irho = 1.0 / Q[0];
  const double p =
      (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);

  f[0] = Q[1];
  f[1] = irho * Q[1] * Q[1] + p;
  f[2] = irho * Q[1] * Q[2];
  f[3] = irho * Q[1] * Q[3];
  f[4] = irho * Q[1] * (Q[4] + p);

  g[0] = Q[2];
  g[1] = irho * Q[2] * Q[1];
  g[2] = irho * Q[2] * Q[2] + p;
  g[3] = irho * Q[2] * Q[3];
  g[4] = irho * Q[2] * (Q[4] + p);
}

void Euler2d::MyEulerSolver::eigenvalues(const double *const Q,
                                         const int normalNonZeroIndex,
                                         double *lambda) {
  // @todo Please implement
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

void Euler2d::MyEulerSolver::adjustedSolutionValues(const double *const x,
                                                    const double w,
                                                    const double t,
                                                    const double dt,
                                                    double *Q) {
  if (tarch::la::equals(t, 0.0)) {
    const double GAMMA = 1.4;

    Q[0] = 1.;
    Q[1] = 0.;
    Q[2] = 0.;
    Q[3] = 0.;
    Q[4] =
        1. / (GAMMA - 1) +
        std::exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)) /
                 (0.05 * 0.05)) *
            1.0e-3;
  }
}

exahype::solvers::Solver::RefinementControl Euler2d::MyEulerSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
//  if (dx[0] > 1./81.) {
//    if (center[0] > 0 && center[0] < 0.3333) {
//      if (center[1] > 0 && center[1] < 0.3333) {
//        //        if (center[0] < 0.5 + std::pow(0.333333, deltaLevel) * .5 &&
//        //            center[0] > 0.5 - std::pow(0.333333, deltaLevel) * .5) {
//        //          if (center[1] < 0.5 + std::pow(0.333333, deltaLevel) * .5 &&
//        //              center[1] > 0.5 - std::pow(0.333333, deltaLevel) * .5) {
//        return exahype::solvers::Solver::RefinementControl::Refine;
//      }
//    }
//  }
//  
//  if (dx[0] > 1./81.) {
//    if (center[0] > 0.33333 && center[0] < 0.66667) {
//      if (center[1] > 0.33333 && center[1] < 0.66667) {
//        //        if (center[0] < 0.5 + std::pow(0.333333, deltaLevel) * .5 &&
//        //            center[0] > 0.5 - std::pow(0.333333, deltaLevel) * .5) {
//        //          if (center[1] < 0.5 + std::pow(0.333333, deltaLevel) * .5 &&
//        //              center[1] > 0.5 - std::pow(0.333333, deltaLevel) * .5) {
//        return exahype::solvers::Solver::RefinementControl::Refine;
//      }
//    }
//  }

//  if (dx[0] > 1./243.) {
//    if (center[0] > 0.66667 && center[0] < 1) {
//      if (center[1] > 0.33333 && center[1] < 1) {
//        //        if (center[0] < 0.5 + std::pow(0.333333, deltaLevel) * .5 &&
//        //            center[0] > 0.5 - std::pow(0.333333, deltaLevel) * .5) {
//        //          if (center[1] < 0.5 + std::pow(0.333333, deltaLevel) * .5 &&
//        //              center[1] > 0.5 - std::pow(0.333333, deltaLevel) * .5) {
//        return exahype::solvers::Solver::RefinementControl::Refine;
//      }
//    }
//  }
  return exahype::solvers::Solver::RefinementControl::Keep;
}
