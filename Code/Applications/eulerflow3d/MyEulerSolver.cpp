#include "MyEulerSolver.h"

Euler3d::MyEulerSolver::MyEulerSolver( int kernelNumber):
  exahype::solvers::Solver("MyEulerSolver",exahype::solvers::Solver::ADER_DG,kernelNumber,5,3+1,exahype::solvers::Solver::GlobalTimeStepping) {
  // @todo Please implement/augment if required
}

int Euler3d::MyEulerSolver::getMinimumTreeDepth() const {
  // @todo Please implement
  return 3;
}

bool Euler3d::MyEulerSolver::hasToAdjustSolution(
    const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t) {
  // @todo Please implement/augment if required
  if (tarch::la::equals(t, 0.0, 1e-15)) {  // @todo precision
    return true;
  }
  return false;
}

void Euler3d::MyEulerSolver::adjustedSolutionValues(const double* const x,
                                                    const double w,
                                                    const double t,
                                                    const double dt,
                                                    double* Q) {
  // Dimensions             = 3
  // Number of variables    = 5
  if (tarch::la::equals(t, 0.0, 1e-15)) {  // @todo precision
    const double GAMMA = 1.4;
    // @todo Please implement
    Q[0] = 1.;
    Q[1] = 0.;
    Q[2] = 0.;
    Q[3] = 0.;
    Q[4] =
        1. / (GAMMA - 1) +
        std::exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5) +
                   (x[2] - 0.5) * (x[2] - 0.5)) /
                 (0.05)) *
            1.0e-3;
    //  Q[4] = 2.5;
  }
}



exahype::solvers::Solver::RefinementControl Euler3d::MyEulerSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::Keep;
}






//************************************************* 
//for FORTRAN kernels the fluxes and eigenvalues 
//have to be implemented in the file ./PDE.f90 
//************************************************* 



