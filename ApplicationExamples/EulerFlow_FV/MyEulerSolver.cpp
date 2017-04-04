#include "MyEulerSolver.h"

// @todo Has to be included by generator
#include "MyEulerSolver_Variables.h"
#include "kernels/KernelUtils.h"

#include "InitialData.h"
#include "tarch/la/Vector.h"

void Euler::MyEulerSolver::init(std::vector<std::string>& cmdlineargs) {
}

bool Euler::MyEulerSolver::useAdjustSolution(
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& dx,
      const double t,
      const double dt) const  {
  return tarch::la::equals(t,0.0); // ? exahype::solvers::Solver::AdjustSolutionValue::Pointwisely : exahype::solvers::Solver::AdjustSolutionValue::No;
}

void Euler::MyEulerSolver::adjustSolution(const double* const x,
                                               const double w,
                                               const double t,
                                               const double dt, double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    explosionProblem(x,Q);
  } 
}

exahype::solvers::Solver::RefinementControl Euler::MyEulerSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void Euler::MyEulerSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
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

  //f.rho ( vars.j()                                 );
  f.rho( tarch::la::Vector<DIMENSIONS,double>(0.0) ); // handle rho with NCP
  f.j   ( irho * outerDot(vars.j(),vars.j()) + p*I );
  f.E   ( irho * (vars.E() + p) * vars.j()         );
}


void Euler::MyEulerSolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int normalNonZero,
    const double* const stateIn,
    double* stateOut) {

    // Compute boundary state.
    explosionProblem(x, stateOut);
}

void Euler::MyEulerSolver::algebraicSource(const double* const Q, double* S) {
  for(int l=0; l<NumberOfVariables; l++) {
    S[l] = 0.0;
  }
}

void Euler::MyEulerSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  for(int l=0; l<NumberOfVariables; l++) {
    BgradQ[l] = 0.0;
  }
}

void Euler::MyEulerSolver::coefficientMatrix(const double* const Q,const int d,double* Bn) {
  for(int l=0; l<NumberOfVariables*NumberOfVariables; l++) {
    Bn[l] = 0.0;
  }
}

