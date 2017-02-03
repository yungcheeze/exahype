#include "MyEulerSolver.h"
#include "MyEulerSolver_Variables.h"
#include "Logo.h"


bool EulerFV::MyEulerSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) {
  return tarch::la::equals(t,0.0);
}


void EulerFV::MyEulerSolver::adjustedSolutionValues(const double* const x,
                                                  const double w,
                                                  const double t,
                                                  const double dt, double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
  // State variables:
  tarch::la::Vector<DIMENSIONS,double> myX( x[0], 1.0-x[1] );
  myX *= static_cast<double>(Image.width);

  tarch::la::Vector<DIMENSIONS,int>    myIntX( myX(0), myX(1) );

  double rho = 0.1;

  if (
    myIntX(0) < static_cast<int>(Image.width)
    &&
    myIntX(1) < static_cast<int>(Image.height)
  ) {
    rho += (
        Image.pixel_data[myIntX(1)*Image.width*3+myIntX(0)*3+0]
      + Image.pixel_data[myIntX(1)*Image.width*3+myIntX(0)*3+1]
      + Image.pixel_data[myIntX(1)*Image.width*3+myIntX(0)*3+2]) / 3.0 / 256.0;
  }

  //if (tarch::la::equals(t, 0.0)) {
    Q[0] = rho;
    Q[1] = 0.0;
    Q[2] = 0.0;
    Q[3] = 0.0;
    Q[4] = 1.0;
  //}
}

exahype::solvers::Solver::RefinementControl EulerFV::MyEulerSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void EulerFV::MyEulerSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
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


void EulerFV::MyEulerSolver::flux(const double* const Q, double** F) {
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


void EulerFV::MyEulerSolver::source(const double* const Q, double* S) {
  Variables source(S);
  source.rho()=0;
  source.E()=0;
  source.j(0,0,0);
}


void EulerFV::MyEulerSolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int normalNonZero,
    const double* const stateInside,
    double* stateOutside) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
}
