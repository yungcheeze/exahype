#include "MyEulerSolver.h"

#include "MyEulerSolver_Variables.h"

#include "LogoExaHyPE.h"


tarch::logging::Log EulerADERDG::MyEulerSolver::_log( "EulerADERDG::MyEulerSolver" );


void EulerADERDG::MyEulerSolver::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

void EulerADERDG::MyEulerSolver::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  if ( tarch::la::equals( t,0.0 ) ) {
    Variables vars(Q);

    tarch::la::Vector<DIMENSIONS,double> myX( x[0] - 0.06, 1.0-x[1] - 0.25 ); // translate
    myX *= static_cast<double>(LogoExaHyPE.width);
    tarch::la::Vector<DIMENSIONS,int>    myIntX( 1.2*myX(0) , 1.2*myX(1) );  // scale

    double Energy = 0.1;

    if (
      myIntX(0) > 0 && myIntX(0) < static_cast<int>(LogoExaHyPE.width)
      &&
      myIntX(1) > 0 && myIntX(1) < static_cast<int>(LogoExaHyPE.height)
    ) {
      Energy += 1.0-LogoExaHyPE.pixel_data[myIntX(1)*LogoExaHyPE.width+myIntX(0)];
    }

    vars.rho() = 1.0;
    vars.E()   = Energy;
    //vars.E()   = 0.1;
    vars.j(0,0,0);
  }
}

void EulerADERDG::MyEulerSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];

  double fi[DIMENSIONS][NumberOfVariables], *F[DIMENSIONS];
  for(int d=0; d<DIMENSIONS; d++) F[d] = fi[d];

  F[normalNonZero] = fluxOut; // This replaces the double pointer at pos normalNonZero by fluxOut.
  flux(stateOut, F);
}


exahype::solvers::Solver::RefinementControl EulerADERDG::MyEulerSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void EulerADERDG::MyEulerSolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  double u_n = Q[d + 1] * irho;
  double c  = std::sqrt(GAMMA * p * irho);

  eigs.rho()=u_n - c;
  eigs.E()  =u_n + c;
  eigs.j(u_n,u_n,u_n);
}


void EulerADERDG::MyEulerSolver::flux(const double* const Q,double** F) {
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




