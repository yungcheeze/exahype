#include "MyEulerSolver.h"

#include "MyEulerSolver_Variables.h"


tarch::logging::Log EulerADERDG::MyEulerSolver::_log( "EulerADERDG::MyEulerSolver" );


void EulerADERDG::MyEulerSolver::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

void EulerADERDG::MyEulerSolver::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  if ( tarch::la::equals( t,0.0 ) ) {
    Variables vars(Q);

    tarch::la::Vector<DIMENSIONS,double> myX( x[0] - 0.06, 1.0-x[1] - 0.25 ); // translate
    myX *= static_cast<double>(Image.width);
    tarch::la::Vector<DIMENSIONS,int>    myIntX( 1.2*myX(0) , 1.2*myX(1) );  // scale

    double Energy = 0.1;

    if (
      myIntX(0) > 0 && myIntX(0) < static_cast<int>(Image.width)
      &&
      myIntX(1) > 0 && myIntX(1) < static_cast<int>(Image.height)
    ) {
      Energy += 1.0-Image.pixel_data[myIntX(1)*Image.width+myIntX(0)];
    }

    vars.rho() = 1.0;
    vars.E()   = Energy;
    vars.j(0,0,0);
  }
}

void EulerADERDG::MyEulerSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 5 + 0

  // @todo Please implement/augment if required
  stateOut[0] = 0.0;
  stateOut[1] = 0.0;
  stateOut[2] = 0.0;
  stateOut[3] = 0.0;
  stateOut[4] = 0.0;

  fluxOut[0] = 0.0;
  fluxOut[1] = 0.0;
  fluxOut[2] = 0.0;
  fluxOut[3] = 0.0;
  fluxOut[4] = 0.0;
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
  // Dimensions                        = 2
  // Number of variables + parameters  = 5 + 0
  
  // @todo Please implement/augment if required
  lambda[0] = 1.0;
  lambda[1] = 1.0;
  lambda[2] = 1.0;
  lambda[3] = 1.0;
  lambda[4] = 1.0;
}


void EulerADERDG::MyEulerSolver::flux(const double* const Q,double** F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 5 + 0
  
  // @todo Please implement/augment if required
  F[0][0] = 0.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;
  
  F[1][0] = 0.0;
  F[1][1] = 0.0;
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;
  
}




