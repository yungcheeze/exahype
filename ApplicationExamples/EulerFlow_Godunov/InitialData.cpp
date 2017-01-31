#include "InitialData.h"
#include "MyEulerSolver_Variables.h"


#include <cmath>

///// 2D /////

#ifdef Dim2
void EulerFVM::rarefactionWave(const double* const x,double* Q) {
  // Adiabatic constant 
  const double GAMMA = 1.4;
  
  Q[0] = 1.;
  // Velocities are set to zero.
  Q[1] = 0.;
  Q[2] = 0.;
  Q[3] = 0.;
  Q[4] =
      1. / (GAMMA -1) +
      std::exp(-((x[0] -0.5) *(x[0] -0.5) + (x[1] -0.5) *(x[1] -0.5) ) /
      (0.05 *0.05)) *
      1.0e-1;
}

void EulerFVM::sodShockTube(const double* const x,double* Q) {
  // Velocities are set to zero.
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  // Jump in density and pressure.
  if(x[0] < 0.5) {
    Q[0] = 1.0;
    Q[4] = 1.0;
  } else {
    Q[0] = 0.125;
    Q[4] = 0.1;
  }
}


void EulerFVM::explosionProblem(const double* const x,double* Q) {

  EulerFVM::MyEulerSolver::Variables var(Q);

  // @todo Das muss alles ins Guidebook

  // Velocities are set to zero (initially). 
/*
  var.v(0) = 0.001;
  var.v(1) = 0.001;
  var.v(2) = 0.001;
*/

/*
  tarch::la::Vector<3,double> tmp;
  tmp = 0.001, 0.001, 0.001;
  var.v( tmp );
*/

  var.v( 0.001, 0.001, 0.001 );

/*
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
*/

  // Circular shaped pressure jump at centre of domain.
  if((x[0] -0.5) *(x[0] -0.5) + (x[1] -0.5) *(x[1] -0.5) < 0.1) {
    var.rho() = 2.0;
    //Q[0] = 1.0;
    var.E() = 1.0;
    //Q[4] = 1.0;

    assertionNumericalEquals( Q[0], 2.0 );
    assertionNumericalEquals( Q[1], 0.001 );
    assertionNumericalEquals( Q[2], 0.001 );
    assertionNumericalEquals( Q[3], 0.001 );
    assertionNumericalEquals( Q[4], 1.0 );
  } else {
    var.rho() = 0.125;
    var.E()   = 0.1;
//    Q[4] = 0.1;

    assertionNumericalEquals( Q[0], 0.125 );
    assertionNumericalEquals( Q[1], 0.001 );
    assertionNumericalEquals( Q[2], 0.001 );
    assertionNumericalEquals( Q[3], 0.001 );
    assertionNumericalEquals( Q[4], 0.1 );
  }
}
#endif

///// 3D /////

#ifdef Dim3
void EulerFVM::rarefactionWave(const double* const x,double* Q) {
  // Adiabatic constant 
  const double GAMMA = 1.4;
  // Density is homogeneous.
  Q[0] = 1.;
  // Velocities are set to zero.
  Q[1] = 0.;
  Q[2] = 0.;
  Q[3] = 0.;
  Q[4] =
      1. / (GAMMA -1) +
      std::exp(-((x[0] -0.5) *(x[0] -0.5) + (x[1] -0.5) *(x[1] -0.5) + (x[2] -0.5) *(x[2] -0.5) ) /
      (0.05 *0.05)) *
      1.0e-1;
}

void EulerFVM::sodShockTube(const double* const x,double* Q) {
  // Velocities are set to zero.
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  // Jump in density and pressure.
  if(x[0] < 0.5) {
    Q[0] = 1.0;
    Q[4] = 1.0;
  } else {
    Q[0] = 0.125;
    Q[4] = 0.1;
  }
}


void EulerFVM::explosionProblem(const double* const x,double* Q) {

  Variables var(Q);


  // Density and velocities are set to zero (initially).
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  // Spherically shaped pressure jump at centre of domain.
  if((x[0] -0.5) *(x[0] -0.5) + (x[1] -0.5) *(x[1] -0.5) + (x[2] -0.5) *(x[2] -0.5) < 0.1) {
   var.rho = 1.0;   // war vorher Q[0] = 1.0;
   Q[4] = 1.0;
  } else {
    Q[0] = 0.125;
    Q[4] = 0.1;
  }
}
#endif

void EulerFVM::initialData(const double* const x,double* Q) {
  //rarefactionWave(x,Q);
  //sodShockTube(x,Q);
  explosionProblem(x,Q);
}
