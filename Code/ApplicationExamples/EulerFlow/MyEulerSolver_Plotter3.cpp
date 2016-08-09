#include "MyEulerSolver_Plotter3.h"
#include "InitialData.h"


MyEulerSolver_Plotter3::MyEulerSolver_Plotter3() {
  // @todo Please insert your code here
}


MyEulerSolver_Plotter3::~MyEulerSolver_Plotter3() {
  // @todo Please insert your code here
}


void MyEulerSolver_Plotter3::startPlotting(double time) {
  this->time = time;
}


void MyEulerSolver_Plotter3::finishPlotting() {
  // @todo Please insert your code here
}


void MyEulerSolver_Plotter3::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  /**
   * This is the plotter for the exact solutions, given
   * as primitive Variables
   **/
  
  double xpos[DIMENSIONS];
  for(int i=0; i<DIMENSIONS; i++) xpos[i] = x[i];

  MovingGauss2D(xpos, outputQuantities, this->time);
  
  //cons2prim(outputQuantities, Q);

  /*
  for (int i=0; i<5; i++){ 
    outputQuantities[i] = Q[i];
  }
  */
}


