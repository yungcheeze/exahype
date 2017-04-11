#include "h5writer.h"


Euler::h5writer::h5writer(MyEulerSolver&  solver) {
  // @todo Please insert your code here
}


Euler::h5writer::~h5writer() {
  // @todo Please insert your code here
}


void Euler::h5writer::startPlotting(double time) {
  // @todo Please insert your code here
}


void Euler::h5writer::finishPlotting() {
  // @todo Please insert your code here
}


void Euler::h5writer::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  // ONLY ONE OUTPUT QUANTITY
  outputQuantities[0] = Q[4]; // energy
}


