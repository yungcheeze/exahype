#include "H5Writer.h"


GRMHD::H5Writer::H5Writer(GRMHDSolver_FV&  solver) {
  // @todo Please insert your code here
}


GRMHD::H5Writer::~H5Writer() {
  // @todo Please insert your code here
}


void GRMHD::H5Writer::startPlotting(double time) {
  // @todo Please insert your code here
}


void GRMHD::H5Writer::finishPlotting() {
  // @todo Please insert your code here
}


void GRMHD::H5Writer::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<25; i++){ 
    outputQuantities[i] = (i<19 ? Q[i] : 0.0);
  }
}


