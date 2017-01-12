#include "ErrorWriter.h"


MHDSolver::ErrorWriter::ErrorWriter(MHDSolver&  solver) {
  // @todo Please insert your code here
}


MHDSolver::ErrorWriter::~ErrorWriter() {
  // @todo Please insert your code here
}


void MHDSolver::ErrorWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void MHDSolver::ErrorWriter::finishPlotting() {
  // @todo Please insert your code here
}


void MHDSolver::ErrorWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<9; i++){ 
    outputQuantities[i] = Q[i];
  }
}


