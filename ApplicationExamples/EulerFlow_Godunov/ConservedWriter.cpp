#include "ConservedWriter.h"


EulerFVM::ConservedWriter::ConservedWriter(MyEulerSolver&  solver) {
  // @todo Please insert your code here
}


EulerFVM::ConservedWriter::~ConservedWriter() {
  // @todo Please insert your code here
}


void EulerFVM::ConservedWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void EulerFVM::ConservedWriter::finishPlotting() {
  // @todo Please insert your code here
}


void EulerFVM::ConservedWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<5; i++){ 
    outputQuantities[i] = Q[i];
  }
}


