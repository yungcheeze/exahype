#include "PrimitivesWriter.h"


EulerADERDG::PrimitivesWriter::PrimitivesWriter(EulerSolver&  solver) {
  // @todo Please insert your code here
}


EulerADERDG::PrimitivesWriter::~PrimitivesWriter() {
  // @todo Please insert your code here
}


void EulerADERDG::PrimitivesWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void EulerADERDG::PrimitivesWriter::finishPlotting() {
  // @todo Please insert your code here
}


void EulerADERDG::PrimitivesWriter::mapQuantities(
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


