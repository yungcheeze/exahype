#include "ComparisonVTKWriter.h"


Euler::ComparisonVTKWriter::ComparisonVTKWriter(MyEulerSolver&  solver) {
  // @todo Please insert your code here
}


Euler::ComparisonVTKWriter::~ComparisonVTKWriter() {
  // @todo Please insert your code here
}


void Euler::ComparisonVTKWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void Euler::ComparisonVTKWriter::finishPlotting() {
  // @todo Please insert your code here
}


void Euler::ComparisonVTKWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  // just for comparing...
  outputQuantities[0] = Q[4]; //energy
}


