#include "ExactPrimitivesWriter.h"

#include "EulerSolver.h"

Euler::ExactPrimitivesWriter::ExactPrimitivesWriter(EulerSolver&  solver) {
  // @todo Please insert your code here
}


Euler::ExactPrimitivesWriter::~ExactPrimitivesWriter() {
  // @todo Please insert your code here
}


void Euler::ExactPrimitivesWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void Euler::ExactPrimitivesWriter::finishPlotting() {
  // @todo Please insert your code here
}


void Euler::ExactPrimitivesWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  EulerSolver::entropyWave(x.data(),timeStamp,outputQuantities);
}


