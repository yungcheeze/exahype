#include "Writers/ConservedWriter.h"


MHDSolver::ConservedWriter::ConservedWriter(MHDSolver&  solver) {}
MHDSolver::ConservedWriter::~ConservedWriter() {}
void MHDSolver::ConservedWriter::startPlotting(double time) {}
void MHDSolver::ConservedWriter::finishPlotting() {}

void MHDSolver::ConservedWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<9; i++){ 
    // Conserved quantities, nothing happens here.
    outputQuantities[i] = Q[i];
  }
}


