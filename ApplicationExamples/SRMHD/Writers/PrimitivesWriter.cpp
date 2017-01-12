#include "Writers/PrimitivesWriter.h"
#include "C2P-MHD.h"

MHDSolver::PrimitivesWriter::PrimitivesWriter(MHDSolver&  solver) {}
MHDSolver::PrimitivesWriter::~PrimitivesWriter() {}
void MHDSolver::PrimitivesWriter::startPlotting(double time) {}
void MHDSolver::PrimitivesWriter::finishPlotting() {}

void MHDSolver::PrimitivesWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  // MHD primitive variable plotter
  int i;
  pdecons2prim_(outputQuantities, Q, &i);
}


