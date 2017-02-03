#include "Writers/ExactPrimitivesWriter.h"
#include "InitialDataAdapter.h"
#include "C2P-MHD.h"

MHDSolver::ExactPrimitivesWriter::ExactPrimitivesWriter(MHDSolver&  solver) {}
MHDSolver::ExactPrimitivesWriter::~ExactPrimitivesWriter() {}
void MHDSolver::ExactPrimitivesWriter::startPlotting(double time) {}
void MHDSolver::ExactPrimitivesWriter::finishPlotting() {}

void MHDSolver::ExactPrimitivesWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  const double *xpos = x.data();
  double Conserved[MHDSolver::MHDSolver::nVar];

  // Caveat: I don't properly use initialdatabyexahypespecfile here as it doesn't pass
  // the time. However, the Alfen Wave is the only exact solution we currently have, anyway

  alfenwave_(xpos, Conserved, &timeStamp);
  int error;
  pdecons2prim_(outputQuantities, Conserved, &error);
}


