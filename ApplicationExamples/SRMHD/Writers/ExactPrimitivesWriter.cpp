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
  // convert tarch::la::Vector to double*  
  double xpos[DIMENSIONS];
  double V[10];
  for(int i=0; i<DIMENSIONS; i++) xpos[i] = x[i];

  // Caveat: I don't properly use initialdatabyexahypespecfile here as it doesn't pass
  // the time. However, the Alfen Wave is the only exact solution we currently have, anyway

  alfenwave_(xpos, V, &timeStamp);
  int i;
  pdecons2prim_(outputQuantities, V, &i);
}


