#include "ErrorWriter.h"

#include "MHDSolver.h"
#include "InitialDataAdapter.h"


MHDSolver::ErrorWriter::ErrorWriter(MHDSolver&  solver) {
}


MHDSolver::ErrorWriter::~ErrorWriter() {
}


void MHDSolver::ErrorWriter::startPlotting(double time) {
}


void MHDSolver::ErrorWriter::finishPlotting() {
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
  const int nDim = MHDSolver::MHDSolver::nDim;
  const int nVar = MHDSolver::MHDSolver::nVar;

  // convert tarch::la::Vector to double*  
  double xpos[nDim];
  for(int i=0; i<nDim; i++) xpos[i] = x[i];

  alfenwave_(xpos, outputQuantities, &timeStamp);
}


