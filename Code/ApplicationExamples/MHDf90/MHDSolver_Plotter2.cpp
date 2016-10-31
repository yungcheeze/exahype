#include "MHDSolver_Plotter2.h"

#include "InitialDataAdapter.h"
#include "MHDSolver.h"


MHDSolver::MHDSolver_Plotter2::MHDSolver_Plotter2(MHDSolver&  solver) {
  // @todo Please insert your code here
}


MHDSolver::MHDSolver_Plotter2::~MHDSolver_Plotter2() {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter2::startPlotting(double time) {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter2::finishPlotting() {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter2::mapQuantities(
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

  // We do NOT output the primitive quantities here as in the original MHD application
  // but the conservatives for CONVENIENCE.

  alfenwave_(xpos, outputQuantities, &timeStamp);
}


