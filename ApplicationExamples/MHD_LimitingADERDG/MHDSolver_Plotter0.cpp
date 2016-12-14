#include "MHDSolver_Plotter0.h"

#include "C2P-MHD.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

MHD::MHDSolver_Plotter0::MHDSolver_Plotter0(exahype::solvers::LimitingADERDGSolver&  solver) {
  // @todo Please insert your code here
}


MHD::MHDSolver_Plotter0::~MHDSolver_Plotter0() {
  // @todo Please insert your code here
}


void MHD::MHDSolver_Plotter0::startPlotting(double time) {
  // @todo Please insert your code here
}


void MHD::MHDSolver_Plotter0::finishPlotting() {
  // @todo Please insert your code here
}


void MHD::MHDSolver_Plotter0::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<9; i++){
    // Conserved variables plotter
    int ierr;
    pdecons2prim_(outputQuantities,Q,&ierr);
  }
}


