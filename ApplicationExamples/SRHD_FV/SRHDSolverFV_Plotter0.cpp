#include "SRHDSolverFV_Plotter0.h"


SRHD::SRHDSolverFV_Plotter0::SRHDSolverFV_Plotter0(SRHDSolverFV&  solver) {
  // @todo Please insert your code here
}


SRHD::SRHDSolverFV_Plotter0::~SRHDSolverFV_Plotter0() {
  // @todo Please insert your code here
}


void SRHD::SRHDSolverFV_Plotter0::startPlotting(double time) {
  // @todo Please insert your code here
}


void SRHD::SRHDSolverFV_Plotter0::finishPlotting() {
  // @todo Please insert your code here
}


void SRHD::SRHDSolverFV_Plotter0::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<5; i++){
    // Conservative Quantities for SRHD_FV
    outputQuantities[i] = Q[i];
  }
}


