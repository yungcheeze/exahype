#include "SRHDSolverFV_Plotter1.h"


SRHD::SRHDSolverFV_Plotter1::SRHDSolverFV_Plotter1(SRHDSolverFV&  solver) {
  // @todo Please insert your code here
}


SRHD::SRHDSolverFV_Plotter1::~SRHDSolverFV_Plotter1() {
  // @todo Please insert your code here
}


void SRHD::SRHDSolverFV_Plotter1::startPlotting(double time) {
  // @todo Please insert your code here
}


void SRHD::SRHDSolverFV_Plotter1::finishPlotting() {
  // @todo Please insert your code here
}


extern "C" {
// Fortran function, implemented in C2P-SRHD.f90
void pdecons2prim_(double* V, double* Q, int* iErr);
}


void SRHD::SRHDSolverFV_Plotter1::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  // SRHD_FV primitive variables
  int iErr;
  pdecons2prim_(outputQuantities, Q, &iErr);
}
