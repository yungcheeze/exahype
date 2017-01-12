#include "Writers/RelativeErrorWriter.h"

#include "MHDSolver.h"
#include "InitialDataAdapter.h"

#include <limits> // double max
#include <cmath> // C++11: std::signbit

MHDSolver::RelativeErrorWriter::RelativeErrorWriter(MHDSolver&  solver) {}
MHDSolver::RelativeErrorWriter::~RelativeErrorWriter() {}
void MHDSolver::RelativeErrorWriter::startPlotting(double time) {}
void MHDSolver::RelativeErrorWriter::finishPlotting() {}

// There is also http://stackoverflow.com/a/25909704 not relying on C++11
template<typename T>
T sign(T arg) {
  // Returns the signum of arg. Works especially with NaNs and zeros.
  // signbit(arg) is true if arg is negative
  return std::signbit(arg) ? T(-1) : T(+1);
}

void MHDSolver::RelativeErrorWriter::mapQuantities(
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
  double exact[nVar];
  for(int i=0; i<nDim; i++) xpos[i] = x[i];

  alfenwave_(xpos, exact, &timeStamp);
  
  // compute the difference
  for (int i=0; i<nVar; i++) {
    outputQuantities[i] = ( Q[i] - exact[i] ) / exact[i];
    // Make NaN to finite value. Basically works like numpy.nan_to_num.
    if(outputQuantities[i] != outputQuantities[i] || outputQuantities[i] == std::numeric_limits<double>::infinity()) {
        outputQuantities[i] = sign(outputQuantities[i]) * std::numeric_limits<double>::max();
    }
  }
}


