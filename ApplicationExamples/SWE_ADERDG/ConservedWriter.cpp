#include "ConservedWriter.h"


SWE::ConservedWriter::ConservedWriter(MySWESolver&  solver) {
  // @todo Please insert your code here
}


SWE::ConservedWriter::~ConservedWriter() {
  // @todo Please insert your code here
}


void SWE::ConservedWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void SWE::ConservedWriter::finishPlotting() {
  // @todo Please insert your code here
}


void SWE::ConservedWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  assertion5( Q[0]>=0, Q[0], Q[1], Q[2], Q[3], Q[4] );
  assertion5( Q[1]>=0, Q[0], Q[1], Q[2], Q[3], Q[4] );
  assertion5( Q[2]>=0, Q[0], Q[1], Q[2], Q[3], Q[4] );
  assertion5( Q[3]>=0, Q[0], Q[1], Q[2], Q[3], Q[4] );
  assertion5( Q[4]>=0, Q[0], Q[1], Q[2], Q[3], Q[4] );

  assertion5( Q[0]<=4, Q[0], Q[1], Q[2], Q[3], Q[4] );
  assertion5( Q[1]<=4, Q[0], Q[1], Q[2], Q[3], Q[4] );
  assertion5( Q[2]<=4, Q[0], Q[1], Q[2], Q[3], Q[4] );
  assertion5( Q[3]<=4, Q[0], Q[1], Q[2], Q[3], Q[4] );
  assertion5( Q[4]<=4, Q[0], Q[1], Q[2], Q[3], Q[4] );

  for (int i=0; i<4; i++){
    outputQuantities[i] = Q[i];
  }
  outputQuantities[4] = Q[3] + Q[0];
}


