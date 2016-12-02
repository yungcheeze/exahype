/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include <cmath> 
#include "pseudoSolver.h"

void PseudoSolver::adjustedSolutionValues(const double* const x,
                                                  const double w,
                                                  const double t,
                                                  const double dt, double* Q) {

  double GAMMA = 1.4;
  Q[0] = 1.;
  Q[1] = 0.;
  Q[2] = 0.;
  Q[3] = 0.;
  Q[4] = 1. / (GAMMA -1) +
        std::exp(-((x[0] -0.5) *(x[0] -0.5) + (x[1] -0.5) *(x[1] -0.5) 
#if DIMENSIONS == 3     
          + (x[2] -0.5) *(x[2] -0.5)
#endif        
        ) /
        (0.5 *0.5
#if DIMENSIONS == 3
          *0.5
#endif        
        )) *
        1.0e-1;
}

int PseudoSolver::getNumberOfVariables() {
  return NVAR;
}

int PseudoSolver::getNodesPerCoordinateAxis() {
  return ORDER +1; //basisSize
}