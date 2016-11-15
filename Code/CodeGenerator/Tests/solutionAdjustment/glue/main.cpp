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

#include <iostream>
#include <iomanip>
#include <cmath>
#include <set>

#include "kernels.h"
#include "pseudoSolver.h"
#include "../../Dependencies/Generic/GaussLegendreQuadrature.h"
#include "../../Dependencies/Optimised/GaussLegendreQuadrature.h"

using namespace std;

int main() {
  cout.precision(16);
  PseudoSolver* solver = new PseudoSolver();
  
  const int luhSize = NVAR*(ORDER+1)
#if DIMENSIONS == 3     
    *(ORDER+1)
#endif    
    *(ORDER+1);
    
  cout << "Parameters: dimension=" << DIMENSIONS << ", order=" << ORDER << ", nVar=" << NVAR << endl;
  
  double* luh_generic = new double[luhSize]();
  double* luh_optimised = new double[luhSize]();
  double* center = new double[DIMENSIONS];
  double* dx = new double[DIMENSIONS];
  double t = 0.0;
  double dt = 0.0;
  
  dx[0] = 1.0;
  center[0] = 0.5;
  dx[1] = 1.0;
  center[1] = 0.5;
#if DIMENSIONS == 3     
  dx[2] = 1.0;
  center[2] = 0.5;
#endif 

  std::set<int> orders;
  orders.insert(ORDER);

  kernels::initGaussLegendreNodesAndWeights(orders);
  kernels::aderdg::optimised::initGaussLegendreNodesAndWeights(orders);

  
  solutionAdjustment_generic<PseudoSolver>( *solver, luh_generic, center, dx, t, dt );
  solutionAdjustment_optimised<PseudoSolver::adjustedSolutionValues>( luh_optimised, center, dx, t, dt );

  double error = 0;
  for(int i=0; i<luhSize; i++) {
    error += fabs((luh_generic[i]-luh_optimised[i]));
  }
  
  cout << "Error: " << error << endl;
  
  cout << "First terms: " << endl << "i | generic | optimised" << endl;
  for(int i=0; i<10; i++) {
    cout << i << " | " << luh_generic[i] << " | "  << luh_optimised[i] << endl;
  }
 
  kernels::freeGaussLegendreNodesAndWeights(orders);
  kernels::aderdg::optimised::freeGaussLegendreNodesAndWeights(orders);
  
  delete solver;
  delete[] luh_generic;
  delete[] luh_optimised;
  delete[] center;
  delete[] dx;
}