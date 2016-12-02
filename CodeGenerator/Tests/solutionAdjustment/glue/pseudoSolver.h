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

#ifndef CODEGENTEST_GEN_PSEUDOSOLVER_H_
#define CODEGENTEST_GEN_PSEUDOSOLVER_H_

class PseudoSolver {
	public: 
    static void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);
    
    int getNumberOfVariables();
    int getNodesPerCoordinateAxis();
};

#endif //CODEGENTEST_GEN_PSEUDOSOLVER_H_