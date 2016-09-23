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

#ifndef _EXAHYPE_TEST_MATRICES_H_
#define _EXAHYPE_TEST_MATRICES_H_


//GaussLegendreQuadrature
//extern double **gaussLegendreNodes;
//extern double **gaussLegendreWeights;
extern double *weights1;
//extern double *weights2;
//extern double *weights3;


// for order 3..8

void initGaussLegendreNodesAndWeights3();
void freeGaussLegendreNodesAndWeights3();

void initGaussLegendreNodesAndWeights4();
void freeGaussLegendreNodesAndWeights4();

void initGaussLegendreNodesAndWeights5();
void freeGaussLegendreNodesAndWeights5();

void initGaussLegendreNodesAndWeights6();
void freeGaussLegendreNodesAndWeights6();

void initGaussLegendreNodesAndWeights7();
void freeGaussLegendreNodesAndWeights7();

void initGaussLegendreNodesAndWeights8();
void freeGaussLegendreNodesAndWeights8();

#endif