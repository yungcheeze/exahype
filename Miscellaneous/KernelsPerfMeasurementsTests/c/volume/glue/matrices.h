#ifndef _EXAHYPE_TEST_MATRICES_H_
#define _EXAHYPE_TEST_MATRICES_H_

//DGMatrices
//extern double *Kxi;
extern double *Kxi_T;
//extern double *iK1;
//extern double *dudx;
extern double *s_m;
//extern double *s_v;
//extern double *tmp_bnd;
//extern double *F0; 
//extern double *FLCoeff;
//extern double *FRCoeff;
//extern double ***equidistantGridProjector1d;
//extern double **** fineGridProjector1d;

//GaussLegendreQuadrature
//extern double **gaussLegendreNodes;
//extern double **gaussLegendreWeights;
//extern double *weights1;
extern double *weights2;
//extern double *weights3;

//DGMatricesGeneric
extern double*** KxiGeneric;

void initDGMatricesGeneric();
void freeDGMatricesGeneric();


// for order 3..8

void initDGMatrices3();
void freeDGMatrices3();
void initGaussLegendreNodesAndWeights3();
void freeGaussLegendreNodesAndWeights3();

void initDGMatrices4();
void freeDGMatrices4();
void initGaussLegendreNodesAndWeights4();
void freeGaussLegendreNodesAndWeights4();

void initDGMatrices5();
void freeDGMatrices5();
void initGaussLegendreNodesAndWeights5();
void freeGaussLegendreNodesAndWeights5();

void initDGMatrices6();
void freeDGMatrices6();
void initGaussLegendreNodesAndWeights6();
void freeGaussLegendreNodesAndWeights6();

void initDGMatrices7();
void freeDGMatrices7();
void initGaussLegendreNodesAndWeights7();
void freeGaussLegendreNodesAndWeights7();

void initDGMatrices8();
void freeDGMatrices8();
void initGaussLegendreNodesAndWeights8();
void freeGaussLegendreNodesAndWeights8();

#endif