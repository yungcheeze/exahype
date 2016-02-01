#include "kernels/GaussLegendreQuadrature.h"


double** kernels::gaussLegendreWeights;

double** kernels::gaussLegendreNodes;

void kernels::initGaussLegendreNodesAndWeights() {
  const int MaxOrder = 9;

  gaussLegendreNodes   = new double*[MaxOrder];
  gaussLegendreWeights = new double*[MaxOrder];

  for (int i=0; i<MaxOrder; i++) {
    gaussLegendreNodes[i]   = new double[i+1];
    gaussLegendreWeights[i] = new double[i+1];
  }

  gaussLegendreWeights[0][0] = 1.0000000000000000;
  gaussLegendreNodes[0][0]   = 0.5000000000000000;

  gaussLegendreWeights[1][0] = 0.5000000000000000;
  gaussLegendreWeights[1][1] = 0.5000000000000000;
  gaussLegendreNodes  [1][0] = 0.2113248654051871;
  gaussLegendreNodes  [1][1] = 0.7886751345948129;

  gaussLegendreWeights[2][0] = 0.2777777777777778;
  gaussLegendreWeights[2][1] = 0.4444444444444444;
  gaussLegendreWeights[2][2] = 0.2777777777777778;
  gaussLegendreNodes  [2][0] = 0.1127016653792583;
  gaussLegendreNodes  [2][1] = 0.5000000000000000;
  gaussLegendreNodes  [2][2] = 0.8872983346207417;

  gaussLegendreWeights[3][0] = 0.1739274225687273;
  gaussLegendreWeights[3][1] = 0.3260725774312732;
  gaussLegendreWeights[3][2] = 0.3260725774312732;
  gaussLegendreWeights[3][3] = 0.1739274225687273;
  gaussLegendreNodes  [3][0] = 0.0694318442029737;
  gaussLegendreNodes  [3][1] = 0.3300094782075719;
  gaussLegendreNodes  [3][2] = 0.6699905217924281;
  gaussLegendreNodes  [3][3] = 0.9305681557970262;
}
