/*
#include "kernels/initialvalues/PDEInitialValues.h"

#include <cmath>

#define GAMMA 1.4

void exahype::pde::PDEInitialValue2d(const double x,const double y,const int nvar,double * value) {
  for (int n=0; n < 5; n++) {
    value[n] = 0;
  }
  value[0] = 1.;
  value[4] = 1./(GAMMA-1) + std::exp(-((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))/(0.05*0.05)) * 1.0e-3;
}
*/
