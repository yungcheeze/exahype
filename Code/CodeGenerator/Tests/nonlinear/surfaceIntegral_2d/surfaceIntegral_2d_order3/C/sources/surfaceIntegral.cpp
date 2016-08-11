// Solve the surface integral 

#include "string.h"
#include "Kernels.h"
#include "DGMatrices.h"
#include "GaussLegendreQuadrature.h"

void kernels::aderdg::optimised::surfaceIntegral( 
  double* restrict lduh, 
  const double* restrict const lFbnd, 
  const double* dx
) {
#ifdef __INTEL_COMPILER
  __assume_aligned(s_m, ALIGNMENT);
  __assume_aligned(FRCoeff, ALIGNMENT);
  __assume_aligned(FLCoeff, ALIGNMENT);
  __assume_aligned(weights1, ALIGNMENT);
  __assume_aligned(weights2, ALIGNMENT);
#endif
  double FRCoeff_s[4] __attribute__((aligned(ALIGNMENT)));
  double FLCoeff_s[4] __attribute__((aligned(ALIGNMENT)));
  // x-direction
  for(int iVar=0;iVar<5;iVar++) {
    for(int jk=0;jk<4;jk++) {
#pragma simd
  for(int it=0;it<4;it++) 
    FRCoeff_s[it] = lFbnd[20+iVar*4+jk] * FRCoeff[it];
#pragma simd
  for(int it=0;it<4;it++) 
    FLCoeff_s[it] = lFbnd[0+iVar*4+jk] * FLCoeff[it];
#pragma simd
    for(int it=0;it<4;it++)
      s_m[it] = weights2[jk]/dx[0]*(FRCoeff_s[it]-FLCoeff_s[it]);
  lduh[jk*20+0+iVar] -= s_m[0];
  lduh[jk*20+5+iVar] -= s_m[1];
  lduh[jk*20+10+iVar] -= s_m[2];
  lduh[jk*20+15+iVar] -= s_m[3];
    }
  }
  // y-direction
  for(int iVar=0;iVar<5;iVar++) {
    for(int i=0;i<4;i++) {
#pragma simd
  for(int it=0;it<4;it++) 
    FRCoeff_s[it] = lFbnd[60+iVar*4+i] * FRCoeff[it];
#pragma simd
  for(int it=0;it<4;it++) 
    FLCoeff_s[it] = lFbnd[40+iVar*4+i] * FLCoeff[it];
        double s = weights1[i]/dx[1];
#pragma simd
        for(int it=0;it<4;it++)
          s_m[it] = s*(FRCoeff_s[it]-FLCoeff_s[it]);
  lduh[0+i*5+iVar] -= s_m[0];
  lduh[20+i*5+iVar] -= s_m[1];
  lduh[40+i*5+iVar] -= s_m[2];
  lduh[60+i*5+iVar] -= s_m[3];
    }
  }
}