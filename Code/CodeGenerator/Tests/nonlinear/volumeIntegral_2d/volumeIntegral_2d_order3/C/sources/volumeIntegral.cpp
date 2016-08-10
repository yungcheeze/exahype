// Solve the volume integral 

#include "Kernels.h"
#include "DGMatrices.h"
#include "GaussLegendreQuadrature.h"
#include <cstring>
#include "asm_volumeIntegral.c"
#include <iostream>

using namespace std;

void kernels::aderdg::optimised::volumeIntegral( 
  double* restrict lduh, 
  const double* restrict const lFhi, 
  const double* dx
) {
  memset(lduh, 0, sizeof(lduh)*80);

  // Assume equispaced mesh, dx[0] == dx[1] == dx[2]
#pragma simd
  for(int it=0;it<16;it++) 
    s_m[it] = weights2[0]/dx[0] * Kxi_T[it];

    gemm_5_4_4_lduh_x(&lFhi[0],&s_m[0],&lduh[0]);
    gemm_5_4_4_lduh_y(&lFhi[128],&s_m[0],&lduh[0]);

#pragma simd
  for(int it=0;it<16;it++) 
    s_m[it] = weights2[1]/dx[0] * Kxi_T[it];

    gemm_5_4_4_lduh_x(&lFhi[32],&s_m[0],&lduh[20]);
    gemm_5_4_4_lduh_y(&lFhi[160],&s_m[0],&lduh[5]);

#pragma simd
  for(int it=0;it<16;it++) 
    s_m[it] = weights2[2]/dx[0] * Kxi_T[it];

    gemm_5_4_4_lduh_x(&lFhi[64],&s_m[0],&lduh[40]);
    gemm_5_4_4_lduh_y(&lFhi[192],&s_m[0],&lduh[10]);

#pragma simd
  for(int it=0;it<16;it++) 
    s_m[it] = weights2[3]/dx[0] * Kxi_T[it];

    gemm_5_4_4_lduh_x(&lFhi[96],&s_m[0],&lduh[60]);
    gemm_5_4_4_lduh_y(&lFhi[224],&s_m[0],&lduh[15]);

}
