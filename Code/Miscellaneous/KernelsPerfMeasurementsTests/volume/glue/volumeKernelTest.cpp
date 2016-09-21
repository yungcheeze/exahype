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

// g++ -std=c++11 -Wall -O3 -DALIGNMENT=64 -Drestrict=__restrict__ -DNDEBUG -march=native volumeKernelTest.cpp
// ---------------
/* 
icpc -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=O2 -restrict -O2 -DNDEBUG -xHost volumeKernelTest.cpp -o test_O2
icpc -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=O3 -restrict -O3 -DNDEBUG -xHost volumeKernelTest.cpp -o test_O3
icpc -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=O2_fast -restrict -O2 -fast -DNDEBUG -xHost volumeKernelTest.cpp -o test_O2_fast
icpc -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=O3_fast -restrict -O3 -fast -DNDEBUG -xHost volumeKernelTest.cpp -o test_O3_fast 
g++  -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=GCC_O2 -Drestrict=__restrict__  -O2 -DNDEBUG -march=native volumeKernelTest.cpp -o test_gcc_O2
g++  -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=GCC_O3 -Drestrict=__restrict__  -O3 -DNDEBUG -march=native volumeKernelTest.cpp -o test_gcc_O3
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <fstream>
#include <mm_malloc.h> //g++

#if defined _OPENMP
#include <omp.h>
#elif defined MPI_INT
#include <mpi.h>
#else
#include <sys/time.h>
#include <time.h>
#endif

#include <ctime>

#include "matrices.h"
#include "optimizedKernels.h"

#include "volumeIntegralNonlinear.cpp"

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)
#define FILENAMESUFFIX STR(FILESUFFIX)

double* Kxi_T;
double** Kxi; //=KxiGeneric[order]
double* s_m;
double* weights2;

using namespace std;

static volatile int iTest;

int main(int argc, char *argv[]) {
    #ifdef __ICC
    cout << "ICC: " << __ICC/100. << endl;
    #else
      #ifdef __GNUC__
      cout << "GCC: " << __GNUC__ << "." <<  __GNUC_MINOR__  << endl;
      #endif
    #endif

    ofstream resultsFile;
    string filename = "results_";
    filename += FILENAMESUFFIX;
    filename += ".dat";
    resultsFile.open(filename);

    clock_t timeStart, timeEnd;
    double t1, t2;
    const int nTests = 1000000;

    std::cout.precision(4);

	constexpr int maxOrder = 8;
	constexpr int nDOF = maxOrder+1;
	constexpr int nVar = NVAR;
  constexpr int nVarPadded = NVARPAD;

	double *lduh = new double[nVar*nDOF*nDOF*nDOF];
	double *lFhi = new double[nVarPadded*nDOF*nDOF*nDOF*nDOF*3];
	memset(lduh, 1.5, nVar*nDOF*nDOF*nDOF*sizeof(double));
	memset(lFhi, 1.5, nVarPadded*nDOF*nDOF*nDOF*nDOF*3*sizeof(double));
	
	


	const double dx[3] = {0.1,0.1,0.1};	

  initDGMatricesGeneric();
    
    cout << "Optimised Kernel order = 3" << endl;
    initDGMatrices3();
    initGaussLegendreNodesAndWeights3();
    timeStart = clock();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegral3(lduh, lFhi, dx);
    }
    timeEnd = clock();
    freeDGMatrices3();
    freeGaussLegendreNodesAndWeights3();
    t1 = double(timeEnd-timeStart) / CLOCKS_PER_SEC;
    //cout << t1 << endl; // in seconds

    cout << "Generic Kernel order = 3" << endl;
    Kxi = KxiGeneric[3];
    timeStart = clock();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegralNonlinear(lduh, lFhi, dx, nVar, 4);
    }
    timeEnd = clock();
    t2 = double(timeEnd-timeStart) / CLOCKS_PER_SEC; 
    //cout << t2 << endl;

    //        order | optimised | generic
    resultsFile << "3 "<<t1<<" "<<t2 << endl;

    // ---------------------------------------------

    cout << "Optimised Kernel order = 4" << endl;
    initDGMatrices4();
    initGaussLegendreNodesAndWeights4();
    timeStart = clock();
    for(iTest=0;iTest<nTests;iTest++) {
      volumeIntegral4(lduh, lFhi, dx);
    }
    timeEnd = clock();
    t1 = double(timeEnd-timeStart) / CLOCKS_PER_SEC; 

    cout << "Generic Kernel order = 4" << endl;
    Kxi = KxiGeneric[4];
    timeStart = clock();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegralNonlinear(lduh, lFhi, dx, nVar, 5);
    }
    timeEnd = clock();
    freeDGMatrices4();
    freeGaussLegendreNodesAndWeights4();
    t2 = double(timeEnd-timeStart) / CLOCKS_PER_SEC;

    
    //        order | optimised | generic
    resultsFile << "4 "<<t1<<" "<<t2 << endl;

    // ---------------------------------------------

    cout << "Optimised Kernel order = 5" << endl;
    initDGMatrices5();
    initGaussLegendreNodesAndWeights5();
    timeStart = clock();
    for(iTest=0;iTest<nTests;iTest++) {
      volumeIntegral5(lduh, lFhi, dx);
    }
    timeEnd = clock();
    freeDGMatrices5();
    freeGaussLegendreNodesAndWeights5();
    t1 = double(timeEnd-timeStart) / CLOCKS_PER_SEC; 

    cout << "Generic Kernel order = 5" << endl;
    Kxi = KxiGeneric[5];
    timeStart = clock();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegralNonlinear(lduh, lFhi, dx, nVar, 6);
    }
    timeEnd = clock();
    t2 = double(timeEnd-timeStart) / CLOCKS_PER_SEC;
    
    //        order | optimised | generic
    resultsFile << "5 "<<t1<<" "<<t2 << endl;	

    // ---------------------------------------------

    cout << "Optimised Kernel order = 6" << endl;
    initDGMatrices6();
    initGaussLegendreNodesAndWeights6();
    timeStart = clock();
    for(iTest=0;iTest<nTests;iTest++) {
      volumeIntegral6(lduh, lFhi, dx);
    }
    timeEnd = clock();
    freeDGMatrices6();
    freeGaussLegendreNodesAndWeights6();
    t1 = double(timeEnd-timeStart) / CLOCKS_PER_SEC;
    

    cout << "Generic Kernel order = 6" << endl;
    Kxi = KxiGeneric[6];
    timeStart = clock();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegralNonlinear(lduh, lFhi, dx, nVar, 7);
    }
    timeEnd = clock();
    t2 = double(timeEnd-timeStart) / CLOCKS_PER_SEC;
    
    //        order | optimised | generic
    resultsFile << "6 "<<t1<<" "<<t2 << endl;	

    // ---------------------------------------------

    cout << "Optimised Kernel order = 7" << endl;
    initDGMatrices7();
    initGaussLegendreNodesAndWeights7();
    timeStart = clock();
    for(iTest=0;iTest<nTests;iTest++) {
      volumeIntegral7(lduh, lFhi, dx);
    }
    timeEnd = clock();
    freeDGMatrices7();
    freeGaussLegendreNodesAndWeights7();
    t1 = double(timeEnd-timeStart) / CLOCKS_PER_SEC;

    cout << "Generic Kernel order = 7" << endl;
    Kxi = KxiGeneric[7];
    timeStart = clock();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegralNonlinear(lduh, lFhi, dx, nVar, 8);
    }
    timeEnd = clock();
    t2 = double(timeEnd-timeStart) / CLOCKS_PER_SEC;
    
    //        order | optimised | generic
    resultsFile << "7 "<<t1<<" "<<t2 << endl;	

    // ---------------------------------------------

    cout << "Optimised Kernel order = 8" << endl;
    initDGMatrices8();
    initGaussLegendreNodesAndWeights8();
    timeStart = clock();
    for(iTest=0;iTest<nTests;iTest++) {
      volumeIntegral8(lduh, lFhi, dx);
    }
    timeEnd = clock();
    freeDGMatrices8();
    freeGaussLegendreNodesAndWeights8();
    t1 = double(timeEnd-timeStart) / CLOCKS_PER_SEC;


    cout << "Generic Kernel order = 8" << endl;
    Kxi = KxiGeneric[8];
    timeStart = clock();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegralNonlinear(lduh, lFhi, dx, nVar, 9);
    }
    timeEnd = clock();
    t2 = double(timeEnd-timeStart) / CLOCKS_PER_SEC;

    
    //        order | optimised | generic
    resultsFile << "8 "<<t1<<" "<<t2 << endl;

    
	resultsFile.close();	
  
  
  freeDGMatricesGeneric();
  
	// clean up
	delete[] lduh;
	delete[] lFhi;

	return 0;    
}

/*
void initMatrices() {
  s_m[0] = 0.000000000000e+00;
  s_m[1] = 0.000000000000e+00;
  s_m[2] = 0.000000000000e+00;
  s_m[3] = 0.000000000000e+00;
  s_m[4] = 0.000000000000e+00;
  s_m[5] = 0.000000000000e+00;
  s_m[6] = 0.000000000000e+00;
  s_m[7] = 0.000000000000e+00;
  s_m[8] = 0.000000000000e+00;
  s_m[9] = 0.000000000000e+00;
  s_m[10] = 0.000000000000e+00;
  s_m[11] = 0.000000000000e+00;
  s_m[12] = 0.000000000000e+00;
  s_m[13] = 0.000000000000e+00;
  s_m[14] = 0.000000000000e+00;
  s_m[15] = 0.000000000000e+00;
  s_m[16] = 0.000000000000e+00;
  s_m[17] = 0.000000000000e+00;
  s_m[18] = 0.000000000000e+00;
  s_m[19] = 0.000000000000e+00;
  s_m[20] = 0.000000000000e+00;
  s_m[21] = 0.000000000000e+00;
  s_m[22] = 0.000000000000e+00;
  s_m[23] = 0.000000000000e+00;
  s_m[24] = 0.000000000000e+00;
  s_m[25] = 0.000000000000e+00;
  s_m[26] = 0.000000000000e+00;
  s_m[27] = 0.000000000000e+00;
  s_m[28] = 0.000000000000e+00;
  s_m[29] = 0.000000000000e+00;
  s_m[30] = 0.000000000000e+00;
  s_m[31] = 0.000000000000e+00;
  s_m[32] = 0.000000000000e+00;
  s_m[33] = 0.000000000000e+00;
  s_m[34] = 0.000000000000e+00;
  s_m[35] = 0.000000000000e+00;
  s_m[36] = 0.000000000000e+00;
  s_m[37] = 0.000000000000e+00;
  s_m[38] = 0.000000000000e+00;
  s_m[39] = 0.000000000000e+00;
  s_m[40] = 0.000000000000e+00;
  s_m[41] = 0.000000000000e+00;
  s_m[42] = 0.000000000000e+00;
  s_m[43] = 0.000000000000e+00;
  s_m[44] = 0.000000000000e+00;
  s_m[45] = 0.000000000000e+00;
  s_m[46] = 0.000000000000e+00;
  s_m[47] = 0.000000000000e+00;
  s_m[48] = 0.000000000000e+00;
  s_m[49] = 0.000000000000e+00;
  s_m[50] = 0.000000000000e+00;
  s_m[51] = 0.000000000000e+00;
  s_m[52] = 0.000000000000e+00;
  s_m[53] = 0.000000000000e+00;
  s_m[54] = 0.000000000000e+00;
  s_m[55] = 0.000000000000e+00;
  s_m[56] = 0.000000000000e+00;
  s_m[57] = 0.000000000000e+00;
  s_m[58] = 0.000000000000e+00;
  s_m[59] = 0.000000000000e+00;
  s_m[60] = 0.000000000000e+00;
  s_m[61] = 0.000000000000e+00;
  s_m[62] = 0.000000000000e+00;
  s_m[63] = 0.000000000000e+00;
  s_m[64] = 0.000000000000e+00;
  s_m[65] = 0.000000000000e+00;
  s_m[66] = 0.000000000000e+00;
  s_m[67] = 0.000000000000e+00;
  s_m[68] = 0.000000000000e+00;
  s_m[69] = 0.000000000000e+00;
  s_m[70] = 0.000000000000e+00;
  s_m[71] = 0.000000000000e+00;
  s_m[72] = 0.000000000000e+00;
  s_m[73] = 0.000000000000e+00;
  s_m[74] = 0.000000000000e+00;
  s_m[75] = 0.000000000000e+00;
  s_m[76] = 0.000000000000e+00;
  s_m[77] = 0.000000000000e+00;
  s_m[78] = 0.000000000000e+00;
  s_m[79] = 0.000000000000e+00;
  s_m[80] = 0.000000000000e+00;
  s_m[81] = 0.000000000000e+00;
  s_m[82] = 0.000000000000e+00;
  s_m[83] = 0.000000000000e+00;
  s_m[84] = 0.000000000000e+00;
  s_m[85] = 0.000000000000e+00;
  s_m[86] = 0.000000000000e+00;
  s_m[87] = 0.000000000000e+00;
  s_m[88] = 0.000000000000e+00;
  s_m[89] = 0.000000000000e+00;
  s_m[90] = 0.000000000000e+00;
  s_m[91] = 0.000000000000e+00;
  s_m[92] = 0.000000000000e+00;
  s_m[93] = 0.000000000000e+00;
  s_m[94] = 0.000000000000e+00;
  s_m[95] = 0.000000000000e+00;
  s_m[96] = 0.000000000000e+00;
  s_m[97] = 0.000000000000e+00;
  s_m[98] = 0.000000000000e+00;
  s_m[99] = 0.000000000000e+00;
  s_m[100] = 0.000000000000e+00;
  s_m[101] = 0.000000000000e+00;
  s_m[102] = 0.000000000000e+00;
  s_m[103] = 0.000000000000e+00;
  s_m[104] = 0.000000000000e+00;
  s_m[105] = 0.000000000000e+00;
  s_m[106] = 0.000000000000e+00;
  s_m[107] = 0.000000000000e+00;
  Kxi_T[0] = -1.255656088098e+00;
  Kxi_T[1] = -4.183978248600e-01;
  Kxi_T[2] = 1.300181658087e-01;
  Kxi_T[3] = -6.548023049355e-02;
  Kxi_T[4] = 4.236031185339e-02;
  Kxi_T[5] = -3.262358480478e-02;
  Kxi_T[6] = 2.916728809238e-02;
  Kxi_T[7] = -3.064117428684e-02;
  Kxi_T[8] = 4.197362432634e-02;
  Kxi_T[9] = 0.000000000000e+00;
  Kxi_T[10] = 0.000000000000e+00;
  Kxi_T[11] = 0.000000000000e+00;
  Kxi_T[12] = 2.010021351292e+00;
  Kxi_T[13] = -5.016657854850e-01;
  Kxi_T[14] = -6.769832419076e-01;
  Kxi_T[15] = 2.692138635353e-01;
  Kxi_T[16] = -1.602985146075e-01;
  Kxi_T[17] = 1.187447077945e-01;
  Kxi_T[18] = -1.039993129635e-01;
  Kxi_T[19] = 1.080391382041e-01;
  Kxi_T[20] = -1.472029988823e-01;
  Kxi_T[21] = 0.000000000000e+00;
  Kxi_T[22] = 0.000000000000e+00;
  Kxi_T[23] = 0.000000000000e+00;
  Kxi_T[24] = -1.294202075032e+00;
  Kxi_T[25] = 1.402699759000e+00;
  Kxi_T[26] = -2.562639281818e-01;
  Kxi_T[27] = -8.239035397637e-01;
  Kxi_T[28] = 3.777470995127e-01;
  Kxi_T[29] = -2.540518734313e-01;
  Kxi_T[30] = 2.124411755333e-01;
  Kxi_T[31] = -2.154851142534e-01;
  Kxi_T[32] = 2.903314666640e-01;
  Kxi_T[33] = 0.000000000000e+00;
  Kxi_T[34] = 0.000000000000e+00;
  Kxi_T[35] = 0.000000000000e+00;
  Kxi_T[36] = 9.350501663232e-01;
  Kxi_T[37] = -8.002227249223e-01;
  Kxi_T[38] = 1.181960666698e+00;
  Kxi_T[39] = -1.131793140867e-01;
  Kxi_T[40] = -9.369707771701e-01;
  Kxi_T[41] = 4.816403690687e-01;
  Kxi_T[42] = -3.644593173890e-01;
  Kxi_T[43] = 3.529618140521e-01;
  Kxi_T[44] = -4.658610418418e-01;
  Kxi_T[45] = 0.000000000000e+00;
  Kxi_T[46] = 0.000000000000e+00;
  Kxi_T[47] = 0.000000000000e+00;
  Kxi_T[48] = -6.759723042189e-01;
  Kxi_T[49] = 5.324612399834e-01;
  Kxi_T[50] = -6.055818291441e-01;
  Kxi_T[51] = 1.047058839436e+00;
  Kxi_T[52] = 4.399672026717e-16;
  Kxi_T[53] = -1.047058839436e+00;
  Kxi_T[54] = 6.055818291441e-01;
  Kxi_T[55] = -5.324612399834e-01;
  Kxi_T[56] = 6.759723042189e-01;
  Kxi_T[57] = 0.000000000000e+00;
  Kxi_T[58] = 0.000000000000e+00;
  Kxi_T[59] = 0.000000000000e+00;
  Kxi_T[60] = 4.658610418418e-01;
  Kxi_T[61] = -3.529618140521e-01;
  Kxi_T[62] = 3.644593173890e-01;
  Kxi_T[63] = -4.816403690687e-01;
  Kxi_T[64] = 9.369707771701e-01;
  Kxi_T[65] = 1.131793140867e-01;
  Kxi_T[66] = -1.181960666698e+00;
  Kxi_T[67] = 8.002227249223e-01;
  Kxi_T[68] = -9.350501663232e-01;
  Kxi_T[69] = 0.000000000000e+00;
  Kxi_T[70] = 0.000000000000e+00;
  Kxi_T[71] = 0.000000000000e+00;
  Kxi_T[72] = -2.903314666640e-01;
  Kxi_T[73] = 2.154851142534e-01;
  Kxi_T[74] = -2.124411755333e-01;
  Kxi_T[75] = 2.540518734313e-01;
  Kxi_T[76] = -3.777470995127e-01;
  Kxi_T[77] = 8.239035397637e-01;
  Kxi_T[78] = 2.562639281818e-01;
  Kxi_T[79] = -1.402699759000e+00;
  Kxi_T[80] = 1.294202075032e+00;
  Kxi_T[81] = 0.000000000000e+00;
  Kxi_T[82] = 0.000000000000e+00;
  Kxi_T[83] = 0.000000000000e+00;
  Kxi_T[84] = 1.472029988823e-01;
  Kxi_T[85] = -1.080391382041e-01;
  Kxi_T[86] = 1.039993129635e-01;
  Kxi_T[87] = -1.187447077945e-01;
  Kxi_T[88] = 1.602985146075e-01;
  Kxi_T[89] = -2.692138635353e-01;
  Kxi_T[90] = 6.769832419076e-01;
  Kxi_T[91] = 5.016657854850e-01;
  Kxi_T[92] = -2.010021351292e+00;
  Kxi_T[93] = 0.000000000000e+00;
  Kxi_T[94] = 0.000000000000e+00;
  Kxi_T[95] = 0.000000000000e+00;
  Kxi_T[96] = -4.197362432634e-02;
  Kxi_T[97] = 3.064117428684e-02;
  Kxi_T[98] = -2.916728809238e-02;
  Kxi_T[99] = 3.262358480478e-02;
  Kxi_T[100] = -4.236031185339e-02;
  Kxi_T[101] = 6.548023049355e-02;
  Kxi_T[102] = -1.300181658087e-01;
  Kxi_T[103] = 4.183978248600e-01;
  Kxi_T[104] = 1.255656088098e+00;
  Kxi_T[105] = 0.000000000000e+00;
  Kxi_T[106] = 0.000000000000e+00;
  Kxi_T[107] = 0.000000000000e+00;
  weights2[0] = 1.651381550887e-03;
  weights2[1] = 3.670517192280e-03;
  weights2[2] = 5.295243737658e-03;
  weights2[3] = 6.346454410738e-03;
  weights2[4] = 6.710000397662e-03;
  weights2[5] = 6.346454410738e-03;
  weights2[6] = 5.295243737658e-03;
  weights2[7] = 3.670517192280e-03;
  weights2[8] = 1.651381550887e-03;
  weights2[9] = 3.670517192280e-03;
  weights2[10] = 8.158439490609e-03;
  weights2[11] = 1.176971074065e-02;
  weights2[12] = 1.410623124142e-02;
  weights2[13] = 1.491428301751e-02;
  weights2[14] = 1.410623124142e-02;
  weights2[15] = 1.176971074065e-02;
  weights2[16] = 8.158439490609e-03;
  weights2[17] = 3.670517192280e-03;
  weights2[18] = 5.295243737658e-03;
  weights2[19] = 1.176971074065e-02;
  weights2[20] = 1.697948376991e-02;
  weights2[21] = 2.035024731670e-02;
  weights2[22] = 2.151597707163e-02;
  weights2[23] = 2.035024731670e-02;
  weights2[24] = 1.697948376991e-02;
  weights2[25] = 1.176971074065e-02;
  weights2[26] = 5.295243737658e-03;
  weights2[27] = 6.346454410738e-03;
  weights2[28] = 1.410623124142e-02;
  weights2[29] = 2.035024731670e-02;
  weights2[30] = 2.439017413386e-02;
  weights2[31] = 2.578732431455e-02;
  weights2[32] = 2.439017413386e-02;
  weights2[33] = 2.035024731670e-02;
  weights2[34] = 1.410623124142e-02;
  weights2[35] = 6.346454410738e-03;
  weights2[36] = 6.710000397662e-03;
  weights2[37] = 1.491428301751e-02;
  weights2[38] = 2.151597707163e-02;
  weights2[39] = 2.578732431455e-02;
  weights2[40] = 2.726450789791e-02;
  weights2[41] = 2.578732431455e-02;
  weights2[42] = 2.151597707163e-02;
  weights2[43] = 1.491428301751e-02;
  weights2[44] = 6.710000397662e-03;
  weights2[45] = 6.346454410738e-03;
  weights2[46] = 1.410623124142e-02;
  weights2[47] = 2.035024731670e-02;
  weights2[48] = 2.439017413386e-02;
  weights2[49] = 2.578732431455e-02;
  weights2[50] = 2.439017413386e-02;
  weights2[51] = 2.035024731670e-02;
  weights2[52] = 1.410623124142e-02;
  weights2[53] = 6.346454410738e-03;
  weights2[54] = 5.295243737658e-03;
  weights2[55] = 1.176971074065e-02;
  weights2[56] = 1.697948376991e-02;
  weights2[57] = 2.035024731670e-02;
  weights2[58] = 2.151597707163e-02;
  weights2[59] = 2.035024731670e-02;
  weights2[60] = 1.697948376991e-02;
  weights2[61] = 1.176971074065e-02;
  weights2[62] = 5.295243737658e-03;
  weights2[63] = 3.670517192280e-03;
  weights2[64] = 8.158439490609e-03;
  weights2[65] = 1.176971074065e-02;
  weights2[66] = 1.410623124142e-02;
  weights2[67] = 1.491428301751e-02;
  weights2[68] = 1.410623124142e-02;
  weights2[69] = 1.176971074065e-02;
  weights2[70] = 8.158439490609e-03;
  weights2[71] = 3.670517192280e-03;
  weights2[72] = 1.651381550887e-03;
  weights2[73] = 3.670517192280e-03;
  weights2[74] = 5.295243737658e-03;
  weights2[75] = 6.346454410738e-03;
  weights2[76] = 6.710000397662e-03;
  weights2[77] = 6.346454410738e-03;
  weights2[78] = 5.295243737658e-03;
  weights2[79] = 3.670517192280e-03;
  weights2[80] = 1.651381550887e-03;
  weights2[81] = 0.000000000000e+00;
  weights2[82] = 0.000000000000e+00;
  weights2[83] = 0.000000000000e+00;

  Kxi[0][0] = -1.260660205839e+00;
  Kxi[0][1] = -4.148927713020e-01;
  Kxi[0][2] = 1.269322840112e-01;
  Kxi[0][3] = -6.242556341322e-02;
  Kxi[0][4] = 3.904733259262e-02;
  Kxi[0][5] = -2.869242542062e-02;
  Kxi[0][6] = 2.397617139031e-02;
  Kxi[0][7] = -2.260981542902e-02;
  Kxi[0][8] = 2.455626988016e-02;
  Kxi[0][9] = -3.422882091683e-02;
  Kxi[1][0] = 2.027220013968e+00;
  Kxi[1][1] = -5.137174273279e-01;
  Kxi[1][2] = -6.663645647310e-01;
  Kxi[1][3] = 2.586902796970e-01;
  Kxi[1][4] = -1.488666530544e-01;
  Kxi[1][5] = 1.051510984479e-01;
  Kxi[1][6] = -8.600065857106e-02;
  Kxi[1][7] = 8.010053769945e-02;
  Kxi[1][8] = -8.638173508757e-02;
  Kxi[1][9] = 1.199851267913e-01;
  Kxi[2][0] = -1.326846946353e+00;
  Kxi[2][1] = 1.425593189709e+00;
  Kxi[2][2] = -2.764647873771e-01;
  Kxi[2][3] = -8.038396407681e-01;
  Kxi[2][4] = 3.558856451654e-01;
  Kxi[2][5] = -2.279530855410e-01;
  Kxi[2][6] = 1.777094388193e-01;
  Kxi[2][7] = -1.612329092336e-01;
  Kxi[2][8] = 1.713638255696e-01;
  Kxi[2][9] = -2.363446367747e-01;
  Kxi[3][0] = 9.843513264005e-01;
  Kxi[3][1] = -8.348390467020e-01;
  Kxi[3][2] = 1.212574259881e+00;
  Kxi[3][3] = -1.436881209529e-01;
  Kxi[3][4] = -9.035733595898e-01;
  Kxi[3][5] = 4.415232856090e-01;
  Kxi[3][6] = -3.106478783087e-01;
  Kxi[3][7] = 2.680707448618e-01;
  Kxi[3][8] = -2.775392562152e-01;
  Kxi[3][9] = 3.780658887103e-01;
  Kxi[4][0] = -7.413087738687e-01;
  Kxi[4][1] = 5.784158377192e-01;
  Kxi[4][2] = -6.463526073133e-01;
  Kxi[4][3] = 1.087886347304e+00;
  Kxi[4][4] = -4.499318311735e-02;
  Kxi[4][5] = -9.925290910988e-01;
  Kxi[4][6] = 5.315862285369e-01;
  Kxi[4][7] = -4.140039734282e-01;
  Kxi[4][8] = 4.085606779488e-01;
  Kxi[4][9] = -5.447221435018e-01;
  Kxi[5][0] = 5.447221435018e-01;
  Kxi[5][1] = -4.085606779488e-01;
  Kxi[5][2] = 4.140039734282e-01;
  Kxi[5][3] = -5.315862285369e-01;
  Kxi[5][4] = 9.925290910988e-01;
  Kxi[5][5] = 4.499318311735e-02;
  Kxi[5][6] = -1.087886347304e+00;
  Kxi[5][7] = 6.463526073133e-01;
  Kxi[5][8] = -5.784158377192e-01;
  Kxi[5][9] = 7.413087738687e-01;
  Kxi[6][0] = -3.780658887103e-01;
  Kxi[6][1] = 2.775392562152e-01;
  Kxi[6][2] = -2.680707448618e-01;
  Kxi[6][3] = 3.106478783087e-01;
  Kxi[6][4] = -4.415232856090e-01;
  Kxi[6][5] = 9.035733595898e-01;
  Kxi[6][6] = 1.436881209529e-01;
  Kxi[6][7] = -1.212574259881e+00;
  Kxi[6][8] = 8.348390467020e-01;
  Kxi[6][9] = -9.843513264005e-01;
  Kxi[7][0] = 2.363446367747e-01;
  Kxi[7][1] = -1.713638255696e-01;
  Kxi[7][2] = 1.612329092336e-01;
  Kxi[7][3] = -1.777094388193e-01;
  Kxi[7][4] = 2.279530855410e-01;
  Kxi[7][5] = -3.558856451654e-01;
  Kxi[7][6] = 8.038396407681e-01;
  Kxi[7][7] = 2.764647873771e-01;
  Kxi[7][8] = -1.425593189709e+00;
  Kxi[7][9] = 1.326846946353e+00;
  Kxi[8][0] = -1.199851267913e-01;
  Kxi[8][1] = 8.638173508757e-02;
  Kxi[8][2] = -8.010053769945e-02;
  Kxi[8][3] = 8.600065857106e-02;
  Kxi[8][4] = -1.051510984479e-01;
  Kxi[8][5] = 1.488666530544e-01;
  Kxi[8][6] = -2.586902796970e-01;
  Kxi[8][7] = 6.663645647310e-01;
  Kxi[8][8] = 5.137174273279e-01;
  Kxi[8][9] = -2.027220013968e+00;
  Kxi[9][0] = 3.422882091683e-02;
  Kxi[9][1] = -2.455626988016e-02;
  Kxi[9][2] = 2.260981542902e-02;
  Kxi[9][3] = -2.397617139031e-02;
  Kxi[9][4] = 2.869242542062e-02;
  Kxi[9][5] = -3.904733259262e-02;
  Kxi[9][6] = 6.242556341322e-02;
  Kxi[9][7] = -1.269322840112e-01;
  Kxi[9][8] = 4.148927713020e-01;
  Kxi[9][9] = 1.260660205839e+00;
}
*/

