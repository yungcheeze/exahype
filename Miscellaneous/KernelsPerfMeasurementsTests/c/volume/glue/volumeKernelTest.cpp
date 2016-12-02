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
#include <chrono>

#if defined _OPENMP
#include <omp.h>
#elif defined MPI_INT
#include <mpi.h>
#else
#include <sys/time.h>
#include <time.h>
#endif

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

double *lduh;
double *lFhi;

using namespace std;
using namespace std::chrono;

static volatile int iTest;

void allocLocal(int nVar, int nVarPadded, int nDOF) {
  lduh = new double[nVar*nDOF*nDOF*nDOF];
	lFhi = new double[nVarPadded*nDOF*nDOF*nDOF*4];
	memset(lduh, 1.5, nVar*nDOF*nDOF*nDOF*sizeof(double));
	memset(lFhi, 1.5, nVarPadded*nDOF*nDOF*nDOF*4*sizeof(double));
}

void freeLocal() {
  delete[] lduh;
	delete[] lFhi;
}

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

    high_resolution_clock::time_point timeStart, timeEnd;
    double t1, t2;
    const double nano2sec = 1000000000.0;
    const int nTests = 1000000;
    
    cout << "nTests=" << nTests << endl;
    resultsFile << "nTests=" << nTests << endl;

	constexpr int maxOrder = 8;
	constexpr int nDOF = maxOrder+1;
	constexpr int nVar = NVAR;
  constexpr int nVarPadded = NVARPAD;
	

	const double dx[3] = {0.1,0.1,0.1};	

  initDGMatricesGeneric();
    
    cout << "Optimised Kernel order = 3" << endl;
    allocLocal(nVar, nVarPadded, 4);
    initDGMatrices3();
    initGaussLegendreNodesAndWeights3();
    timeStart = high_resolution_clock::now();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegral3(lduh, lFhi, dx);
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    freeDGMatrices3();
    freeGaussLegendreNodesAndWeights3();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    //cout << t1 << endl; // in seconds

    cout << "Generic Kernel order = 3" << endl;
    allocLocal(nVar, nVarPadded, 4);
    Kxi = KxiGeneric[3];
    timeStart = high_resolution_clock::now();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegralNonlinear(lduh, lFhi, dx, nVar, 4);
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec; 
    //cout << t2 << endl;

    //        order | optimised | generic
    resultsFile << "3 "<<t1<<" "<<t2 << endl;

    // ---------------------------------------------

    cout << "Optimised Kernel order = 4" << endl;
    allocLocal(nVar, nVarPadded, 5);
    initDGMatrices4();
    initGaussLegendreNodesAndWeights4();
    timeStart = high_resolution_clock::now();
    for(iTest=0;iTest<nTests;iTest++) {
      volumeIntegral4(lduh, lFhi, dx);
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    freeDGMatrices4();
    freeGaussLegendreNodesAndWeights4();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec; 

    cout << "Generic Kernel order = 4" << endl;
    allocLocal(nVar, nVarPadded, 5);
    Kxi = KxiGeneric[4];
    timeStart = high_resolution_clock::now();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegralNonlinear(lduh, lFhi, dx, nVar, 5);
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    
    //        order | optimised | generic
    resultsFile << "4 "<<t1<<" "<<t2 << endl;

    // ---------------------------------------------

    cout << "Optimised Kernel order = 5" << endl;
    allocLocal(nVar, nVarPadded, 6);
    initDGMatrices5();
    initGaussLegendreNodesAndWeights5();
    timeStart = high_resolution_clock::now();
    for(iTest=0;iTest<nTests;iTest++) {
      volumeIntegral5(lduh, lFhi, dx);
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    freeDGMatrices5();
    freeGaussLegendreNodesAndWeights5();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec; 

    cout << "Generic Kernel order = 5" << endl;
    allocLocal(nVar, nVarPadded, 6);
    Kxi = KxiGeneric[5];
    timeStart = high_resolution_clock::now();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegralNonlinear(lduh, lFhi, dx, nVar, 6);
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    
    //        order | optimised | generic
    resultsFile << "5 "<<t1<<" "<<t2 << endl;	

    // ---------------------------------------------

    cout << "Optimised Kernel order = 6" << endl;
    allocLocal(nVar, nVarPadded, 7);
    initDGMatrices6();
    initGaussLegendreNodesAndWeights6();
    timeStart = high_resolution_clock::now();
    for(iTest=0;iTest<nTests;iTest++) {
      volumeIntegral6(lduh, lFhi, dx);
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    freeDGMatrices6();
    freeGaussLegendreNodesAndWeights6();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    

    cout << "Generic Kernel order = 6" << endl;
    allocLocal(nVar, nVarPadded, 7);
    Kxi = KxiGeneric[6];
    timeStart = high_resolution_clock::now();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegralNonlinear(lduh, lFhi, dx, nVar, 7);
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    
    //        order | optimised | generic
    resultsFile << "6 "<<t1<<" "<<t2 << endl;	

    // ---------------------------------------------

    cout << "Optimised Kernel order = 7" << endl;
    allocLocal(nVar, nVarPadded, 8);
    initDGMatrices7();
    initGaussLegendreNodesAndWeights7();
    timeStart = high_resolution_clock::now();
    for(iTest=0;iTest<nTests;iTest++) {
      volumeIntegral7(lduh, lFhi, dx);
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    freeDGMatrices7();
    freeGaussLegendreNodesAndWeights7();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;

    cout << "Generic Kernel order = 7" << endl;
    allocLocal(nVar, nVarPadded, 8);
    Kxi = KxiGeneric[7];
    timeStart = high_resolution_clock::now();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegralNonlinear(lduh, lFhi, dx, nVar, 8);
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    
    //        order | optimised | generic
    resultsFile << "7 "<<t1<<" "<<t2 << endl;	

    // ---------------------------------------------

    cout << "Optimised Kernel order = 8" << endl;
    allocLocal(nVar, nVarPadded, 9);
    initDGMatrices8();
    initGaussLegendreNodesAndWeights8();
    timeStart = high_resolution_clock::now();
    for(iTest=0;iTest<nTests;iTest++) {
      volumeIntegral8(lduh, lFhi, dx);
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    freeDGMatrices8();
    freeGaussLegendreNodesAndWeights8();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;


    cout << "Generic Kernel order = 8" << endl;
    allocLocal(nVar, nVarPadded, 9);
    Kxi = KxiGeneric[8];
    timeStart = high_resolution_clock::now();
    for(iTest=0;iTest<nTests;iTest++) {
        volumeIntegralNonlinear(lduh, lFhi, dx, nVar, 9);
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    
    //        order | optimised | generic
    resultsFile << "8 "<<t1<<" "<<t2 << endl;

    
	resultsFile.close();	
  
  
  freeDGMatricesGeneric();

	return 0;    
}
