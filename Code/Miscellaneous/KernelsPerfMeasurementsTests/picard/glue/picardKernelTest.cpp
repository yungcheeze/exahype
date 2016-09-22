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


// g++ -std=c++11 -Wall -DALIGNMENT=64 volumeKernelTest.cpp
// g++ -std=c++11 -Wall -DALIGNMENT=64 -Drestrict=__restrict__ -DNDEBUG -march=native picardKernelTest.cpp
// icpc -std=c++11 -Wall -DALIGNMENT=64 -restrict -O3 -DNDEBUG -xHost picardKernelTest.cpp

/* 
icpc -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=O2 -restrict -O2 -DNDEBUG -xHost picardKernelTest.cpp -o test_O2
icpc -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=O3 -restrict -O3 -DNDEBUG -xHost picardKernelTest.cpp -o test_O3
icpc -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=O2_fast -restrict -O2 -fast -DNDEBUG -xHost picardKernelTest.cpp -o test_O2_fast
icpc -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=O3_fast -restrict -O3 -fast -DNDEBUG -xHost picardKernelTest.cpp -o test_O3_fast 
g++  -std=c++11 -DALIGNMENT=64 -DFILESUFFIX=GCC_O2 -Drestrict=__restrict__  -O2 -DNDEBUG -march=native picardKernelTest.cpp -o test_gcc_O2
g++  -std=c++11 -DALIGNMENT=64 -DFILESUFFIX=GCC_O3 -Drestrict=__restrict__  -O3 -DNDEBUG -march=native picardKernelTest.cpp -o test_gcc_O3
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

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)
#define FILENAMESUFFIX STR(FILESUFFIX)

#include "matrices.h"

#include "optimizedKernels.h"
#include "spaceTimePredictorNonlinear.cpph" //generic C

double* Kxi_T;
double** KXI;
double* Kxi;
double* s_m;
double* weights1;
double* weights2;
double *weights3;
double *F0;
double *iK1;

using namespace std;
using namespace std::chrono;

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
    const int nTests = 1000;
    
    cout << "nTests=" << nTests << endl;
    resultsFile << "nTests=" << nTests << endl;
    
    constexpr int maxOrder = 8;
    constexpr int nDOF = maxOrder+1;
    constexpr int nVar = NVAR;
    constexpr int nVarPad = NVARPAD;

    double *luh = new double[nVar*nDOF*nDOF*nDOF*nDOF];
    double *lQi = new double[nVarPad*nDOF*nDOF*nDOF*nDOF];
    double *lFi = new double[nVarPad*nDOF*nDOF*nDOF*nDOF*3*2];
    memset(luh, 1.5, nVar*nDOF*nDOF*nDOF*nDOF*sizeof(double));
    memset(lQi, 1.5, nVarPad*nDOF*nDOF*nDOF*nDOF*sizeof(double));
    memset(lFi, 1.5, nVarPad*nDOF*nDOF*nDOF*nDOF*3*2*sizeof(double));

    double dt = 0.01;
    const double dx[3] = {0.1,0.1,0.1};

    initDGMatricesGeneric();

    cout << "Optimised Kernel order = 3" << endl;
    initDGMatrices3();
    initGaussLegendreNodesAndWeights3();
    KXI = KxiGeneric[3];
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        picardLoop3(lQi, lFi, luh, dx, dt);
    }
    timeEnd = high_resolution_clock::now();
    freeDGMatrices3();
    freeGaussLegendreNodesAndWeights3();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;

    cout << "Generic Kernel order = 3" << endl;
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        aderPicardLoopNonlinear(luh, dt, &dx[0], nVar, 4, lQi, lFi);
    }
    timeEnd = high_resolution_clock::now();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec; 

    //        order | optimised | generic
    resultsFile << "3 "<<t1<<" "<<t2 << endl;

	// ---------------------------------------------

    cout << "Optimised Kernel order = 4" << endl;
    initDGMatrices4();
    initGaussLegendreNodesAndWeights4();
    KXI = KxiGeneric[4];
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        picardLoop4(lQi, lFi, luh, dx, dt);
    }
    timeEnd = high_resolution_clock::now();
    freeDGMatrices4();
    freeGaussLegendreNodesAndWeights4();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec; 

    cout << "Generic Kernel order = 4" << endl;
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        aderPicardLoopNonlinear(luh, dt, &dx[0], nVar, 5, lQi, lFi);
    }
    timeEnd = high_resolution_clock::now();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;


    //        order | optimised | generic
    resultsFile << "4 "<<t1<<" "<<t2 << endl;

    // ---------------------------------------------

    cout << "Optimised Kernel order = 5" << endl;
    initDGMatrices5();
    initGaussLegendreNodesAndWeights5();
    KXI = KxiGeneric[5];
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        picardLoop5(lQi, lFi, luh, dx, dt);
    }
    timeEnd = high_resolution_clock::now();
    freeDGMatrices5();
    freeGaussLegendreNodesAndWeights5();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec; 

    cout << "Generic Kernel order = 5" << endl;
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        aderPicardLoopNonlinear(luh, dt, &dx[0], nVar, 6, lQi, lFi);
    }
    timeEnd = high_resolution_clock::now();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    
    //        order | optimised | generic
    resultsFile << "5 "<<t1<<" "<<t2 << endl;	

	// ---------------------------------------------

    cout << "Optimised Kernel order = 6" << endl;
    initDGMatrices6();
    initGaussLegendreNodesAndWeights6();
    KXI = KxiGeneric[6];
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        picardLoop6(lQi, lFi, luh, dx, dt);
    }
    timeEnd = high_resolution_clock::now();
    freeDGMatrices6();
    freeGaussLegendreNodesAndWeights6();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;


    cout << "Generic Kernel order = 6" << endl;
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        aderPicardLoopNonlinear(luh, dt, &dx[0], nVar, 7, lQi, lFi);
    }
    timeEnd = high_resolution_clock::now();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;

    //        order | optimised | generic
    resultsFile << "6 "<<t1<<" "<<t2 << endl;	

	// ---------------------------------------------

    cout << "Optimised Kernel order = 7" << endl;
    initDGMatrices7();
    initGaussLegendreNodesAndWeights7();
    KXI = KxiGeneric[7];
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        picardLoop7(lQi, lFi, luh, dx, dt);
    }
    timeEnd = high_resolution_clock::now();
    freeDGMatrices7();
    freeGaussLegendreNodesAndWeights7();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;

    cout << "Generic Kernel order = 7" << endl;
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        aderPicardLoopNonlinear(luh, dt, &dx[0], nVar, 8, lQi, lFi);
    }
    timeEnd = high_resolution_clock::now();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    
    //        order | optimised | generic
    resultsFile << "7 "<<t1<<" "<<t2 << endl;	

	// ---------------------------------------------

    cout << "Optimised Kernel order = 8" << endl;
    initDGMatrices8();
    initGaussLegendreNodesAndWeights8();
    KXI = KxiGeneric[8];
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        picardLoop8(lQi, lFi, luh, dx, dt);
    }
    timeEnd = high_resolution_clock::now();
    freeDGMatrices8();
    freeGaussLegendreNodesAndWeights8();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;


    cout << "Generic Kernel order = 8" << endl;
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        aderPicardLoopNonlinear(luh, dt, &dx[0], nVar, 9, lQi, lFi);
    }
    timeEnd = high_resolution_clock::now();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;


    //        order | optimised | generic
    resultsFile << "8 "<<t1<<" "<<t2 << endl;

  resultsFile.close();

	// clean up
  freeDGMatricesGeneric();
	delete[] luh;
  delete[] lQi;
  delete[] lFi;

	return 0;    
}
