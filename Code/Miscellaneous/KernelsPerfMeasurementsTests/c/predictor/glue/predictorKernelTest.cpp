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


// icpc -std=c++11 -Wall -DALIGNMENT=64 -restrict -O3 -DNDEBUG -xHost predictorKernelTest.cpp

/* 
icpc -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=O2 -restrict -O2 -DNDEBUG -xHost predictorKernelTest.cpp -o test_O2
icpc -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=O3 -restrict -O3 -DNDEBUG -xHost predictorKernelTest.cpp -o test_O3
icpc -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=O2_fast -restrict -O2 -fast -DNDEBUG -xHost predictorKernelTest.cpp -o test_O2_fast
icpc -std=c++11 -Wall -DALIGNMENT=64 -DFILESUFFIX=O3_fast -restrict -O3 -fast -DNDEBUG -xHost predictorKernelTest.cpp -o test_O3_fast 
g++  -std=c++11 -DALIGNMENT=64 -DFILESUFFIX=GCC_O2 -Drestrict=__restrict__  -O2 -DNDEBUG -march=native predictorKernelTest.cpp -o test_gcc_O2
g++  -std=c++11 -DALIGNMENT=64 -DFILESUFFIX=GCC_O3 -Drestrict=__restrict__  -O3 -DNDEBUG -march=native predictorKernelTest.cpp -o test_gcc_O3
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
#elsec
#include <sys/time.h>
#include <time.h>
#endif


#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)
#define FILENAMESUFFIX STR(FILESUFFIX)

#include "matrices.h"

#include "optimizedKernels.h"

#include "spaceTimePredictorNonlinear.cpph" //generic C

//global
double* weights1;

//local
double *lQi;
double *lFi;
double *lQhi;
double *lFhi;

using namespace std;
using namespace std::chrono;

void allocLocal(int nVar, int nVarPadded, int nDOF) {
  lQi = new double[nVarPadded*nDOF*nDOF*nDOF*nDOF];
  lFi = new double[nVarPadded*nDOF*nDOF*nDOF*nDOF*4];
  lQhi = new double[nVarPadded*nDOF*nDOF*nDOF];
  lFhi = new double[nVarPadded*nDOF*nDOF*nDOF*4];
  memset(lQi, 1.5, nVarPadded*nDOF*nDOF*nDOF*nDOF*sizeof(double));
  memset(lFi, 1.5, nVarPadded*nDOF*nDOF*nDOF*nDOF*4*sizeof(double));
  memset(lQhi, 1.5, nVarPadded*nDOF*nDOF*nDOF*sizeof(double));
  memset(lFhi, 1.5, nVarPadded*nDOF*nDOF*nDOF*4*sizeof(double));
}

void freeLocal() {
  delete[] lQi;
  delete[] lFi;
  delete[] lQhi;
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
    int nDOF;
    constexpr int nVar = NVAR;
    constexpr int nVarPad = NVARPAD;

    double dt = 0.01;
    const double dx[3] = {0.1,0.1,0.1};

    

	// ---------------------------------------------

    nDOF = 4;
    cout << "Optimised Kernel order = 3" << endl;
    allocLocal(nVar, nVarPad, nDOF); 
    initGaussLegendreNodesAndWeights3();
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
      predictor3(lQhi, lFhi, lQi, lFi); 
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    //cout << t1 << endl; // in seconds


    cout << "Generic Kernel order = 3" << endl;
    allocLocal(nVar, nVarPad, nDOF);
    
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        aderPredictorNonlinear(lQi, lFi, nVar, nDOF, lQhi, &lFhi[0], &lFhi[nVarPad*nDOF*nDOF*nDOF], &lFhi[2*nVarPad*nDOF*nDOF*nDOF]);
    }
    timeEnd = high_resolution_clock::now();
    freeGaussLegendreNodesAndWeights3();
    freeLocal();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec; 
    //cout << t2 << endl;

    //        order | optimised | generic
    resultsFile << "3 "<<t1<<" "<<t2 << endl;

	// ---------------------------------------------

    nDOF = 5;
    cout << "Optimised Kernel order = 4" << endl;
    allocLocal(nVar, nVarPad, nDOF);
    initGaussLegendreNodesAndWeights4();
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        predictor4(lQhi, lFhi, lQi, lFi); 
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    //cout << t1 << endl; // in seconds


    cout << "Generic Kernel order = 4" << endl;
    allocLocal(nVar, nVarPad, nDOF);
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        aderPredictorNonlinear(lQi, lFi, nVar, nDOF, lQhi, &lFhi[0], &lFhi[nVarPad*nDOF*nDOF*nDOF], &lFhi[2*nVarPad*nDOF*nDOF*nDOF]);
    }
    timeEnd = high_resolution_clock::now();
    freeGaussLegendreNodesAndWeights4();
    freeLocal();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec; 
    //cout << t2 << endl;

    //        order | optimised | generic
    resultsFile << "4 "<<t1<<" "<<t2 << endl;

	// ---------------------------------------------

    nDOF = 6;
    cout << "Optimised Kernel order = 5" << endl;
    allocLocal(nVar, nVarPad, nDOF);
    initGaussLegendreNodesAndWeights5();
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        predictor5(lQhi, lFhi, lQi, lFi); 
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    //cout << t1 << endl; // in seconds


    cout << "Generic Kernel order = 5" << endl;
    allocLocal(nVar, nVarPad, nDOF);
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        aderPredictorNonlinear(lQi, lFi, nVar, nDOF, lQhi, &lFhi[0], &lFhi[nVarPad*nDOF*nDOF*nDOF], &lFhi[2*nVarPad*nDOF*nDOF*nDOF]);
    }
    timeEnd = high_resolution_clock::now();
    freeGaussLegendreNodesAndWeights5();
    freeLocal();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec; 
    //cout << t2 << endl;

    //        order | optimised | generic
    resultsFile << "5 "<<t1<<" "<<t2 << endl;

	// ---------------------------------------------

    nDOF = 7;
    cout << "Optimised Kernel order = 6" << endl;
    allocLocal(nVar, nVarPad, nDOF);
    initGaussLegendreNodesAndWeights6();
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        predictor6(lQhi, lFhi, lQi, lFi); 
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    //cout << t1 << endl; // in seconds


    cout << "Generic Kernel order = 6" << endl;
    allocLocal(nVar, nVarPad, nDOF);
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        aderPredictorNonlinear(lQi, lFi, nVar, nDOF, lQhi, &lFhi[0], &lFhi[nVarPad*nDOF*nDOF*nDOF], &lFhi[2*nVarPad*nDOF*nDOF*nDOF]);
    }
    timeEnd = high_resolution_clock::now();
    freeGaussLegendreNodesAndWeights6();
    freeLocal();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec; 
    //cout << t2 << endl;

    //        order | optimised | generic
    resultsFile << "6 "<<t1<<" "<<t2 << endl;

	// ---------------------------------------------

    nDOF = 8;
    cout << "Optimised Kernel order = 7" << endl;
    allocLocal(nVar, nVarPad, nDOF);
    initGaussLegendreNodesAndWeights7();
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        predictor7(lQhi, lFhi, lQi, lFi); 
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    //cout << t1 << endl; // in seconds


    cout << "Generic Kernel order = 7" << endl;
    allocLocal(nVar, nVarPad, nDOF);
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        aderPredictorNonlinear(lQi, lFi, nVar, nDOF, lQhi, &lFhi[0], &lFhi[nVarPad*nDOF*nDOF*nDOF], &lFhi[2*nVarPad*nDOF*nDOF*nDOF]);
    }
    timeEnd = high_resolution_clock::now();
    freeGaussLegendreNodesAndWeights7();
    freeLocal();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec; 
    //cout << t2 << endl;

    //        order | optimised | generic
    resultsFile << "7 "<<t1<<" "<<t2 << endl;

	// ---------------------------------------------

    nDOF = 9;
    cout << "Optimised Kernel order = 8" << endl;
    allocLocal(nVar, nVarPad, nDOF);
    initGaussLegendreNodesAndWeights8();
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        predictor8(lQhi, lFhi, lQi, lFi); 
    }
    timeEnd = high_resolution_clock::now();
    freeLocal();
    t1 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec;
    //cout << t1 << endl; // in seconds


    cout << "Generic Kernel order = 8" << endl;
    allocLocal(nVar, nVarPad, nDOF);
    timeStart = high_resolution_clock::now();
    for(int iTest=0;iTest<nTests;iTest++) {
        aderPredictorNonlinear(lQi, lFi, nVar, nDOF, lQhi, &lFhi[0], &lFhi[nVarPad*nDOF*nDOF*nDOF], &lFhi[2*nVarPad*nDOF*nDOF*nDOF]);
    }
    timeEnd = high_resolution_clock::now();
    freeGaussLegendreNodesAndWeights8();
    freeLocal();
    t2 = duration_cast<nanoseconds>(timeEnd-timeStart).count() / nano2sec; 
    //cout << t2 << endl;

    //        order | optimised | generic
    resultsFile << "8 "<<t1<<" "<<t2 << endl;

	// ---------------------------------------------

  resultsFile.close();

	return 0;  
}
