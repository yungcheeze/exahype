#ifndef __TIMESERIES_EULER__
#define __TIMESERIES_EULER__

/**
 * TimeSeriesReductions is a C++ class to be used in UserOnTheFlyPostProcessing classes
 * to compute reductions of fields and write them out as ASCII files. The reductions are
 * 
 *   l1: L^1 norm of field (l1 = sum of all cells)
 *   l2: L^2 norm of field (l2 = sqrt(sum of cell[i]^2))
 *   max: Maximum value of field (equal to L^\inf norm)
 *   min: Minimum value of field
 *   avg: Average value of field
 * 
 * The code is not very enthusiastic, so the following shortcomings are known:
 * 
 *   1) Code will not work with several MPI ranks. This needs communication
 *      to determine a single value over all ranks.
 *
 * (c) ExaHyPE 2016
 *
 **/
#include <vector>
#include <array>
#include <stdio.h>
#include <errno.h>
#include <string.h>

#include "tarch/parallel/Node.h"

using namespace std;

class TimeSeriesReductions {
    enum index { tidx=0, time=1, l1=2, l2=3, max=4, min=5, avg=6, LEN=7 };
    const char * const colnames[LEN]  = { "plotindex ", "time ", "l1norm ", "l2norm ", "max ", "min ", "avg " };
    const char * const colformat[LEN] = { "%.0f\t",        "%e\t",     "%e\t",       "%e\t",       "%e\t",    "%e\t",    "%e\t"  };
    double data[LEN];
    int avgcnt;
    FILE* asc;

public:
    TimeSeriesReductions(const char* filename) {
	if(!shallIRun()) return;
	    
	asc = fopen(filename, "w");
	if(asc == NULL) {
		fprintf(stderr, "ASCII writer: Could not open output file '%s'\n", filename);//, strerror(errno));
		exit(-1);
	}
	
	// print the file header line
	// fprintf(asc, "# ");
	for(int i=0; i<LEN; i++) {
	fputs(colnames[i], asc);
	}
	fprintf(asc, "\n");

        // zero everything
        for(int i=0; i<LEN; i++)
            data[i] = 0;
    }

    void initRow(double current_time) {
	if(!shallIRun()) return;
        data[time] = current_time;
        data[l1] = data[l2] = 0;
        data[max] = 0;
        data[min] = std::numeric_limits<double>::max();
        data[avg] = 0;
        avgcnt = 0;
    }

    void addValue(double val, double dx) {
        if(!shallIRun()) return;
	// dx is the scaling (volume form) as computed by peano::la::volume,
        // compare the calculation in MyEulerSolver_Plotter1.cpp 
        data[l1] += abs(val) * dx;
        data[l2] += val * val * dx;
        data[max] = std::max( data[max], val);
        data[min] = std::min( data[min], val);
        data[avg] += val;
        avgcnt++;
    }

    void finishRow() {
        if(!shallIRun()) return;
        data[tidx]++;
        data[l2] = sqrt(data[l2]);
        data[avg] = data[avg] / avgcnt;
    }

    void writeRow() {
        if(!shallIRun()) return;
        finishRow();

	for(int i=0; i<LEN; i++) {
            fprintf(asc, colformat[i], data[i]);
        }
        fprintf(asc, "\n");
	fflush(asc); // write out this line immediately
    }
    
    inline bool shallIRun() {
#ifdef Parallel
	return tarch::parallel::Node::getInstance().isGlobalMaster();
#else
	return true;
#endif
    }
};


#endif /* __TIMESERIES_EULER__ */
