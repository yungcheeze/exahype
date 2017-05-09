#ifndef __INITIAL_DATA_ADAPTER_GRMHD__
#define __INITIAL_DATA_ADAPTER_GRMHD__

// Interface to initial data
struct idobj {
	virtual void Interpolate(const double* x, double t, double* Q) = 0;
};

extern idobj *id; // a global accessible pointer

// a function in InitialData.cpp which prepares the ID, both accessible
// from a pure ADERDG, pure FV or limiting application
void prepare_id();

// C ID:
#include "InitialData/PizzaTOV.h"

// FORTRAN ID:
extern "C" {

// FORTRAN functions called by C
void initialdata_(const double* x, const double* const t, double* Q);

// only initialdata_ is used, no more.

void initialalfenwave_(const double* x,const double* const t,  double* Q);
void initialblast_(const double* x, const double* const t, double* Q);
void initialorsagtang_(const double* x,const double* const t,  double* Q);
void initialrotor_(const double* x, const double* const t, double* Q);

// Exact solutions in FORTRAN
//void alfenwave_(const double* x, double* Q, const double* /* scalar */ t);

void initialaccretiondisc_(const double* x, const double* const t,  double* Q);
void initialaccretiondisc3d_(const double* x, const double* const t, double* Q);

}/* extern "C" */

// Fortran'ish C++ access
#include <stdio.h>
struct fortranid : public idobj {
	void Interpolate(const double* x, double t, double* Q) {
		initialdata_(x, &t, Q);
	}
};


#endif /* __INITIAL_DATA_ADAPTER_GRMHD__ */
