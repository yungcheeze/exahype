#include "ComputeGlobalIntegrals.h"
#include "Primitives.h"
#include "InitialData.h"

Euler::ComputeGlobalIntegrals::ComputeGlobalIntegrals(MyEulerSolver&  solver)
{
	// open all the reductions
	assert( 5 == nVar );
	
	conserved[0] = new TimeSeriesReductions("output/dens.asc");
	conserved[1] = new TimeSeriesReductions("output/sconx.asc");
	conserved[2] = new TimeSeriesReductions("output/scony.asc");
	conserved[3] = new TimeSeriesReductions("output/sconz.asc");
	conserved[4] = new TimeSeriesReductions("output/ener.asc");
	
	primitives[0] = new TimeSeriesReductions("output/rho.asc"); // V[0]=Q[0]
	primitives[1] = new TimeSeriesReductions("output/velx.asc");
	primitives[2] = new TimeSeriesReductions("output/vely.asc");
	primitives[3] = new TimeSeriesReductions("output/velz.asc");
	primitives[4] = new TimeSeriesReductions("output/press.asc");
	
	errors[0] = new TimeSeriesReductions("output/error-rho.asc");
	errors[1] = new TimeSeriesReductions("output/error-velx.asc");
	errors[2] = new TimeSeriesReductions("output/error-vely.asc");
	errors[3] = new TimeSeriesReductions("output/error-velz.asc");
	errors[4] = new TimeSeriesReductions("output/error-press.asc");
	
	statistics = new TimeSeriesReductions("output/volform.asc");
}


Euler::ComputeGlobalIntegrals::~ComputeGlobalIntegrals() {
	// delete all reductions.
	// @todo  the deletes are missing here
}


void Euler::ComputeGlobalIntegrals::startPlotting(double time) {
	this->time = time;
	for(int i=0; i<nVar; i++) {
		conserved[i]->initRow(time);
		primitives[i]->initRow(time);
		errors[i]->initRow(time);
	}
	statistics->initRow(time);
}


void Euler::ComputeGlobalIntegrals::finishPlotting() {
	for(int i=0; i<nVar; i++) {
		conserved[i]->writeRow();
		primitives[i]->writeRow();
		errors[i]->writeRow();
	}
	statistics->writeRow();
}


void Euler::ComputeGlobalIntegrals::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
	// make sure this plotter has no output associated
	assertion( outputQuantities == nullptr );

	const double NumberOfLagrangePointsPerAxis = Euler::MyEulerSolver::order + 1;
	//const double NumberOfUnknownsPerGridPoint = nVar;

	// volume form for integration
	double scaling = tarch::la::volume(sizeOfPatch* (1.0/NumberOfLagrangePointsPerAxis));
	statistics->addValue(scaling, 1);

	// reduce the conserved quantities
	for (int i=0; i<nVar; i++)
		conserved[i]->addValue( Q[i], scaling );

	// reduce the primitive quantities
	double V[nVar];
	cons2prim(V, Q);
	for(int i=0; i<nVar; i++)
		primitives[i]->addValue( V[i], scaling );

	// now do the convergence test, as we have exact initial data
	double ExactCons[nVar];
	double ExactPrim[nVar];
	const double *xpos = x.data();
	
	idfunc(xpos, ExactCons, time); // Sven, this returns the conserved quantities
	cons2prim(ExactPrim, ExactCons);

	double localError[nVar];
	for(int i=0; i<nVar; i++) {
		localError[i] = abs(V[i] - ExactPrim[i]);
		errors[i]->addValue( localError[i], scaling );
	}
}


