/**
 * This is the plotter abused to compute global spacetime integrals.
 * This is actually used to compare with the exact solution for
 * convergence tests.
 **/

#include "GRMHDBenchmarkNormWriter.h"
#include "Fortran/PDE.h"
#include "kernels/GaussLegendreQuadrature.h"
#include <cmath>

#include "kernels/aderdg/generic/c/sizes.cpph"

GRMHD::GRMHDBenchmarkNormWriter::GRMHDBenchmarkNormWriter(exahype::solvers::LimitingADERDGSolver&  solver)
	: GRMHDBenchmarkNormWriter() { plotForADERSolver = true; }

GRMHD::GRMHDBenchmarkNormWriter::GRMHDBenchmarkNormWriter(GRMHD::GRMHDSolver_ADERDG&  solver)
	: GRMHDBenchmarkNormWriter() { plotForADERSolver = true; }

GRMHD::GRMHDBenchmarkNormWriter::GRMHDBenchmarkNormWriter(GRMHD::GRMHDSolver_FV&  solver)
	: GRMHDBenchmarkNormWriter() { plotForADERSolver = false; }

GRMHD::GRMHDBenchmarkNormWriter::GRMHDBenchmarkNormWriter() :
	evolvedPrimitives("evolved-"),
	errorPrimitives("error-")
{
	evolvedPrimitives.add(0, "dens");
	evolvedPrimitives.add(1, "sconx");
	evolvedPrimitives.add(2, "scony");
	evolvedPrimitives.add(3, "sconz");
	evolvedPrimitives.add(4, "ener");
	evolvedPrimitives.add(5, "bx");
	evolvedPrimitives.add(6, "by");
	evolvedPrimitives.add(7, "bz");
	evolvedPrimitives.add(8, "psi");
	
	errorPrimitives.add(0, "rho");
	errorPrimitives.add(1, "velx");
	errorPrimitives.add(2, "vely");
	errorPrimitives.add(3, "velz");
	errorPrimitives.add(4, "press");
	errorPrimitives.add(5, "bx");
	errorPrimitives.add(6, "by");
	errorPrimitives.add(7, "bz");
	errorPrimitives.add(8, "psi");
}


GRMHD::GRMHDBenchmarkNormWriter::~GRMHDBenchmarkNormWriter() {
	// delete all the TimeSeriesReductions. Not that important
	// as this is a kind of a Singleton object.
}


void GRMHD::GRMHDBenchmarkNormWriter::startPlotting(double time) {
	evolvedPrimitives.startRow(time);
	errorPrimitives.startRow(time);
}


void GRMHD::GRMHDBenchmarkNormWriter::finishPlotting() {
	evolvedPrimitives.finishRow();
	errorPrimitives.finishRow();
}


void GRMHD::GRMHDBenchmarkNormWriter::mapQuantities(
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


	double dV;
	if(plotForADERSolver) {
		const int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
		dV = kernels::ADERDGVolume(order, sizeOfPatch, pos);
	} else {
		const int patchSize = GRMHD::AbstractGRMHDSolver_FV::PatchSize;
		dV = tarch::la::volume(sizeOfPatch)/patchSize; // correct is probably (patchSize+1)
	}

	// statistics.addValue(dV, 1);

	// reduce the primitive quantities
	double V[nVar];
	int err;
	pdecons2prim_(V, Q, &err);
	evolvedPrimitives.addValue(V, dV);

	// now do the convergence test, as we have exact initial data
	double ExactCons[nVar];
	double ExactPrim[nVar];
	const double *xpos = x.data();
	
	alfenwave_(xpos, &timeStamp, ExactCons);
	pdecons2prim_(ExactPrim, ExactCons, &err);
	
	double localError[nVar];
	for(int i=0; i<nVar; i++) {
		localError[i] = std::abs(V[i] - ExactPrim[i]);
	}
	
	errorPrimitives.addValue(localError, dV);
}


