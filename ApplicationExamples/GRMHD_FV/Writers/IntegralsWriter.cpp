/**
 * This is the plotter abused to compute global spacetime integrals.
 * This is actually used to compare with the exact solution for
 * convergence tests.
 **/

#include "Writers/IntegralsWriter.h"
#include "Writers/TimeSeriesReductions.h"
#include "Fortran/C2P-GRMHD.h"
#include "Fortran/InitialData.h"
#include "kernels/GaussLegendreQuadrature.h"
#include <cmath>

GRMHD::IntegralsWriter::IntegralsWriter(GRMHDSolver&  solver) :
	conserved("output/cons-"),
	primitives("output/prim-"),
	errors("output/error-"),
	statistics("output/volform.asc"),
	masschange("output/massdt.asc")
{
	conserved.add(0, "dens");
	conserved.add(1, "sconx");
	conserved.add(2, "scony");
	conserved.add(3, "sconz");
	conserved.add(4, "ener");
	conserved.add(5, "bx");
	conserved.add(6, "by");
	conserved.add(7, "bz");
	conserved.add(8, "psi");
		
	primitives.add(0, "rho"); // V[0]=Q[0]
	primitives.add(1, "velx");
	primitives.add(2, "vely");
	primitives.add(3, "velz");
	primitives.add(4, "press");
	primitives.add(5, "bx");
	primitives.add(6, "by");
	primitives.add(7, "bz");
	primitives.add(8, "psi");
	
	errors.add(0, "rho");
	errors.add(1, "velx");
	errors.add(2, "vely");
	errors.add(3, "velz");
	errors.add(4, "press");
	errors.add(5, "bx");
	errors.add(6, "by");
	errors.add(7, "bz");
	errors.add(8, "psi");
}


GRMHD::IntegralsWriter::~IntegralsWriter() {
	// delete all the TimeSeriesReductions. Not that important
	// as this is a kind of a Singleton object.
}


void GRMHD::IntegralsWriter::startPlotting(double time) {
	conserved.startRow(time);
	primitives.startRow(time);
	errors.startRow(time);
	statistics.startRow(time);
	masschange.startRow(time);
}


void GRMHD::IntegralsWriter::finishPlotting() {
	conserved.finishRow();
	primitives.finishRow();
	errors.finishRow();
	statistics.finishRow();
	masschange.finishRow();
}


void GRMHD::IntegralsWriter::mapQuantities(
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

	// volume form for integration
	const double NumberOfLagrangePointsPerAxis = GRMHD::AbstractGRMHDSolver::PatchSize + 1;
	double scaling = tarch::la::volume(sizeOfPatch* (1.0/NumberOfLagrangePointsPerAxis));
	statistics.addValue(scaling, 1);

	// Mass Accretion Rate

	// We start to compute the accretion rate at the r_excision
	double rmin = 2.1;
	// And we stop integration at some specific detector (r_max)
	double rmax = 2.2;
	double r = norm2(x);

	if(r > rmin && r < rmax) {
	  double mdot;
	  double vx;
      	  double vy;
	  double vz;
	}


	
	// reduce the conserved quantities
	conserved.addValue(Q, scaling);

	// reduce the primitive quantities
	double V[nVar];
	int err;
	pdecons2prim_(V, Q, &err);
	primitives.addValue(V, scaling);

	// now do the convergence test, as we have exact initial data
	double ExactCons[nVar];
	double ExactPrim[nVar];
	const double *xpos = x.data();
	
	initialdata_(xpos, &timeStamp, ExactCons);
	pdecons2prim_(ExactPrim, ExactCons, &err);
	
	double localError[nVar];
	for(int i=0; i<nVar; i++) {
		localError[i] = std::abs(V[i] - ExactPrim[i]);
	}
	
	errors.addValue(localError, scaling);
}


