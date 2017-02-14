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

GRMHD::IntegralsWriter::IntegralsWriter(GRMHDSolver&  solver) {
	// open all the reductions
	assert( 9 == numReduced );
	
	conserved[0] = new TimeSeriesReductions("output/dens.asc");
	conserved[1] = new TimeSeriesReductions("output/sconx.asc");
	conserved[2] = new TimeSeriesReductions("output/scony.asc");
	conserved[3] = new TimeSeriesReductions("output/sconz.asc");
	conserved[4] = new TimeSeriesReductions("output/ener.asc");
	conserved[5] = new TimeSeriesReductions("output/bx.asc");
	conserved[6] = new TimeSeriesReductions("output/by.asc");
	conserved[7] = new TimeSeriesReductions("output/bz.asc");
	conserved[8] = new TimeSeriesReductions("output/psi.asc");
		
	primitives[0] = new TimeSeriesReductions("output/rho.asc"); // V[0]=Q[0]
	primitives[1] = new TimeSeriesReductions("output/velx.asc");
	primitives[2] = new TimeSeriesReductions("output/vely.asc");
	primitives[3] = new TimeSeriesReductions("output/velz.asc");
	primitives[4] = new TimeSeriesReductions("output/press.asc");
	primitives[5] = new TimeSeriesReductions("output/prim-bx.asc");
	primitives[6] = new TimeSeriesReductions("output/prim-by.asc");
	primitives[7] = new TimeSeriesReductions("output/prim-bz.asc");
	primitives[8] = new TimeSeriesReductions("output/prim-psi.asc");
	
	errors[0] = new TimeSeriesReductions("output/error-rho.asc");
	errors[1] = new TimeSeriesReductions("output/error-velx.asc");
	errors[2] = new TimeSeriesReductions("output/error-vely.asc");
	errors[3] = new TimeSeriesReductions("output/error-velz.asc");
	errors[4] = new TimeSeriesReductions("output/error-press.asc");
	errors[5] = new TimeSeriesReductions("output/error-bx.asc");
	errors[6] = new TimeSeriesReductions("output/error-by.asc");
	errors[7] = new TimeSeriesReductions("output/error-bz.asc");
	errors[8] = new TimeSeriesReductions("output/error-psi.asc");
	
	statistics = new TimeSeriesReductions("output/volform.asc");
}


GRMHD::IntegralsWriter::~IntegralsWriter() {
	// delete all the TimeSeriesReductions. Not that important
	// as this is a kind of a Singleton object.
}


void GRMHD::IntegralsWriter::startPlotting(double time) {
	for(int i=0; i<numReduced; i++) {
		conserved[i]->initRow(time);
		primitives[i]->initRow(time);
		errors[i]->initRow(time);
	}
	statistics->initRow(time);
}


void GRMHD::IntegralsWriter::finishPlotting() {
	for(int i=0; i<numReduced; i++) {
		conserved[i]->writeRow();
		primitives[i]->writeRow();
		errors[i]->writeRow();
	}
	statistics->writeRow();
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

	const double NumberOfLagrangePointsPerAxis = GRMHD::AbstractGRMHDSolver::Order + 1;
	//const double NumberOfUnknownsPerGridPoint = nVar;
	
	// Gauss-Legendre weights from pos argument
	double wx = kernels::gaussLegendreWeights[GRMHD::AbstractGRMHDSolver::Order][pos[0]];
	double wy = kernels::gaussLegendreWeights[GRMHD::AbstractGRMHDSolver::Order][pos[1]];
	double wz = 1;
	#ifdef Dim3
	wz = kernels::gaussLegendreWeights[GRMHD::AbstractGRMHDSolver::Order][pos[2]];
	#endif

	// volume form for integration
	double scaling = tarch::la::volume(sizeOfPatch);
	statistics->addValue(scaling, 1);

	// reduce the conserved quantities
	for (int i=0; i<numReduced; i++)
		conserved[i]->addValue( Q[i], scaling );

	// reduce the primitive quantities
	double V[nVar];
	int err;
	pdecons2prim_(V, Q, &err);
	for(int i=0; i<numReduced; i++)
		primitives[i]->addValue( V[i], scaling );

	// now do the convergence test, as we have exact initial data
	double ExactCons[nVar];
	double ExactPrim[nVar];
	const double *xpos = x.data();
	
	initialdata_(xpos, ExactCons);
	pdecons2prim_(ExactPrim, ExactCons, &err);
	
	double localError[nVar];
	for(int i=0; i<numReduced; i++) {
		localError[i] = std::fabs(V[i] - ExactPrim[i]);
		//localError[i] = abs(Q[i] - ExactCons[i]);
		errors[i]->addValue( localError[i], scaling*wx*wy*wz );
		
		/*if(i==0) {
			printf("RHO; t=%e x=(%e,%e) num=%14e exact=%14e error=%e\n", timeStamp, x[0],x[1], V[0], ExactPrim[0], localError[0]);
		}*/
    // Uncomment for debugging reasons
//    std::cout << "V_ana["<<i<<"]="<<ExactPrim[i]<<",V_h["<<i<<"]="<<V[i]<< std::endl;
//    std::cout << "Q_ana["<<i<<"]="<<ExactCons[i]<<",Q_h["<<i<<"]="<<Q[i] <<std::endl;
	}
}


