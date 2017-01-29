#include "ComputeGlobalIntegralsLegendre.h"
#include "Primitives.h"
#include "InitialData.h"

#include "kernels/GaussLegendreQuadrature.h"

Euler::ComputeGlobalIntegralsLegendre::ComputeGlobalIntegralsLegendre(MyEulerSolver&  solver) {
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


Euler::ComputeGlobalIntegralsLegendre::~ComputeGlobalIntegralsLegendre() {
	// delete all reductions.
	// @todo  the deletes are missing here
}


void Euler::ComputeGlobalIntegralsLegendre::startPlotting(double time) {
	this->time = time;
	for(int i=0; i<nVar; i++) {
		conserved[i]->initRow(time);
		primitives[i]->initRow(time);
		errors[i]->initRow(time);
	}
	statistics->initRow(time);
}


void Euler::ComputeGlobalIntegralsLegendre::finishPlotting() {
	for(int i=0; i<nVar; i++) {
		conserved[i]->writeRow();
		primitives[i]->writeRow();
		errors[i]->writeRow();
	}
	statistics->writeRow();
}


void Euler::ComputeGlobalIntegralsLegendre::mapQuantities(
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
  double scaling = tarch::la::volume(sizeOfPatch);
  statistics->addValue(scaling, 1);

  // Gauss-Legendre weights from pos argument
  double wx = kernels::gaussLegendreWeights[Euler::MyEulerSolver::order][pos[0]];
  double wy = kernels::gaussLegendreWeights[Euler::MyEulerSolver::order][pos[1]];
  double wz = 1;
  #ifdef Dim3
  wz = kernels::gaussLegendreWeights[Euler::MyEulerSolver::order][pos[2]];
  #endif

	// reduce the conserved quantities
	for (int i=0; i<nVar; i++)
		conserved[i]->addValue( Q[i], scaling );

	// reduce the primitive quantities
	double V[nVar];
	cons2prim(V, Q);
	for(int i=0; i<nVar; i++)
		primitives[i]->addValue( V[i], scaling );

	// now do the convergence test, as we have exact initial data
	double Exact[nVar];
	// TODO: Need a way to access _data in tarch::la::Vector.
	double xpos[DIMENSIONS];
	for(int i=0; i<DIMENSIONS; i++) xpos[i] = x[i];
	
	idfunc(xpos, Exact, time);
	double ExactPrim[nVar];
	cons2prim(ExactPrim, Exact);

	// Uncomment for debugging reasons
//	std::cout << "x="<<x.toString()<<"J="<<scaling << std::endl;
//  std::cout << "wx="<<wx<<",wy="<<wy<<",wz="<<wz << std::endl;

	double localError[nVar];
	for(int i=0; i<nVar; i++) {
		localError[i] = abs(V[i] - ExactPrim[i]);
		errors[i]->addValue( localError[i], scaling*wx*wy*wz );

		// Uncomment for debugging reasons
//		std::cout << "V_ana["<<i<<"]="<<ExactPrim[i]<<",V_h["<<i<<"]="<<V[i]<< std::endl;
//		std::cout << "Q_ana["<<i<<"]="<<Exact[i]<<",Q_h["<<i<<"]="<<Q[i] <<std::endl;
	}
}


