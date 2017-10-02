/**
 * Initial data for testing the new C++ GRMHD PDE.
 **/

#include "InitialData.h"

#include <cmath>
constexpr double pi = M_PI;

using namespace std;


AlfenWave::AlfenWave(const double* const x, const double t, double* Q) : VacuumInitialData(Q) {
		// Computes the AlfenWave conserved variables (Q) at a given time t.
	// Use it ie. with t=0 for initial data
	// Use it for any other time ie. for comparison
	
	/*
	! GRID FOR ALFENWAVE:
	!     dimension const                = 2
	!     width                          = 1.0, 0.3
	!     offset                         = 0.0, 0.0
	!     end-time                       = 2.1
	!
	!  maximum-mesh-size              = 0.04
	*/

	constexpr double time_offset = 1.0;
	constexpr double gamma = GRMHD::Parameters::gamma;

	double eta  = 1.;
	double B0   = 1.0;
	double hh = 1.0 + gamma / ( gamma - 1.0) * p0 / rho0;
	double tempaa = rho0 * hh + B0*B0 * ( 1.0 + eta*eta);
	double tempab = 2.0 * eta * B0*B0 / tempaa;
	double tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab*tempab));
	double va2 = B0*B0 / ( tempaa * tempac);
	double vax = sqrt(va2);

	// c2p-invariant: Magnetic field
	DFOR(i) Bmag.up(i) = eta * B0 * sin(2*pi*( x[0] - vax*(t-time_offset)));
	Bmag.up(0) = B0;

	// primitive variables
	DFOR(i) vel.up(i) = - vax * Bmag.up(i) / B0;
	vel.up(0) = 0.0;
	
	//if(DIMENSIONS == 2) { Bmag.up(2) = vel.up(2) = 0; }

	//Prim2Cons(Q,V);
	Dens = rho;
	tau = press;
	DFOR(i) Si.lo(i) = vel.up(i);
} // AlfenWave


