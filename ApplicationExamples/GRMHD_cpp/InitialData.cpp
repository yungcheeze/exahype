/**
 * Initial data for testing the new C++ GRMHD PDE.
 **/

#include "PDE/PDE.h"

#include <cmath>
constexpr double pi = M_PI;

using namespace std;
using namespace GRMHD;

struct InitialState : public Conserved::Shadow, public ADMBase::Shadow, public Primitives::Stored {
	InitialState(double* Q) : Conserved::Shadow(Q), ADMBase::Shadow(Q), Primitives::Stored() {}
};

/// Defaults to vacuum
struct VacuumInitialData : public InitialState {
	static constexpr double rho0 = 1.;
	static constexpr double p0   = 1.;
	
	void flatSpace() {
		// ADM base: Flat space
		alpha = 1.0;
		DFOR(i) beta.up(i) = 0;
		SYMFOR(i,j) gam.lo(i,j) = sym::delta(i,j);
	}
	
	VacuumInitialData(double* Q) : InitialState(Q) {
		// c2p invariants:
		DFOR(i) Bmag.up(i) = 0;
		phi = 0; // damping term
		
		// primitives:
		rho = rho0;
		press = p0;
		DFOR(i) vel.up(i) = 0;
		
		// ADM
		flatSpace();
		
		// Don't forget to call the Cons2Prim afterwards.
	}
};

struct AlfenWave : public VacuumInitialData {
	AlfenWave(const double* const x, const double t, double* Q) : VacuumInitialData(Q) {
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
		constexpr double gamma = 1.0; // TODO

		double eta  = 1.;
		double B0   = 1.0;
		double hh = 1.0 + gamma / ( gamma - 1.0) * p0 / rho0;
		double tempaa = rho0 * hh + B0*B0 * ( 1.0 + eta*eta);
		double tempab = 2.0 * eta * B0*B0 / tempaa;
		double tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab*tempab));
		double va2 = B0*B0 / ( tempaa * tempac);
		double vax = sqrt(va2);

		// c2p-invariant: Magnetic field
		Bmag.up(1) = B0;
		Bmag.up(2) = eta * B0 * cos(2*pi*( x[1] - vax*(t-time_offset)));
		Bmag.up(3) = eta * B0 * sin(2*pi*( x[1] - vax*(t-time_offset)));

		// primitive variables
		vel.up(1) = 0.0;
		vel.up(2) = - vax * Bmag.up(2) / B0;
		vel.up(3) = - vax * Bmag.up(3) / B0;

		Prim2Cons(Q,V);
	}
}; // AlfenWave


