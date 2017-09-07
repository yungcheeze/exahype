
#include "PDE.h"
#include <cmath>

using namespace GRMHD;

void GRMHD::Prim2Cons::perform() {
	/**
	* W: Lorentz factor
	* enth: Enthalpy (h)
	**/

	// The hydro exact known prim2cons
	Dens = rho * W;
	DFOR(i) Si.lo(i) = rho * enth * W*W * vel.up(i); // danger. Looks wrong.
	tau = rho*W * (enth*W - 1) - press;
}

void GRMHD::Prim2Cons::prepare() {
	// We start with vel.up computed and assume vel.lo ALSO computed.
	
	// Determine v^2 = v_i * v^i from  an existing vel.lo.
	double VelVel = 0; DFOR(i) VelVel += vel.lo(i) * vel.up(i);
	
	W = 1. / std::sqrt(1.0 - VelVel);
	enth = 5; // Enthalpy (h) TODO
	// I think that enth = ww= w * lf^2 = (rho+gamma/(gamma-1)*p)*lf^2
}

void GRMHD::Prim2Cons::copyFullStateVector() {
	copy_c2p_invariant(Q + Conserved::Indices::c2p_invariant_start);
	copy_admvars(Q + ADMBase::Indices::adm_start);
}
