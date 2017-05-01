
#include "AbstractGRMHDSolver_ADERDG.h"
#include "GRMHDSolver_ADERDG_Variables.h"
#include "InitialData/InitialData.h"
#include "Fortran/PDE.h"
#include <stdlib.h>
#include <stdio.h>

#ifndef PIZZATOV_AVAILABLE

pizzatov::pizzatov() {
	printf("Cannot call Pizza as not compiled with -DPIZZA_TOV_AVAILABLE");
	abort();
}
void pizzatov::Interpolate(const double* x, double t, double* Q) {}

#else /* PIZZATOV_AVAILABLE */

#include "pizza_tovfront/pointwise_tov.h"

using namespace Pizza;
using namespace Pizza::TOV;

constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;

pizzatov::pizzatov() {
	tov = new pointwise_tov();
	
	// here is the place to set some TOV parameters,
	// if needed. Such as the EOS.
	
	std::string EOS_BU =
		"name      = BU (Polytropic)\n"
		"type      = polytrope\n"
		"poly_n    = 1\n"
		"poly_rmd  = 6.176e+18\n"
	;
	tov->read_eos_from_string(EOS_BU);
	
	// set mapping to GRMHD application
	// rho:1,vel:3,E:1,B:3,psi:1,lapse:1,shift:3,gij:6
	GRMHD::AbstractGRMHDSolver_ADERDG::VariableShortcuts vars;

	// GRMHD Material parameters
	tov->adm.alp = vars.lapse;
	tov->adm.shift1 = vars.shift + 0;
	tov->adm.shift2 = vars.shift + 0;
	tov->adm.shift3 = vars.shift + 0;
	tov->adm.gxx = vars.gij + 0;
	tov->adm.gxy = vars.gij + 1;
	tov->adm.gxz = vars.gij + 2;
	tov->adm.gyy = vars.gij + 3;
	tov->adm.gyz = vars.gij + 4;
	tov->adm.gzz = vars.gij + 5;
	
	// GRMHD actual computation variables
	// Note that PizzaTOV sets primitives and we have to do the cons2prim afterwards.
	tov->hydro.rho = vars.rho;
	tov->hydro.velx = vars.vel + 0;
	tov->hydro.vely = vars.vel + 1;
	tov->hydro.velz = vars.vel + 2;
	tov->hydro.press = vars.E;
	
	// do the main computation (takes no time)
	tov->compute_star();
}

void pizzatov::Interpolate(const double* x, double t, double* Q) {
	double V[nVar];
	tov->initial_data(x, V);
	pdeprim2cons_(Q, V);
}

#endif  /* PIZZATOV_AVAILABLE */
