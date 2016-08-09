/**
 * Initial Data for the EulerFlow.
 * 
 * As we don't have parameters for that in the moment, we use environment
 * variables to distuingish at run time which initial data to load.
 * 
 **/

#include "InitialData.h"
#include "Primitives.h"
#include "GeneratedConstants.h"
#include "PastaMatrix.h"

#include <stdlib.h>
#include <math.h>
#include <cstring>


#include "tarch/logging/Log.h"


using namespace std;

/**
 * Source: MD ADERDG F90 CODE
 * ShuVortex produces primitive Variables which have to pass Prim2Con afterwards
 **/
void ShuVortex2D(const double* const  x, double* V, double t=0.0) {
	static const double epsilon = 5.0;
	static const double pi = acos(-1.0);

	double r = sqrt(SQ(x[0]-t-5.)+SQ(x[1]-t-5.));
	double du = epsilon/2./pi*exp(0.5*(1.-r*r))*(5. - x[1] + t);
	double dv = epsilon/2./pi*exp(0.5*(1.-r*r))*(x[0]  - 5.- t);
	double dTemp = -(eos_gamma-1.)*SQ(epsilon)/8./eos_gamma/SQ(pi)*exp(1.-r*r);
	double drho = pow(1.+dTemp, 1./(eos_gamma-1.))-1.;
	double dp   = pow(1.+dTemp, eos_gamma/(eos_gamma-1.))-1.;

	V[0] = 1. + drho;
	V[1] = 1. + du;
	V[2] = 1. + dv;
	V[3] = 0.0;
	V[4] = 1. + dp;
}

/**
 * MovingGauss2D is a moving gaussian matter distribution where it is simple
 * to give an analytic result.
 * 
 * 
 **/
void MovingGauss2D(const double* const x, double* V, double t=0.0) {
	Pasta::vec2 xvec(x);
	Pasta::vec2 v0({ 0.5, 0 });
	//Pasta::vec2 v0({ 0.0, 0 });
	Pasta::vec2 x0({ 0.5, 0.5 });
	double width = 0.25;
	
	V[0] = 0.5 + 0.2 * exp(- (xvec - x0 - v0*t).norm() / pow(width, MY_DIMENSIONS)); // rho
	V[1] = v0(0);
	V[2] = v0(1);
	V[3] = 0.;
	V[4] = 1.; // pressure
}


/**
 * DiffusingGauss is not a consistent solution of Euler's equations, but instead
 * a perturbation which immediately changes rho, vx, vy and vz. It corresponds
 * to consistent initial data with roughly
 *   rho = exp( - (r - v0*t)**2)
 *   vx  = v0*t * exp(-y**2)  // actually even more complicated, it is more
 *   vy  = v0*t * exp(-x**2)  // a velocity on a kind of ring of radius v0*t
 *   E   = p/(gamma-1) + rho/2 * v**2 = p/(gamma-1) + alpha*exp(- (r-v0*t)**2)
 *   p   = 1
 * However, it is much more complicated to write the closed form solution instead
 * of the perturbation approach. However, the closed form solution allows to
 * specify the solution at any time while while the initial perturbation form
 * does *not* allow to specify the solution.
 * 
 * Attention: If the prim2con in the Primitives.h is used, for some reason it
 * surpresses this solution and rho=1, p=1, E=p/(gamma-1) is the resolut on the
 * whole grid. I am not sure why this happens, it seems to be an error.
 **/
void DiffusingGauss(const double* const  x, double* Q) {
	Pasta::vec2 xvec(x);
	Pasta::vec2 x0({ 0.5, 0.5 });
	double width = 0.05;

	Q[0] = 1.;
	Q[1] = 0.;
	Q[2] = 0.;
	Q[3] = 0.;
	Q[4] = 1./(eos_gamma - 1) + exp(- (xvec - x0).norm() / pow(width, MY_DIMENSIONS) ) * 2;
}


void InitialData(const double* const  x, double* Q, double t) {
        static tarch::logging::Log _log( "" );
        static bool wroteAboutInitialData(false);

	const char* default_id = "DiffusingGauss";
	const char* id = getenv("EXAHYPE_INITIALDATA");
	if(!id) {
          if(!wroteAboutInitialData) logInfo( "InitialData(double*,double,double)", "Using default ID");
	  id = default_id;
	}

	//logInitialData("Have read '%s'\n", id);
	std::string sid(id);
	if(sid == "ShuVortex") {
                if(!wroteAboutInitialData) logInfo( "InitialData(double*,double,double)", "Loading ShuVortex Initial Data");
		// ShuVortex gives us primitive data
                double V[MY_NUMBER_OF_VARIABLES];
		ShuVortex2D(x, V, t);
                prim2con(Q, V);
	} else if(sid == "MovingGauss2D") {
                if(!wroteAboutInitialData) logInfo( "InitialData(double*,double,double)", "Loading moving Gauss");
		double V[MY_NUMBER_OF_VARIABLES];
		MovingGauss2D(x, V, t);
		prim2con(Q, V);
	} else if(sid == "DiffusingGauss") {
                if(!wroteAboutInitialData) logInfo( "InitialData(double*,double,double)", "Loading diffusing Gauss Initial Data");
		// default:
		DiffusingGauss(x, Q);
	} else {
                logError( "InitialData(double*,double,double)", "Do not understand requested Initial Data key");
		exit(-42);
	}
	wroteAboutInitialData = true;
}
