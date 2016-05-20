#ifndef __SRHD_APPLICATION_H__
#define __SRHD_APPLICATION_H__

const int nVar = 5; // this should be given by ExaHype.

// Ideal EOS properties
const double eos_gamma = 5./3.;

// math:
#include <math.h>
// C polyfills. Feel free to replace with C++98 <algorithms.h> 
inline double max(double a, double b) { return a>b ? a : b; }
inline double abs(double a) { return a<0 ? (-a) : a; }


// con2prims.cpp:

bool cons2prim(double* V, const double *Q);
bool prim2cons(const double* V, double *Q);

double rtsafe_c2p_rmhd1(double& X1, double& X2, double& XACC, double& gam, double& d, double& e, double& s2, double& b2, double& sb2, double& w, bool& failed);
void func_c2p_rmhd1(double& x, double& f, double& df, double& gam, double& d, double& e, double& s2, double& b2, double& sb2, double& w);


/************************** Syntactic sugar ********************************/

// labels for indices into Q and V.
#define DEFINE_LABELS \
	const int rho = 0; \
	const int vx = 1; \
	const int vy = 2; \
	const int vz = 3; \
	const int p = 4
/* Usage:
     function foo(double* Q, double* V) {
        DEFINE_LABELS;
        Q[rho] = V[p];
     }
*/


	
// direct references deep into X  (might be V or Q).	
#define SHORTHANDS(X)   double &rho=X[0], &vx=X[1], &vy=X[2], &vz=X[3], &p=X[4]
/* Usage:
   function foo(double* Q) {
      SHORTHANDS(Q);
      rho = p*p + vx;
   }
   Note: Of course this will probably produce -Wunused-but-set-variable warnings.
*/
	
	
const bool FAILED = true;
const bool SUCCESS = false;

#endif /*__SRHD_APPLICATION_H__*/