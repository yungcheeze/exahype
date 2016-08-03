#ifndef __EULER_PRIMITIVES__
#define __EULER_PRIMITIVES__

// Ideal equation of state's Gamma
static const double eos_gamma = 5./3.;

// helper: Square of a.
inline double SQ(double a) { return a*a; }

// another pseudo vector helper
inline double SQ3(const double* const v, int start=1) {
	double res=0;
	for(int i=start;i<3; i++) res += SQ(v[i]);
	return res;
}

/** The reason for this being in a header is pure lazyness and triviality of the C2P for Euler */

inline void cons2prim(double* V, const double* Q) {
	// Euler c2p
	double p = (eos_gamma-1)*(Q[4] - 0.5/Q[0] * SQ3(Q) );
	
	V[0] = Q[0]; // fluid density
	V[1] = Q[1] / Q[0]; // fluid velocities
	V[2] = Q[2] / Q[0];
	V[3] = Q[3] / Q[0];
	V[4] = p; // fluid pressure
}

inline void prim2con(double* Q, const double* V) {
	Q[0] = V[1]; // fluid density
	Q[1] = V[1] * V[1]; // momentum
	Q[2] = V[1] * V[2];
	Q[3] = V[1] * V[3];
	// total energy = internal energy + kinetic energy
	Q[4] = V[4] / (eos_gamma-1) + 0.5*V[0] * SQ3(V);
}



#endif /* __EULER_PRIMITIVES__ */