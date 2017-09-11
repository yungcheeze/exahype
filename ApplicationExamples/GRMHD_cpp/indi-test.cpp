#include <stdio.h>
#include "PDE/PDE.h"
#include "InitialData.h"
#include <cmath>
#include <algorithm>
using namespace std;

#include "Fortran/PDE.h"

// Fortran
extern "C" {
void alfenwave_(const double* x,const double* const t,  double* Q);
//void pdeprim2cons_(double* /* OUT */ Q, double* /* IN */ V);
//void pdecons2prim_(double* /* OUT */ V, double* /* IN */ Q, int* iErr);
}

int vacuumtest() {
	double x[3] = {0};
	double t = 0;
	constexpr int nCheck = 19;
	constexpr int nSize = 30;
	constexpr double eps = 1e-7;
	double QC[nSize];

	VacuumInitialData Q(QC);
	for(int j=0; j<nCheck; j++) {
		printf("QC[%d] = %e\n", j, QC[j]);
	}
}

int main() {
	double x[3] = {0};
	double t = 0;
	double dx = 0.1;
	int N = 1;
	
	constexpr int nCheck = 19;
	constexpr int nSize = 30;
	constexpr double eps = 1e-7;
	double QC[nSize], QF[nSize], VC[nSize], VF[nSize];
	double VVC[nSize], VVF[nSize]; // Cons2Prim test
	fill_n(QC, nSize, 137);
	fill_n(QF, nSize, 137);
	fill_n(VC, nSize, 137);
	fill_n(VF, nSize, 137);
	InitialState q(QC); // NaN default
	
	for(int i=0; i<N; i++) {
		x[1] += dx;
		AlfenWave id(x,t,VC);
		alfenwave_(x,&t,VF);
		
		// Prim check: works:
		printf("x=[%f %f %f]\n", x[0], x[1], x[2]);
		for(int j=0; j<nCheck; j++) {
			double difference = std::abs(VC[j] - VF[j]);
			bool differ = ( difference > eps );
			printf("V[%2d]: VC=%15e VF=%15e %s\t(diff=%e)\n", j, VC[j], VF[j], differ ? "differ" : "SAME", difference);
		}
		
		// understand whether raising works correctly
		//Hydro::Conserved::ConstShadowExtendable
		//public Magneto::ConstShadowExtendable, 
		/*
		ADMBase::Full adm(VC);
		Hydro::Primitives::ShadowExtendable prim(VC);
		Lo<vec::stored3> myvel;
		myvel.lo = 12345;
		*/
		
		//myvel.lo=0; DFOR(i) CONTRACT(j) myvel.lo(i) += prim.vel.up(j) * adm.gam.lo(i,j);
		//adm.gam.lower_vec(prim.vel, myvel);
		/*
		for(int i=0;i<6;i++) printf("gmunu.lo(%d)=%e\n",i,adm.gam.lo.data[i]);
		DFOR(i) DFOR(j) printf("g.lo(%d,%d)=%e\n", i,j, adm.gam.lo(i,j));
		for(int i=0;i<3; i++) printf(" prim.vel.up(%d)=%e\n",i,prim.vel.up(i));
		for(int i=0;i<3;i++) printf("myvel.vel.lo(%d)=%e\n",i,myvel.lo(i));
		
		return -1;
		*/
		
		// Do the c2p:
		pdeprim2cons_(QF, VF);
		GRMHD::Prim2Cons(QC, VC).copyFullStateVector();
		
		// do again the comparison
		printf("Prim2Cons:\n");
		for(int j=0; j<nCheck; j++) {
			double difference = std::abs(QC[j] - QF[j]);
			bool differ = ( difference > eps );
			printf("Q[%2d]: QC=%15e QF=%15e %s\t(diff=%e)\n", j, QC[j], QF[j], differ ? "differ" : "SAME", difference);
		}
		
		// do the conversion back
		int iErr;
		pdecons2prim_(VVF, QF, &iErr);
		GRMHD::Cons2Prim(VVC, QC).copyFullStateVector();
		
		// Compare again:
		printf("Cons2Prim:\n");
		for(int j=0; j<nCheck; j++) {
			double difference = std::abs(VVC[j] - VVF[j]);
			bool differ = ( difference > eps );
			printf("VV[%2d]: VVC=%15e VVF=%15e %s\t(diff=%e)\n", j, VVC[j], VVF[j], differ ? "differ" : "SAME", difference);
		}
		
		// Fluxes in one test direction
		constexpr int dir = 0;
		double FC[3][nSize], FF[3][nSize];
		double *FCp[3] = {FC[0],FC[1],FC[2]};
		// to detect non set values
		std::fill_n(FC[0], 3*nSize, 123456789);
		std::fill_n(FF[0], 3*nSize, 123456789);
		
		PDE(QC).flux(FCp);
		pdeflux_(FF[0], FF[1], FF[2], QF);
		
		printf("Fluxes in direction dir=%i\n", dir);
		for(int j=0; j<nCheck; j++) {
			double difference = std::abs(FC[dir][j] - FF[dir][j]);
			bool differ = ( difference > eps );
			printf("F[%d][%2d]: CF=%15e FF=%15e %s\t(diff=%e)\n", dir, j, FC[dir][j], FF[dir][j], differ ? "differ" : "SAME", difference);
		}
	}
}
