/**
 * A C++ variant of the GRMHD equations for ExaHyPE.
 * Read PDE.h for further documentation.
 * Written at 2017-08-06 by SvenK.
 **/

#include "PDE.h"
#include <cmath> // NAN

/******************************************************************************
 * Status of the PDE system:
 *   It is now equal to the schemes of BHAC and ECHO, i.e. it seems to be valid
 *   for the cowling approximation but not neccessarily for a dynamical
 *   spacetime.
 * 
 * Important TODO: Stripping of the determinant of the metric is missing.
 *****************************************************************************/

#define S(x) printf(#x " = %e\n", x);
#define SI(x) S(x(0));S(x(1));S(x(2));

void GRMHD::FluxBase::prepare() {
	// Sij is the 3-Energy-Momentum tensor: We only need S^i_j in the flux.
	SYMFOR(i,j) Sij.ul(i,j) = Si.up(i)*vel.lo(j) + ptot*delta(i,j) - Bmag.up(i)*Bmag.lo(j)/WW - BmagVel * vel.up(i) * Bmag.lo(j);
	
	// Zeta is the transport velocity (curly V in BHAC paper)
	DFOR(k) zeta.up(k) = alpha*vel.up(k) - beta.up(k);
}

void GRMHD::FluxBase::compute(Flux& flux, int k) {
	// F^k: flux in direction k.

	flux.Dens = Dens * zeta.up(k);
	DFOR(i) flux.Si.lo(i) = alpha*Sij.ul(k,i) - beta.up(k)*Si.lo(i);
	flux.tau = alpha*(Si.up(k) - vel.up(k)*Dens) - beta.up(k)*tau;
	DFOR(j) flux.Bmag.up(j) = zeta.up(k)*Bmag.up(j) - zeta.up(j)*Bmag.up(k);
	
	// Constraint damping contributions:
	DFOR(j) flux.Bmag.up(j) -= Bmag.up(k)*beta.up(j);
	flux.phi = alpha*Bmag.up(k) - phi*beta.up(k);
}

void GRMHD::Fluxes::zeroMaterialFluxes() {
	DFOR(k) {
		F[k].alpha = 0;
		F[k].beta.up = 0;
		F[k].gam.lo = 0;
	}
}

GRMHD::PDE& GRMHD::PDE::RightHandSide(/* const double* const Q, */const double* const gradQ_Data, double* Source_data) {
	GRMHDSystem::Shadow source(Source_data);
	const Gradients grad(gradQ_Data+0,gradQ_Data+nVar,gradQ_Data+2*nVar);
	return RightHandSide(grad,source);
}
			
GRMHD::PDE& GRMHD::PDE::RightHandSide(const double* const Qx, const double* const Qy, const double* const Qz, double* fusedSource_Data) {
	GRMHDSystem::Shadow source(fusedSource_Data);
	const Gradients grad(Qx,Qy,Qz);
	return RightHandSide(grad,source);
}
	
GRMHD::PDE& GRMHD::PDE::RightHandSide(const Gradients& grad, Source& source) {
	// This is the fusedSource = -NCP + AlgebraicSource
	
	// Sij is the 3-Energy-Momentum tensor. We need S^{ij} and S^i_j in the NCP.	
	UpSym<sym::stored<3>, sym::stored<3>> Sij;
	SYMFOR(i,j) Sij.up(i,j) = Si.up(i)*vel.up(j) + gam.up(i,j)*ptot - Bmag.up(i)*Bmag.up(j)/WW - BmagVel * vel.up(i) * Bmag.up(j);
	Sij.ul=0; SYMFOR(i,j) CONTRACT(k) Sij.ul(i,j) += Sij.up(i,k) * gam.lo(j,k); // Sij^i_j = Sij^{ik} gam_{jk}

	// source for D
	source.Dens = 0;
	
	// U is the total energy seen by Eulerian observer
	// D is the density
	// tau is the rescaled energy density.
	// We evolve tau for an improved low density region.
	double U = tau + Dens;

	// S_i
	DFOR(i) {
		source.Si.lo(i) = - U * grad.lo(i).alpha;
		CONTRACT(k) source.Si.lo(i) += Si.lo(k) * grad.lo(i).beta.up(k);
		CONTRACT2(l,m) source.Si.lo(i) += alpha/2. * Sij.up(l,m) * grad.lo(i).gam.lo(l,m);
	}
	
	// tau
	CONTRACT(k) source.tau = -Si.up(k) * grad.lo(k).alpha;
	// further BHAC terms for tau:
	CONTRACT2(i,j) source.tau += Sij.ul(j,i) * grad.lo(j).beta.up(i);
	CONTRACT3(i,k,j) source.tau += 1./2. * Sij.up(i,k)*beta.up(j) * grad.lo(j).gam.lo(i,k);
	
	//CONTRACT3(k,i,j) { printf("%d,%d,%d: ",k,i,j); S(grad.lo(k).gam.lo(i,j)); }
	//S(source.tau);
	
	// Starting from here:
	// Divergence Cleaning/Constraint damping sources for B^j and Phi.
	// In case no cleaning would be applied, source.Bmag == source.phi == 0.
	
	// B^i magnetic field
	DFOR(i) {
		source.Bmag.up(i) = 0;
		CONTRACT(k) source.Bmag.up(i) -= Bmag.up(k) * grad.lo(k).beta.up(i);
		CONTRACT(j) source.Bmag.up(i) -= alpha * gam.up(i,j) * grad.lo(j).phi;
	}
	
	// phi
	source.phi = - alpha * damping_term_kappa * phi; // algebraic source
	CONTRACT(k) source.phi += Bmag.up(k) * grad.lo(k).alpha;
	CONTRACT(k) source.phi -= phi * grad.lo(k).beta.up(k);
	CONTRACT3(k,l,m) source.phi += alpha/2 * gam.up(l,m) * beta.up(k) * grad.lo(k).gam.lo(l,m);
	
	return *this; // chainable
}

/*
void GRMHD::PDE::algebraicSource(double* Source_data) {
	ConservedVariables Source(Source_data);
	// we could put the following term here:
	// Source.phi = - alpha * kappa * phi; // algebraic source
	// All other terms should be zero.
}

void GRMHD::PDE::fusedSource(const double* const gradQ_Data, double* Source) {
	// Variables Source(Source_data, nVar);
	// 1. Compute NCP
	nonConservativeProduct(gradQ_Data, Source);
	// 2. Flip sign
	NVARS(m) Source[m] = - Source[m];
	// 3. Add the algebraic source terms:
	// CONTRACT2(l,m) Source.tau += lapse * Sij.up(l,m) * Kextr.lo(l,m);
	// Source.phi = - lapse * kappa * phi;
	// 4. Done.
}
*/

void GRMHD::PDE::eigenvalues(const double* const Q, const int d, double* lambda) {
	//std::fill_n(lambda, nVar, 1.0);
	NVARS(m) lambda[m] = 1.0;
}
