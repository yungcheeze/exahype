/**
 * A C++ variant of the GRMHD equations for ExaHyPE.
 * Read PDE.h for further documentation.
 * Written at 2017-08-06 by SvenK.
 **/

#include "PDE.h"
#include "really.h"

void GRMHD::PDE::prepare() {
	// NEED to prepare:
	// Generic for Sij and preparation:
	DFOR(i) CONTRACT(j) Bmag.lo(i) = gam.lo(i,j) * Bmag.up(j);
	DFOR(i) CONTRACT(j) vel.lo(i)  = gam.lo(i,j) * vel.up(j);
	// S^i is needed in both flux and ncp
	DFOR(i) CONTRACT(j) Si.up(i)   = gam.up(i,j) * Si.lo(j);

	// TODO Clarify: Do we store beta up or beta lo
	
	WW = SQ(Dens/rho); // W^2
	BmagBmag = 0; CONTRACT(k) BmagBmag += Bmag.lo(k)*Bmag.up(k); // B^2
	BmagVel = 0;  CONTRACT(j) BmagVel  += Bmag.up(j)*vel.lo(j); // B^j * v_j
	ptot = press + 0.5*(BmagBmag/WW + SQ(BmagVel)); // total pressure incl magn. field
}

void GRMHD::PDE::flux(double** Fluxes) {
	// Fluxes Shape: F[nDim][nVar]
	
	// Sij is the 3-Energy-Momentum tensor: We only need S^i_j in the flux.
	Ul<sym::stored<3>> Sij;
	DFOR(i) DFOR(j) Sij.ul(i,j) = Si.up(i)*vel.lo(j) + ptot*delta(i,j) - Bmag.up(i)*Bmag.lo(j)/WW - BmagVel * vel.up(i) * Bmag.lo(j);
	
	// Zeta is the transport velocity (curly V in BHAC paper)
	Up<vec::stored<3>> zeta; DFOR(k) zeta.up(k) = alpha*vel.up(k) - beta.up(k);

	DFOR(k) { // F^k: flux in direction k.
		ConservedVariables Flux(Fluxes[k]);
		Flux.Dens = Dens * zeta.up(k);
		DFOR(i) Flux.Si.lo(i) = alpha*Sij.ul(k,i) - beta.up(k)*Si.lo(i);
		Flux.tau = alpha*(Si.up(k) - vel.up(k)*Dens) - beta.up(k)*tau;
		DFOR(j) Flux.Bmag.up(j) = zeta.up(k)*Bmag.up(j) - zeta.up(j)*Bmag.up(k) - Bmag.up(k)*beta.up(j);
		Flux.phi = alpha*Bmag.up(k) - phi*beta.up(k);
	}
}

void GRMHD::PDE::nonConservativeProduct(const double* const gradQ_Data, double* BgradQ_Data) {
	// Shape: BgradQ(nVar)
	ConservedVariables BgradQ(BgradQ_Data);
	const Really::const_shadow gradQ(gradQ_Data, nDim, nVar);
	
	// Sij is the 3-Energy-Momentum tensor. We need S^{ij} and S^i_j in the NCP.	
	UpSym<sym::stored<3>, sym::stored<3>> Sij;
	DFOR(i) DFOR(j) Sij.up(i,j) = Si.up(i)*vel.up(j) + gam.up(i,j)*ptot - Bmag.up(i)*Bmag.up(j)/WW - BmagVel * vel.up(i) * Bmag.up(j);
	Sij.ul_from_up(gam);

	// D
	BgradQ.Dens = 0;
	
	// U is the total energy seen by Eulerian observer
	// D is the density
	// tau is the rescaled energy density.
	// We evolve tau for an improved low density region.
	double U = tau + Dens;

	// S_i
	DFOR(i) {
		BgradQ.Si.lo(i) = U * gradQ(i, Qi.lapse);
		CONTRACT(k) BgradQ.Si.lo(i) -= Si.lo(k) * gradQ(i, Qi.shift.lo(k));
		CONTRACT2(l,m) BgradQ.Si.lo(i) -= alpha/2. * Sij.up(l,m) * gradQ(i, Qi.gam.lo(l,m));
	}
	
	// tau
	CONTRACT(k) BgradQ.tau = Si.up(k) * gradQ(k, Qi.lapse); // Alejandros term
	// further BHAC terms for tau:
	CONTRACT2(i,j) BgradQ.tau -= Sij.ul(j,i) * gradQ(j, Qi.shift.lo(i)); // check: shift_lo!? shift_up!!!
	CONTRACT3(i,k,j) BgradQ.tau -= alpha/2. * Sij.up(i,k)*beta.up(k) * gradQ(j, Qi.gam.lo(i,k));
	
	// B^i magnetic field
	DFOR(i) {
		BgradQ.Bmag.up(i) = 0;
		CONTRACT(k) BgradQ.Bmag.up(i) += Bmag.up(k) * gradQ(k, Qi.shift.lo(i));
		CONTRACT(j) BgradQ.Bmag.up(i) += alpha * gam.up(i,j) * gradQ(j, Qi.Phi);
	}
	
	// phi
	const double kappa = 5;
	BgradQ.phi = - alpha * kappa * phi; // algebraic source
	CONTRACT(k) BgradQ.phi -= Bmag.up(k) * gradQ(k, Qi.lapse);
	CONTRACT(k) BgradQ.phi += phi * gradQ(k, Qi.shift.lo(k));
	                             //   CHECK: Bmag.up(k) or shift.up(k)???
	CONTRACT3(k,l,m) BgradQ.phi += alpha/2 * Bmag.up(k) * gam.up(l,m) * gradQ(k, Qi.gam.lo(l,m));
}

/*
void GRMHD::PDE::algebraicSource(double* Source_data) {
	ConservedVariables Source(Source_data);
	// we could put the following term here:
	// Source.phi = - alpha * kappa * phi; // algebraic source
}
*/

void GRMHD::PDE::fusedSource(/* const double* const Q, */const double* const gradQ_Data, double* Source) {
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

void GRMHD::PDE::eigenvalues(const double* const Q, const int d, double* lambda) {
	NVARS(m) lambda[m] = 1.0;
}

int main() {
	return 0;
}
