/**
 * A C++ variant of the GRMHD equations for ExaHyPE.
 * 
 * This formulation uses the Tensish library for the Tensor abstraction.
 * We mimic the ExaHyPE notation of state vector abstractions with structures.
 * We copy all neccessary parts for practical purposes (scalar, vector, tensor
 * notation) as well as distinguation between conserved, primitive, material
 * parameters.
 *
 * Written at 2017-08-06 by SvenK.
 **/
#ifndef GRMHD_PDE_CPP
#define GRMHD_PDE_CPP

#include "tensish.cpph"

constexpr int nVar = 10; //GRMHD::Whatever::NumberOfVariables;
#define NVARS(i) for(int i=0; i<nVar; i++)

namespace GRMHD {
	class State;
	using sym::delta; // Kronecker symbol
	
	/**
	 * These structs store the index positions of fields, i.e. Q.Dens=0.
	 * They are similiar to the VariableShortcuts which are generated
	 * by the ExaHyPE java toolkit. However, by defining our own we are
	 * independent and can make use of the functional notation.
	 **/
	struct ConservedIndices { // aka VariableShortcuts
		// Conservatives
		static constexpr int   Dens = 0;
		static constexpr int   Si_lo(int i) { return 1+i; }
		static constexpr int   tau = 4;
		// Conserved (C2P invariant)
		static constexpr int   Bmag_up(int i) { return 5+i; }
		static constexpr int   Phi = 8;
		// material parameters
		static constexpr int   lapse = 9;
		static struct { static constexpr int lo(int i) { return 10+i; } } shift; // TODO: really lo?
		static struct { static constexpr int lo(int i, int j) { return 13+sym::index(i,j); } } gam; // TODO: Really lo? 
	};
	
	struct PrimitiveIndices { // no correspondence in ExaHyPE
		static constexpr int rho = 0;
		// TODO: check if up
		static constexpr int vec_up = 1;
		static constexpr int press = 4;
	};
	
	struct PDE; // forward
	struct Primitives;
	struct ADMBase;
	template<typename vector, typename vector_lo, typename vector_up, typename scalar> struct Conserved;
	
	typedef GRMHD::Conserved<
		const double* const,
		ConstLo<vec::stored<3>,vec::const_shadow>,
		ConstUp<vec::const_shadow, vec::stored<3>>,
		const double> ConservedConstants;
	
	typedef GRMHD::Conserved<
		double* const,
		StoredLo<vec::stored<3>, vec::shadow>,
		StoredUp<vec::shadow, vec::stored<3>>,
		double> ConservedVariables;
} // namespace GRMHD

// Variables store the actual (conserved) variables in a shorter/better readable/friendlier
// manner than the exahype generated variables header file
template<typename vector, typename vector_lo, typename vector_up, typename scalar>
struct GRMHD::Conserved {
	vector Q; // const double* const
	const GRMHD::ConservedIndices Qi;
	
	// Conserved Scalars
	scalar &Dens, &tau, &phi; // const double
	
	// Conserved vectors
	vector_lo Si;
	vector_up Bmag;
	
	Conserved(vector Q_) :
		Q(Q_),
		Qi(),
		// Conserved scalars
		Dens (Q[Qi.Dens]),
		tau  (Q[Qi.tau]),
		phi  (Q[Qi.Phi]),
		// Conserved vectors
		Si   (Q+Qi.Si_lo(0)),
		Bmag (Q+Qi.Bmag_up(0))
		{}
	
	void toPrimitives(double* V); // C2P Cons2Prim operation
};

struct GRMHD::Primitives { // are always read only in our formulation, no need for P2C after ID
	double V[nVar];
	const GRMHD::PrimitiveIndices Vi;
	
	// Primitive Scalars
	const double &rho, &press;
	
	// Primitive vectors
	ConstUp<vec::const_shadow, vec::stored<3>> vel;
	
	Primitives(const double* const Q) :
		Vi(),
		// Primitive Scalars
		rho  (V[Vi.rho]),
		press(V[Vi.press]),
		// Primitive vectors
		vel  (V+Vi.vec_up)
		{
			//Cons2Prim(Q, V); // TODO include
		}
	
	void toConserved(double* Q); // P2C Prims2Cons operation
};

/*
// You could define this if you'd like:
struct GRMHD::HydroBase : public GRMHD::ConservedVariables, public GRMHD::PrimitiveVariables {
	HydroBase(const double* const Q) :
		ConservedVariables(Q),
		PrimitiveVariables(Q) {}
};
*/

class GRMHD::ADMBase { // Material Parameters, always read only
	const GRMHD::ConservedIndices Ai;
public:
	const double &alpha;
	// (Conserved) Material parameter vector
	ConstLo<vec::stored<3>, vec::const_shadow> beta;
	
	// (Conserved) Material parameter tensor
	metric3 gam;
	
	ADMBase(const double* const Q) :
		Ai(),
		alpha(Q[Ai.lapse]),
		beta (Q+Ai.shift.lo(0)),
		gam  (Q+Ai.gam.lo(0,0))
		{}
	
	void complete(/* some vector */) {
		// should provide this.
	}
};


struct GRMHD::PDE : public GRMHD::ConservedConstants, public GRMHD::Primitives, public GRMHD::ADMBase {
	// 3-Energy momentum tensor
	// Full<sym::stored<3>, sym::stored<3>> Sij;
	
	PDE(const double*const Q_) :
		ConservedConstants(Q_),
		Primitives(Q_),
		ADMBase(Q_)
		{ prepare(); }

	// Quantities for computing the energy momentum tensor
	double WW, BmagBmag, BmagVel, ptot;
	void prepare(); // prepare that shit
	
	void flux(/* const double* const Q, */ double** Fluxes);
	void nonConservativeProduct(/* const double* const Q, */ const double* const gradQ_Data, double* BgradQ_Data);
	void algebraicSource(/* const double* const Q, */ double* Source_data);
	void fusedSource(/* const double* const Q, */const double* const gradQ_Data, double* Source_data);
	
	// currently trivial:
	static void eigenvalues(const double* const Q, const int d,double* lambda);
};


#endif /* GRMHD_PDE_CPP */
