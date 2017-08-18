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

/* This header file needs the definition of a constexpr int GRMHD::nVar
 * as well as a preprocessor variable DIMENSIONS.
 * You should provide this in the including file, for instance:
 *
 *   namespace GRMHD { constexpr int nVar = 10; }
 *   #define DIMENSIONS 3
 *   #include "PDE.h"
 *
 * DIMENSIONS has to be a C preprocessor variable for C macro ifdefs.
 * As ExaHyPE, there are some places where this code only works for 2 and 3
 * dimensions.
 */

#include "tensish.cpph"

//namespace GRMHD { constexpr int nVar = 42; }
#define NVARS(i) for(int i=0; i < GRMHD::nVar; i++)

namespace GRMHD {
	using namespace tensish; 
	using sym::delta; // Kronecker symbol
	
	/**
	 * These structs store the index positions of fields, i.e. Q.Dens=0.
	 * They are similiar to the VariableShortcuts which are generated
	 * by the ExaHyPE java toolkit. However, by defining our own we are
	 * independent and can make use of the functional notation.
	 **/
	struct ConservedIndices {
		// Conservatives
		static constexpr int   Dens = 0;
		static constexpr int   Si_lo = 1;
		static constexpr int   tau = 4;
		// Conserved (C2P invariant)
		static constexpr int   Bmag_up = 5;    /* length: 3 */
		static constexpr int   Phi = 8;
		// material parameters
		static constexpr int   lapse = 9;
		static constexpr int   shift_lo = 10;  /* length: 3 */
		static constexpr int   gam_lo = 13;    /* length: 6 */
		static constexpr int   detg = 19;
	};
	
	struct PrimitiveIndices { // no correspondence in ExaHyPE
		static constexpr int rho = 0;
		static constexpr int vec_up = 1; /* sic */
		static constexpr int press = 4;
	};
	
	/**
	* Storage for the conserved hydro variables, i.e. the tuple
	*   Q = (Dens, S_i, tau, B^i, phi).
	* The structure allows addressing the vector Q with names.
	* The template allows to specify types for the scalars and 3-vectors inside Q.
	* See the typedefs below for usable instances of the "Conserved" structure.
	*
	* In principle, this structure mimics the ExaHyPE generated variables structure
	* but in a much more usable manner (exploiting all the substructure of Q).
	**/
	template<typename vector, typename vector_lo, typename vector_up, typename scalar>
	struct Conserved {
		vector Q; // const double* const
		const GRMHD::ConservedIndices Qi;
		
		// Conserved Scalars
		scalar &Dens, &tau, &phi; // const double
		
		// Conserved vectors
		vector_lo Si; /* sic */
		vector_up Bmag; /* sic */
		
		constexpr Conserved(vector Q_) :
			Q(Q_),
			Qi(),
			// Conserved scalars
			Dens (Q[Qi.Dens]),
			tau  (Q[Qi.tau]),
			phi  (Q[Qi.Phi]),
			// Conserved vectors
			Si   (Q+Qi.Si_lo),
			Bmag (Q+Qi.Bmag_up)
			{}
		
		void toPrimitives(double* V); // C2P Cons2Prim operation
	};

	/**
	 * Access the Conserved Variables Q fully read-only without access to anything
	 * beyond what's inside Q (no further storage allocated except pointers to Q).
	 * That is, you cannot access S^i or B_i as we don't compute it. Therefore, we
	 * say this structure is just "shadowing" Q. You might want to use it like
	 * 
	 *    const ConservedConstShadow Q(Q_);
	 *    double foo = Q.Bmag.up(1) * Q.Dens;
	 * 
	 **/
	typedef GRMHD::Conserved<
		const double* const,
		const new_Lo<vec::const_shadow, const double* const>,
		const new_Up<vec::const_shadow, const double* const>,
		const double> ConservedConstShadow;

	/**
	 * The variable version of the ConservedConstShadow: Allows to change all
	 * entries of Q by their aliases, as well as directly changing the vectors, i.e.
	 *
	 *    ConservedVariableShadow Q(Q_);
	 *    Q.Bmag.up(2) = Q.Dens + 15;
	 * 
	 * But not:
	 * 
	 *    Q.Bmag.lo(2) = 15; 
	 * 
	 * That is, we dont' allow setting (or even getting) B_i or S^j. This can be
	 * helpful to avoid errors in the formulation.
	 **/
	typedef GRMHD::Conserved<
		double* const,
		new_Lo<vec::shadow, double* const>,
		new_Up<vec::shadow, double* const>,
		double> ConservedVariableShadow;
		
	/**
	 * Access to the full set of conserved variables in a read only fashion which
	 * however allows you to retrieve B_i and S^j, i.e. we have storage for these
	 * variables. These two vectors are the only non-constant members of the
	 * data structure as you need to fill them after construction.
	 **/
	typedef GRMHD::Conserved<
		const double* const,
		ConstLo<vec::stored<3>,vec::const_shadow>,
		ConstUp<vec::const_shadow, vec::stored<3>>,
		const double> ConservedConstFull;

	/**
	 * The variable version of ConservedConstFull: You can change all conserved
	 * variables as well access the B_i and S^j. Note that
	 * 
	 *   Q.Bmag.lo(2) = 15;
	 * 
	 * will not change the shadowed double* Q_.
	 **/
	typedef GRMHD::Conserved<
		double* const,
		StoredLo<vec::stored<3>, vec::shadow>,
		StoredUp<vec::shadow, vec::stored<3>>,
		double> ConservedVariableFull;
	

	/**
	 * Prmitive Variables.
	 * 
	 * They should only store the hydro variables.
	 * 
	 **/
	struct Primitives { // are always read only in our formulation, no need for P2C after ID
		double V[nVar];
		const GRMHD::PrimitiveIndices Vi;
		
		// Primitive Scalars
		const double &rho, &press;
		
		// Primitive vectors
		ConstUp<vec::const_shadow, vec::stored<3>> vel; /* vel is vec_up */
		
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

	template<typename vector_up, typename metric_lo>
	class ADMBase { // Material Parameters, always read only
		const GRMHD::ConservedIndices Ai;
	public:
		const double &alpha, &detg; // Scalars: Lapse, Determinant of g_ij
		vector_up beta; // Shift vector: (Conserved) Material parameter vector
		metric_lo gam;  // 3-Metric: (Conserved) Material parameter tensor
		
		constexpr ADMBase(const double* const Q) :
			Ai(),
			alpha(Q[Ai.lapse]),
			detg (Q[Ai.detg]),
			beta (Q+Ai.shift_lo),
			gam  (Q+Ai.gam_lo)
			{}
		
		void complete(/* some vector */) {
			// should (could) provide this.
		}
	};
	
	typedef GRMHD::ADMBase<
		// ConstUp<vec::const_shadow, vec::stored<3>>, // if you ever need beta_lo
		new_Up<vec::const_shadow, const double* const>, // if you never need beta_lo
		metric3> FullADMBase;
	
	// A version of the ADMVariables where beta_lo and the upper metric
	// are not recovered.
	typedef GRMHD::ADMBase<
		const new_Up<vec::const_shadow, const double* const>,
		const new_Lo<sym::const_shadow, const double* const>
		> BasicADMBase;


	// A type storing the gradients of the conserved vector in one direction.
	// Since it does not make sense to
	// retrieve the lower/upper components of gradients, we don't even allocate storage or
	// provide conversion strategies.
	struct Gradient : public BasicADMBase, public ConservedConstShadow {
		constexpr Gradient(const double* const Q) : BasicADMBase(Q), ConservedConstShadow(Q) {}
	};

	// Gradients in each direction: partial_i * something
	struct Gradients {
		const Gradient dir[DIMENSIONS];
		// Access an element in some direction
		const Gradient& operator[](int i) const { return dir[i]; }
		// Access an element, indicating that it's partial_i
		const Gradient& lo(int i) const { return dir[i]; }
		
		Gradients(const double* const gradQ, const int nVar) :
			#if DIMENSIONS == 3
			dir{gradQ+0, gradQ+nVar, gradQ+2*nVar}
			#elif DIMENSIONS == 2
			dir{gradQ+0, gradQ+nVar}
			#endif
			// function body:
			{}
	};

	/**
	 * The PDE structure wraps ("shadows") the read-only conserved vector Q in order to provide
	 * the possibility to compute the flux and the source terms. It inherits the conserved
	 * variables, primitive variables and ADM variables for convenient access when implementing
	 * the source and the flux.
	 *
	 * It used to provide storage for the 3-Energy momentum tensor in the past, but now this is
	 * stored directly in the methods because they need different forms (S_ij, S^ij or S^i_j)
	 * which we address differently in the current formalism.
	 **/
	struct PDE : public ConservedConstFull, public Primitives, public FullADMBase {
		// 3-Energy momentum tensor
		// Full<sym::stored<3>, sym::stored<3>> Sij;
		
		// Constraint damping constant
		const double damping_term_kappa = 5;
		
		PDE(const double*const Q_) :
			ConservedConstFull(Q_),
			Primitives(Q_),
			ADMBase(Q_)
			{ prepare(); }

		// Quantities for computing the energy momentum tensor
		double WW, BmagBmag, BmagVel, ptot;
		
		/// Prepares the Conserved B_i, S^i as well as the quantities neccessary to compute
		/// The energy momentum tensor.
		void prepare();
		
		/// This is the algebraic conserved flux. Chainable.
		void flux(/* const double* const Q, */ double** Fluxes);
		
		/// This is the fusedSource, but we just call it RightHandSide because it is on the RHS
		/// of the PDE. Chainable.
		void RightHandSide(/* const double* const Q, */const double* const gradQ_Data, double* Source_data);
		
		// You can recover these functions by just reworking what's in the fusedSource.
		//void nonConservativeProduct(/* const double* const Q, */ const double* const gradQ_Data, double* BgradQ_Data);
		//void algebraicSource(/* const double* const Q, */ double* Source_data);
		
		// currently trivial. Since we don't need to access Q, we provide static
		// access in order to avoid unneccessary boilerplate work.
		static void eigenvalues(const double* const Q, const int d,double* lambda);
	};

} // namespace GRMHD
#endif /* GRMHD_PDE_CPP */
