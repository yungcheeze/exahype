/* Outsourcing from GRMHDSolver_*.cpp */
#ifndef GENERIC_BOUNDARY_CONDITIONS_WHICH_SHOULD_GO_SOMEWHERE_CENTRAL
#define GENERIC_BOUNDARY_CONDITIONS_WHICH_SHOULD_GO_SOMEWHERE_CENTRAL

/** Helper class to parse key value strings */
class StringMapView {
   std::string base;
   public:
    StringMapView(std::string _base) : base(_base) {}
    std::string getValueAsString(const std::string& key) const {
	std::size_t startIndex = base.find(key);
	startIndex = base.find(":", startIndex);
	std::size_t endIndex = base.find_first_of("}, \n\r", startIndex + 1);
	return base.substr(startIndex + 1, endIndex - startIndex - 1);
    }
    bool isValueValidString(const std::string& key) const {
	std::size_t startIndex = base.find(key);
	return (startIndex != std::string::npos);
    }
  };

/*
 * status of this file:
 * 
 * Is ugly with the ADERDG/FV basically replicated 90% of the code.
 *
 */

#define ADERDG_BOUNDARY_SIGNATURE \
  const double* const x,const double t,const double dt,const int faceIndex,const int d,\
  const double * const fluxIn,const double* const stateIn,\
  double *fluxOut, double* stateOut
#define ADERDG_BOUNDARY_CALL \
  x, t, dt, faceIndex, d, fluxIn, stateIn, fluxOut, stateOut

#define FV_BOUNDARY_SIGNATURE \
    const double* const x, const double t,const double dt,\
    const int faceIndex, const int d,\
    const double* const stateIn, double* stateOut
#define FV_BOUNDARY_CALL \
   x, t, dt, faceIndex, d, stateIn, stateOut


template <typename SolverType>
class ADERDGBoundaryConditions {
public:
	// face indices: 0 x=xmin 1 x=xmax, 2 y=ymin 3 y=ymax 4 z=zmin 5 z=zmax
	// corresponding 0-left, 1-right, 2-front, 3-back, 4-bottom, 5-top
	static constexpr int EXAHYPE_FACE_LEFT = 0;
	static constexpr int EXAHYPE_FACE_RIGHT = 1;
	static constexpr int EXAHYPE_FACE_FRONT = 2;
	static constexpr int EXAHYPE_FACE_BACK = 3;
	static constexpr int EXAHYPE_FACE_BOTTOM = 4;
	static constexpr int EXAHYPE_FACE_TOP = 5;
	
	SolverType* solver;
	
	typedef ADERDGBoundaryConditions<SolverType> me;
	typedef void (me::*boundarymethod)(ADERDG_BOUNDARY_SIGNATURE);
	
	// this could be improved by using boundarymethod side[6]
	// and having apply defined just as  side[faceIndex](...);.
	boundarymethod left, right, front, back, bottom, top;
	
	ADERDGBoundaryConditions(SolverType* _solver) : solver(_solver),
		left(nullptr), right(nullptr), front(nullptr), back(nullptr), bottom(nullptr), top(nullptr) {}
	
	bool allFacesDefined() {
		return (left != nullptr && right != nullptr && front != nullptr && back != nullptr && bottom != nullptr && top != nullptr);
	}
	
	/// Apply by looking up the saved boundarymethods
	void apply(ADERDG_BOUNDARY_SIGNATURE) {
		switch(faceIndex) {
			case EXAHYPE_FACE_LEFT:   (this->*left  )(ADERDG_BOUNDARY_CALL); break;
			case EXAHYPE_FACE_FRONT:  (this->*front )(ADERDG_BOUNDARY_CALL); break;
			case EXAHYPE_FACE_BOTTOM: (this->*bottom)(ADERDG_BOUNDARY_CALL); break;
			case EXAHYPE_FACE_RIGHT:  (this->*right )(ADERDG_BOUNDARY_CALL); break;
			case EXAHYPE_FACE_BACK:   (this->*back  )(ADERDG_BOUNDARY_CALL); break;
			case EXAHYPE_FACE_TOP:    (this->*top   )(ADERDG_BOUNDARY_CALL); break;
			default: throw std::runtime_error("Inconsistent face index");
		}
	}
	
	// just to highlight that apply() is the main method:
	void operator()(ADERDG_BOUNDARY_SIGNATURE) { apply(ADERDG_BOUNDARY_CALL); }
	
	/***************************************************************************
	 * Boundary Methods.
	 * 
	 * To add your own, do the following:
	 * 
	 *  1) Implement it with a method void yourName(ADERDG_BOUNDARY_SIGNATURE)
	 *  2) Add yourName to the parseFromString thing
	 * 
	 ***************************************************************************/
	
	#define BC_FROM_SPECFILE_SIDE(name) \
		if(constants.isValueValidString(#name)) { \
			std::string val = constants.getValueAsString(#name);\
			name = parseFromString(val);\
			if(name==nullptr) {\
				logError("setFromSpecFile", "Boundary condition at " #name ": Invalid method: '" << val << "'");\
			} else {\
				logInfo("setFromSpecFile", #name " boundary: " << val);\
			}\
		} else {\
			logError("setFromSpecFile", "Boundary condition at " #name " is not valid according to specfile parser.");\
		}

	/**
	 * MapView shall have methods
	 *    bool isValueValidString(std::string)
	 *    std:.string getvalueAsString(std::string)
	 * such as the Parser::ParserView class has.
	 * 
	 * @returns true in case of success, false otherwise
	 **/
	template<typename MapView>
	bool setFromSpecFile(MapView constants) {
		tarch::logging::Log _log("BoundaryCondition");
		BC_FROM_SPECFILE_SIDE(left);
		BC_FROM_SPECFILE_SIDE(right);
		BC_FROM_SPECFILE_SIDE(front);
		BC_FROM_SPECFILE_SIDE(back);
		BC_FROM_SPECFILE_SIDE(bottom);
		BC_FROM_SPECFILE_SIDE(top);
		return allFacesDefined();
	}

	/**
	 * assign a boundary method by some string value
	 * @returns pinter if value could be correctly parsed, in case of error nullptr.
	 **/
	static boundarymethod parseFromString(const std::string& value) {
		boundarymethod target=nullptr;
		if(value == "zero")
			target = &me::zero;
		if(value == "reflective" || value == "refl")
			target = &me::reflective;
		if(value == "copy" || value == "outflow")
			target = &me::outflow;
		if(value == "exact")
			target = &me::exact;
		return target;
		//return !(target==nullptr);
		//throw std::exception("Invalid boundary method"); // Todo: return a usable error, or std::optional<BoundaryMethod> or so.
	}

	static constexpr int nVar = SolverType::NumberOfVariables;
	static constexpr int nDim = DIMENSIONS;
	static constexpr int basisSize = SolverType::Order + 1;
	static constexpr int order = SolverType::Order;

	void zero(ADERDG_BOUNDARY_SIGNATURE) {
		std::memset(stateOut, 0, nVar * sizeof(double));
		std::memset(fluxOut, 0, nVar * sizeof(double));
	}
	
	/// Reflection/Hard wall BC
	void reflective(ADERDG_BOUNDARY_SIGNATURE) {
		// Reflection BC
		for(int m=0; m<nVar; m++) {
			stateOut[m] = stateIn[m];
			fluxOut[m] = -fluxIn[m];
		}
	}
	
	/// Outflow / Copy BC
	void outflow(ADERDG_BOUNDARY_SIGNATURE) {
		for(int m=0; m<nVar; m++) {
			stateOut[m] = stateIn[m];
			fluxOut[m] = +fluxIn[m];
		}
	}
	
	/// Exact time integrated BC
	void exact(ADERDG_BOUNDARY_SIGNATURE) {
		double Qgp[nVar], Fs[nDim][nVar], *F[nDim];
		for(int d=0; d<nDim; d++) F[d] = Fs[d];
		zero(ADERDG_BOUNDARY_CALL); // zeroise stateOut, fluxOut
		for(int i=0; i < basisSize; i++)  { // i == time
			const double weight = kernels::gaussLegendreWeights[order][i];
			const double xi = kernels::gaussLegendreNodes[order][i];
			double ti = t + xi * dt;

			solver->adjustPointSolution(x, weight/* not sure, just filling*/, ti, dt, Qgp);
			//id->Interpolate(x, ti, Qgp);
			solver->flux(Qgp, F);
			
			for(int m=0; m < nVar; m++) {
				stateOut[m] += weight * Qgp[m];
				fluxOut[m] += weight * Fs[d][m];
			}
		}
	}
};


template <typename SolverType>
class FVBoundaryConditions {
public:
	// face indices: 0 x=xmin 1 x=xmax, 2 y=ymin 3 y=ymax 4 z=zmin 5 z=zmax
	// corresponding 0-left, 1-right, 2-front, 3-back, 4-bottom, 5-top
	static constexpr int EXAHYPE_FACE_LEFT = 0;
	static constexpr int EXAHYPE_FACE_RIGHT = 1;
	static constexpr int EXAHYPE_FACE_FRONT = 2;
	static constexpr int EXAHYPE_FACE_BACK = 3;
	static constexpr int EXAHYPE_FACE_BOTTOM = 4;
	static constexpr int EXAHYPE_FACE_TOP = 5;
	
	SolverType* solver;
	
	typedef FVBoundaryConditions<SolverType> me;
	typedef void (me::*boundarymethod)(FV_BOUNDARY_SIGNATURE);
	
	// this could be improved by using boundarymethod side[6]
	// and having apply defined just as  side[faceIndex](...);.
	boundarymethod left, right, front, back, bottom, top;
	
	FVBoundaryConditions(SolverType* _solver) : solver(_solver),
		left(nullptr), right(nullptr), front(nullptr), back(nullptr), bottom(nullptr), top(nullptr) {}
	
	bool allFacesDefined() {
		return (left != nullptr && right != nullptr && front != nullptr && back != nullptr && bottom != nullptr && top != nullptr);
	}
	
	/// Apply by looking up the saved boundarymethods
	void apply(FV_BOUNDARY_SIGNATURE) {
		switch(faceIndex) {
			case EXAHYPE_FACE_LEFT:   (this->*left  )(FV_BOUNDARY_CALL); break;
			case EXAHYPE_FACE_FRONT:  (this->*front )(FV_BOUNDARY_CALL); break;
			case EXAHYPE_FACE_BOTTOM: (this->*bottom)(FV_BOUNDARY_CALL); break;
			case EXAHYPE_FACE_RIGHT:  (this->*right )(FV_BOUNDARY_CALL); break;
			case EXAHYPE_FACE_BACK:   (this->*back  )(FV_BOUNDARY_CALL); break;
			case EXAHYPE_FACE_TOP:    (this->*top   )(FV_BOUNDARY_CALL); break;
			default: throw std::runtime_error("Inconsistent face index");
		}
	}
	
	// just to highlight that apply() is the main method:
	void operator()(FV_BOUNDARY_SIGNATURE) { apply(FV_BOUNDARY_CALL); }
	
	/***************************************************************************
	 * Boundary Methods.
	 * 
	 * To add your own, do the following:
	 * 
	 *  1) Implement it with a method void yourName(FV_BOUNDARY_SIGNATURE)
	 *  2) Add yourName to the parseFromString thing
	 * 
	 ***************************************************************************/
	
	// similar as ADERDg
	template<typename MapView>
	bool setFromSpecFile(MapView constants) {
		tarch::logging::Log _log("FV-BoundaryCondition");
		BC_FROM_SPECFILE_SIDE(left);
		BC_FROM_SPECFILE_SIDE(right);
		BC_FROM_SPECFILE_SIDE(front);
		BC_FROM_SPECFILE_SIDE(back);
		BC_FROM_SPECFILE_SIDE(bottom);
		BC_FROM_SPECFILE_SIDE(top);
		return allFacesDefined();
	}

	/**
	 * assign a boundary method by some string value
	 * @returns pinter if value could be correctly parsed, in case of error nullptr.
	 **/
	static boundarymethod parseFromString(const std::string& value) {
		boundarymethod target=nullptr;
		if(value == "zero")
			target = &me::zero;
		if(value == "reflective" || value == "refl")
			target = &me::reflective;
		if(value == "copy" || value == "outflow")
			target = &me::outflow;
		if(value == "exact")
			target = &me::exact;
		return target;
		//return !(target==nullptr);
		//throw std::exception("Invalid boundary method"); // Todo: return a usable error, or std::optional<BoundaryMethod> or so.
	}
	
	static constexpr int nVar = SolverType::NumberOfVariables;
	static constexpr int nDim = DIMENSIONS;

	void zero(FV_BOUNDARY_SIGNATURE) {
		std::memset(stateOut, 0, nVar * sizeof(double));
		//std::memset(fluxOut, 0, nVar * sizeof(double));
	}
	
	/// Reflection/Hard wall BC
	void reflective(FV_BOUNDARY_SIGNATURE) {
		// Reflection BC
		for(int m=0; m<nVar; m++) {
			stateOut[m] = stateIn[m];
			//fluxOut[m] = -fluxIn[m];
			// dafuq, how we get negative fluxes in FV?
		}
	}
	
	/// Outflow / Copy BC
	void outflow(FV_BOUNDARY_SIGNATURE) {
		for(int m=0; m<nVar; m++) {
			stateOut[m] = stateIn[m];
			//fluxOut[m] = +fluxIn[m];
		}
	}
	
	/// Exact time integrated BC
	void exact(FV_BOUNDARY_SIGNATURE) {
		double snan = std::numeric_limits<double>::signaling_NaN();
		solver->adjustSolution(x, snan/* not sure, just filling*/, t, dt, stateOut);
	}
};

#endif /* GENERIC_BOUNDARY_CONDITIONS_WHICH_SHOULD_GO_SOMEWHERE_CENTRAL */
