#ifndef _OH_REALLY_COME_ON
#define _OH_REALLY_COME_ON

/**
 * Really is a template free compile time library to provide most of FORTRANs n-dimensional
 * arrays without providing a fully-fledged tensor library.
 * It is optimized for short compile times and no runtime overhead.
 * 
 * This copy of Really is a stripped down version which only contains the n-dimensional
 * "shadow" data container but no lazy-evaluated algebra (MATMUL, etc.). That means that
 * all FORTRAN-ness have been removed.
 * 
 * The name of "REALLY" comes from the fact that FORTRAN calls doubles "REAL".
 * 
 * Written by SvenK in 2017.
 **/
namespace Really {
	
	/******************************************************************************
	 * 
	 * INDEXING AND PRIMITIVE BASICS
	 * 
	 ******************************************************************************/
	
	struct range {
		const int a, b;
		range(const int a_, const int b_) : a(a_), b(b_) {}
	};
	
	typedef range slice;

	/**
	 * An index class is just a fancy way of writing an array.
	 * 
	 **/
	struct index {
		#define INDICES(type, suffix) \
			type i0 suffix,\
			type i1 suffix,\
			type i2 suffix,\
			type i3 suffix,\
			type i4 suffix,\
			type i5 suffix
		#define INDICES_LIST INDICES(,)
		
		#define INDICES_INT INDICES(int,)
		#define INDICES_EQ1 INDICES(int,=1)
		
		const int INDICES_LIST;
		// todo: define with j.
		index(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) :
			i0(j0), i1(j1), i2(j2), i3(j3), i4(j4), i5(j5) {}
	};
	
	/**
	* This is a single successor class for the idx2, idx3, idx4, idx5, idx6 classes.
	* It works basically like idx6. If you work with less than 6 dimensions, nothing
	* will change and the compiler optimizes away the unused variables.
	* 
	* In contrast to the shape classes below, this new shape class can also compute
	* the inverse shape.
	* 
	* shape was index in kernelutils.
	**/
	struct shape {
		// add further indices if you want to support more than maximal 6 indices.
		const int i0, i1, i2, i3, i4, i5; // Lenghts
		const int b0, b1, b2, b3, b4, b5; // Basis
		const int size;
		
		constexpr shape(int j0=1, int j1=1, int j2=1, int j3=1, int j4=1, int j5=1) :
			i0(j0), i1(j1), i2(j2), i3(j3), i4(j4), i5(j5),
			b0(i1 * i2 * i3 * i4 * i5),
			b1(i2 * i3 * i4 * i5),
			b2(i3 * i4 * i5),
			b3(i4 * i5),
			b4(i5),
			b5(1),
			size(j0*j1*j2*j3*j4*j5) {}
		
		/**
		* Compute a single shape ("supershape", global shape, ...) from the tuples.
		**/
		int get(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) const {
			return j0 * b0 + j1 * b1 + j2 * b2 + j3 * b3 + j4 * b4 + j5 * b5;
		}
		
		/**
		* Inverse shape: Get the shape tuple for a given shape. This is the inverse
		* function of get(), so such a call will not do anything:
		*    myidx.rev(myidx.get(a,b,c...), a, b, c, ...)
		**/
		void rev(int pos, int& j0, int& j1, int& j2, int& j3, int& j4, int& j5) const {
			// The algorithm exploits the short notation for int a,b that
			// a/b == std::floor(a/b). The c's are the contributions at position
			// i while the j's are the digits for position i. If you want to
			// increase the overall amount of positions, extend logically.
			j0 = (pos) / b0;
			j1 = (pos-j0*b0) / b1;
			j2 = (pos-j0*b0-j1*b1) / b2;
			j3 = (pos-j0*b0-j1*b1-j2*b2) / b3;
			j4 = (pos-j0*b0-j1*b1-j2*b2-j3*b3) / b4;
			j5 = (pos-j0*b0-j1*b1-j2*b2-j3*b3-j4*b4) / b5;
		}
		
		/* As there are no reference default values, ie. we cannot write
		*   void rev(int pos, int& j0=nullptr) or
		*   void rev(int pos, int& j1=int(), ...)
		* we do it with pointers. Now you can write instead
		*   int i, j, k;
		*   myidx.rev(position, &i, &j, &k);
		* which is much more convenient than tracking the non-used indices
		* for your own.
		*/
		void rev(int pos, int* const j0=nullptr, int* const j1=nullptr, int* const j2=nullptr, int* const j3=nullptr, int* const j4=nullptr, int* const j5=nullptr) const {
			int a0=0, a1=0, a2=0, a3=0, a4=0, a5=0; // default storages
			rev(pos, j0?*j0:a0, j1?*j1:a1, j2?*j2:a2, j3?*j3:a3, j4?*j4:a4, j5?*j5:a5);
		}
		
		// syntactic sugar:
		int operator()(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) const {
			return get(j0, j1, j2, j3, j4, j5);
		}
		
		/**
		* Checks if the given indices are in the valid range. Works also if assertions are
		* not enabled. Can be handy to check access to variables.
		**/
		bool check(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) const {
			return (j0 < i0) && (j1 < i1) && (j2 < i2) && (j3 < i3) && (j4 < i4) && (j5 < i5);
		}
		
		/**
		* Allows to print shape tuples.
		**/
		/*
		static std::string strIndex(int min=0, int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0);
		
		/// Some shape to string
		std::string getStr(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) const;
		std::string revStr(int pos) const;

		/// Object to string
		std::string toString() const;
		*/
	};
	
	//i0, i1, i2, i3, i4, i5 // usw.
	#define INDEX_LOOP(idx) \
		for(int i0=0; i0<idx.i0; i0++)\
		for(int i1=0; i1<idx.i1; i1++)\
		for(int i2=0; i2<idx.i2; i2++)\
		for(int i3=0; i3<idx.i3; i3++)\
		for(int i4=0; i4<idx.i4; i4++)\
		for(int i5=0; i5<idx.i5; i5++)  // no \ at end

	#define IDX_LOOP(I, idx) \
		for(idx I(0); i0<idx.i0; i0++)\
		for(int i1=0; i1<idx.i1; i1++)\
		for(int i2=0; i2<idx.i2; i2++)\
		for(int i3=0; i3<idx.i3; i3++)\
		for(int i4=0; i4<idx.i4; i4++)\
		for(int i5=0; i5<idx.i5; i5++)  // no \ at end

	/**
	* The same as array, just without local storage but instead shadowing data which is not owned.
	* This can be handy for accessing big arrays stored somewhere else. A tiny example is
	* 
	*	double storage[5*7*9];
	*	dshadow convenient(storage, 5, 7, 9);
	*	convenient(2,3,8) = 10;
	*      double val = convenient(2,3,8);
	*
	**/
	struct shadow  {
		typedef shadow me;
		double *data;
		const shape idx;
		
		shadow(double* storage, const shape& i) : data(storage), idx(i) {}
		shadow(double* storage, int j0=1, int j1=1, int j2=1, int j3=1, int j4=1, int j5=1) : data(storage), idx(j0,j1,j2,j3,j4,j5) {}
		
		#define GETTER(getter_name) \
			double& getter_name(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) const {\
				return data[idx(j0, j1, j2, j3, j4, j5)];\
			}
		
		// Three ways to access something
		GETTER(get);
		GETTER(operator());
		//GETTER(operator[]);
		
		// we skip returning *this in order to maintain better extendability
		// could be done with another argument
		#define SETTER(setter_name, rhs_argument, rhs) \
			shadow&	 setter_name (rhs_argument) {\
				INDEX_LOOP(idx) {\
					get(INDICES_LIST) = rhs; \
				}\
				return *this;\
			}

		SETTER(set, const shadow& rhs, rhs(INDICES_LIST) );
		SETTER(set, const double& rhs, rhs ); // a particular scalar

		// non-inheriting class methods:
		// Copy consructor and assignments
		#define SETTER_EQ(typename) \
			typename& operator=(const typename& rhs) { return (typename&)set(rhs); }\
			typename& operator=(const double& rhs) { return (typename&)set(rhs); }
		
		SETTER_EQ(shadow);
		shadow(const shadow& rhs) { set(rhs); }
		
		me& operator+(const me& b) {
			INDEX_LOOP(idx) {
				get(INDICES_LIST) += b(INDICES_LIST);
			}
			return *this;
		}
	};
	
	// Const read only version
	struct const_shadow {
		typedef shadow me;
		const double* const data;
		const shape idx;
		
		const_shadow(const double* const storage, const shape& i) : data(storage), idx(i) {}
		const_shadow(const double* const storage, int j0=1, int j1=1, int j2=1, int j3=1, int j4=1, int j5=1) : data(storage), idx(j0,j1,j2,j3,j4,j5) {}
		
		#define CONST_GETTER(getter_name) \
			const double& getter_name(int j0=0, int j1=0, int j2=0, int j3=0, int j4=0, int j5=0) const {\
				return data[idx(j0, j1, j2, j3, j4, j5)];\
			}
		
		// Three ways to access something
		CONST_GETTER(get);
		CONST_GETTER(operator());
		//GETTER(operator[]);
	};

} // namespace

#endif /* _OH_REALLY_COME_ON */
