#ifndef EXAHYPE_PLOTTERS_SLICING_BASE
#define EXAHYPE_PLOTTERS_SLICING_BASE

namespace exahype {
  namespace plotters {
    struct Slicer;
  }
}

#include "peano/utils/Dimensions.h"
#include "tarch/la/Vector.h"
#include <iostream>

/**
 * Abstract slicing base class for <b>patchwise</b> slicing.
 **/
struct exahype::plotters::Slicer {
	typedef tarch::la::Vector<DIMENSIONS, double> dvec;
	typedef tarch::la::Vector<DIMENSIONS, int> ivec;
	
	virtual bool isPatchActive(const dvec& offsetOfPatch, const dvec& sizeOfPatch) const = 0;
	
	/**
	 * Kind of a registry: Make an educated guess which slicer is fitting.
	 * As Slicer is abstract, this allocates a subclass.
	 **/
	Slicer* bestFromSelectionQuery(const std::string& select);

};

#endif /* EXAHYPE_PLOTTERS_SLICING_BASE */
