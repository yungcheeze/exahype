#ifndef EXAHYPE_PLOTTERS_SLICING_REGION
#define EXAHYPE_PLOTTERS_SLICING_REGION

namespace exahype {
  namespace plotters {
    struct RegionSlicer;
  }
}

#include "exahype/plotters/slicing/Slicer.h"
#include "tarch/la/Vector.h"
#include <string>

struct exahype::plotters::RegionSlicer : public exahype::plotters::Slicer {
	typedef tarch::la::Vector<DIMENSIONS, double> dvec;
	typedef tarch::la::Vector<DIMENSIONS, int> ivec;
	
	const dvec  _regionOfInterestLeftBottomFront, _regionOfInterestRightTopBack;
	
	RegionSlicer(const dvec& regionOfInterestLeftBottomFront, const dvec& regionOfInterestRightTopBack) :
		_regionOfInterestLeftBottomFront(regionOfInterestLeftBottomFront),
		_regionOfInterestRightTopBack(regionOfInterestRightTopBack) {}
		
	static RegionSlicer* fromSelectionQuery(const std::string& select);
	
	bool isPatchActive(const dvec& offsetOfPatch, const dvec& sizeOfPatch) const override {
		return
			tarch::la::allSmaller(_regionOfInterestLeftBottomFront,offsetOfPatch+sizeOfPatch)
			&&
			tarch::la::allGreater(_regionOfInterestRightTopBack,offsetOfPatch);
	}
};

#endif /* EXAHYPE_PLOTTERS_SLICING_REGION */ 
