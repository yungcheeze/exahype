// Implementing all of the slicers in a single file to decrease
// compilation time, for the moment.

#include "Slicer.h"
#include "RegionSlicer.h"
#include "CartesianSlicer.h"

#include "exahype/Parser.h"

using namespace exahype::plotters;


CartesianSlicer::CartesianSlicer(const dvec& _req, const ivec& _active, int _baseDim) : 
		targetDim(_baseDim - tarch::la::sum(_active)),
		baseDim(_baseDim),
		req(_req),
		active(_active),
		activeAxes(-1),
		runningAxes(-1) {
		
	for(int i=0; i<DIMENSIONS; i++) {
		activeAxes(i) = disabled;
		runningAxes(i) = disabled;
		
		// This algorithm is crazy. Needed a lot of debugging with standalone
		// examples, but now its tested for DIM<=3.
		
		for(int j=i; j<DIMENSIONS; j++) { // forward check for actives
			if(active(j)) {
				activeAxes(i)=j;
				for(int k=0; k<i; k++) { // backward check if not already included
					if(activeAxes(k)==j)
						activeAxes(i)=disabled;
				}
				if(activeAxes(i)!=disabled)
					break;
			}
		}
		
		for(int j=i; j<DIMENSIONS; j++) { // forward check for actives
			if(!active(j)) {
				runningAxes(i)=j;
				for(int k=0; k<i; k++) { // backward check if not already included
					if(runningAxes(k)==j)
						runningAxes(i)=disabled;
				}
				if(runningAxes(i)!=disabled)
					break;
			}
		}
	}
}
	
CartesianSlicer* CartesianSlicer::fromSelectionQuery(const std::string& select) {
	dvec r;
	r(0) = Parser::getValueFromPropertyString(select, "x");
	r(1) = Parser::getValueFromPropertyString(select, "y");
	#if DIMENSIONS==3
	r(2) = Parser::getValueFromPropertyString(select, "z");
	#endif
	ivec ron;
	// NaN means the property was not present in the select string
	for(int i=0; i<DIMENSIONS; i++) { ron(i) = (r(i)!=r(i)) ? 0 : 1; }
	
	return new CartesianSlicer(r, ron);
}


std::ostream& operator<<(std::ostream &s,const CartesianSlicer& c) {
  s << "CartesianSlicer, Reducing Dimension " << c.baseDim << " to " << c.targetDim << ":\n";
  s << "   req = " << c.req << "\n";
  s << "   active = " << c.active << "\n";
  s << "   activeAxes = " << c.activeAxes << "\n";
  s << "   runningAxes = " << c.runningAxes << "\n";
  return s;
}


RegionSlicer* RegionSlicer::fromSelectionQuery(const std::string& select) {
	dvec regionOfInterestLeftBottomFront, regionOfInterestRightTopBack;
	double x;
	x = Parser::getValueFromPropertyString( select, "left" );
	regionOfInterestLeftBottomFront(0) = x!=x ? -std::numeric_limits<double>::max() : x; // "-", min
	x = Parser::getValueFromPropertyString( select, "bottom" );
	regionOfInterestLeftBottomFront(1) = x!=x ? -std::numeric_limits<double>::max() : x; // "-", min
	#if DIMENSIONS==3
	x = Parser::getValueFromPropertyString( select, "front" );
	regionOfInterestLeftBottomFront(2) = x!=x ? -std::numeric_limits<double>::max() : x; // "-", min
	#endif
	
	x = Parser::getValueFromPropertyString( select, "right" );
	regionOfInterestRightTopBack(0) = x!=x ? std::numeric_limits<double>::max() : x;
	x = Parser::getValueFromPropertyString( select, "top" );
	regionOfInterestRightTopBack(1) = x!=x ? std::numeric_limits<double>::max() : x;
	#if DIMENSIONS==3
	x = Parser::getValueFromPropertyString( select, "back" );
	regionOfInterestRightTopBack(2) = x!=x ? std::numeric_limits<double>::max() : x;
	#endif
	
	return new RegionSlicer(regionOfInterestLeftBottomFront, regionOfInterestRightTopBack);
}

Slicer* Slicer::bestFromSelectionQuery(const std::string& select) {
	// TODO: Implement by looking at which property strings are there.
	return RegionSlicer::fromSelectionQuery(select);
}

	
