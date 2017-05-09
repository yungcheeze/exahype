#ifndef EXAHYPE_PLOTTERS_SLICING_CARTESIAN
#define EXAHYPE_PLOTTERS_SLICING_CARTESIAN

namespace exahype {
  namespace plotters {
    struct CartesianSlicer;
  }
}

#include <algorithm>
#include "exahype/plotters/slicing/Slicer.h"
#include "tarch/la/Vector.h"
#include <iostream>

/**
 * A small auxilliary class to simplify cartesian plotting, ie.
 *  a) slicing on planes parallel to the xy, xz or yz plane
 *  b) slicing on lines parallel to the x, y or z axis.
 **/
struct exahype::plotters::CartesianSlicer : public exahype::plotters::Slicer {
	typedef tarch::la::Vector<DIMENSIONS, double> dvec;
	typedef tarch::la::Vector<DIMENSIONS, int> ivec;
	
	static constexpr int disabled = -1;
	const int targetDim; ///< The computed lower dimension. Typically 1 or 2.
	const int baseDim; ///< Actually DIMENSIONS. Typically 2 or 3.
	const dvec req; ///< The requested abscissa in each axis, for instance [NaN,NaN,42] for z=42 and [0,0,NaN] for x=x0, y=y0
	const ivec active; ///< (effective) boolean determining wether this axis is not NaN, for instance [0,0,1] for z=z0 and [1,1,0] for x=x0, y=y0
	ivec activeAxes; ///< A vector (starting from 0) indicating the active axis, for instance [2,-1,-1] for z=z0 and [0,1,-1] for x=x0, y=y0
	ivec runningAxes; ///< A vector indicating the free axis indices, for instance [0,1,-1] for z=z0 and [2,-1,-1] for x=x0, y=y0
	
	CartesianSlicer(const dvec& _req, const ivec& _active, int _baseDim=DIMENSIONS);
	static CartesianSlicer* fromSelectionQuery(const std::string& select);
	
	/// The inverse of active
	int running(int d) const { return active(d) ? 0 : 1; }

	/**
	 * Coarse patch selection criterion, as in all VTK plotters.
	 **/
	bool isPatchActive(const dvec& offsetOfPatch, const dvec& sizeOfPatch) const override {
		for(int axis=0; axis<baseDim; axis++) {
			if(active(axis) && (
			  (offsetOfPatch(axis)+sizeOfPatch(axis) < req(axis)) || // upper right bound smaller than requested coordinate
			  (offsetOfPatch(axis) <= req(axis))                     // lowe left bound smaller than requested coordinate
			)) {
				return false; // patch does not touch req(axis)
			}
		}
		return true;
	}

	/**
	 * Project point onto the slice, ie onto the 2D plane or onto a 1d line.
	 *
	 * The projection is not the shorted distance to the plane/line but a projection
	 * in terms of the coordinate axis, ie. replacing the coordinates. I didn't find
	 * a better name for this...
	 **/
	dvec project(dvec point) {
		for(int i=0; i<DIMENSIONS; i++) {
			if(active(i)) {
				point(i) = req(i);
			}
		}
		return point;
	}
	
	/**
	 * Project index onto 2D plane or 1D line in a way that it lives afterwards on
	 * the object.
	 **/
	ivec project(ivec index) {
		ivec ret(0);
		for(int i=0; i<DIMENSIONS; i++) {
			if(running(i)) {
				ret(i) = index(i);
			}
		}
		return ret;
	}
}; // class CartesianSlicer

std::ostream& operator<<(std::ostream &s,const exahype::plotters::CartesianSlicer& c);

#endif /* EXAHYPE_PLOTTERS_SLICING_CARTESIAN */
