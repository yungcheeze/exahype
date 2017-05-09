/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 *
 * @authors: Sven Koeppel
 **/

#ifndef _EXAHYPE_PLOTTERS_CARPET_HDF5_WRITER_
#define _EXAHYPE_PLOTTERS_CARPET_HDF5_WRITER_

#include "exahype/plotters/Plotter.h"
#include "kernels/KernelUtils.h" // idx::kernels

namespace exahype {
  namespace plotters {
    class CarpetHDF5Writer;
    class CarpetHDF5Slicer;
  }
}

#ifdef HDF5 // Only if H5 is present

// HDF5 library, only available if HDF5 is on the path
#include "H5Cpp.h"

/**
 * <h2>Writing CarpetHDF5 files which are compatible to Cactus/EinsteinToolkit</h2>
 * 
 * This writer produces files similar to the CarpetHDF5 file format. This file format is
 * produced by the https://www.carpetcode.org/ Carpet code for http://cactuscode.org/ Cactus
 * by the http://einsteintoolkit.org/ Einsteintoolkit.
 *
 * By using this writer, users can use their Cactus postprocessing tools seamless with ExaHyPE
 * output without even knowing. There are also readers included into Visit and Amira for
 * opening the CarpetHDF5 file format.
 *
 * <h3>Output format variants</h3>
 * 
 * The writer has two major flags to trigger the way how H5 files are produced:
 * 
 * <strong> oneFilePerTimestep </strong> mimics the way how ExaHyPE produces VTK files: There is
 * always (at least) one file per timestep. This is <em>not</em> the way how the CarpetHDF5 file
 * format works and you will have to join the HDF5 files with tools like h5join (shipped with
 * Carpet/Cactus) in order to read these files with ordinary writers. Nevertheless this is helpful
 * for debugging or immediately opening single timestep snapshots in Visit.
 * 
 * <strong> allUnknownsInOneFile </strong> allows you to reduce the number of H5 files by putting
 * the (vector of) all written unknowns in a single file. This also mimics how ExaHyPE's VTK files
 * look like. In contrast, in Cactus one typically has one hdf5 file per group/per physical field.
 * This allows easily to copy only the interesting fields from a supercomputer.
 *
 * <h3>General limitations</h3>
 * 
 * <strong> No global MPI writing </strong> Currently, each MPI rank produces its own file. In
 * contrast, in Cactus there is by default a single file accross all ranks. Again, files can be
 * easily merged with standard cactus tools.
 *
 * <strong> Structure of array vs. array of structures </strong>
 * The CarpetHDF5 file format directly resembles the closest way how to dump Cactus memory into files.
 * Cactus stores structures of arrays, for instance the four hydro variables {rho,velx,vely,velz,eps}
 * are stored as four big arrays in Cactus. In contrast, ExaHyPE stores arrays of structures, hence
 * in principle one big array where at each point there are the four variables. For this very reason,
 * producing CarpetHDF5 files comes with a lot of overhead, ie. a needless amount of over and over
 * repeated metadata and way too much components.
 *
 * By no means this is the fault of the HDF5 binary table format but really the way how CarpetHDF5
 * works or can be "maximally streched" to also fit ExaHyPE in.
 * 
 * <strong>Always vertex-centered data</strong>
 * As the Finite Differencing code Cactus uses by default vertex-centered data and stores them as
 * such, the data layout only supports vertex data and not cell-centered data. Maybe the actual
 * CarpetHDF5 format supports also cell-centered representation.
 * 
 * <h3>How to build the CarpetHDF5 plotters into your release</h3>
 * 
 * By default, the HDF5 plotters are excluded from compiling. If you enable it in your spec file, the code
 * will stop at startup, throwing a message that this plotter is not supported.
 * 
 * In order to use these plotters, define "HDF5" and provide the HDF5 serial C++ header path and
 * libraries. This can be done by adding to your project
 *
    PROJECT_CFLAGS+=-DHDF5
    PROJECT_CFLAGS+=-I/usr/include/hdf5/serial -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -Wdate-time -D_FORTIFY_SOURCE=2 -fstack-protector-strong
    PROJECT_LFLAGS+=-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -l hdf5_cpp /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a
    
 * This snippet works on a standard Ubuntu system with hdf5 developer libraries installed.
 * I extracted them from "h5cc -show" and "pkg-config --libs hdf5"
 * 
 * <h3>Connection to other plotters in ExaHyPE</h3>
 *
 * Actually this plotter is quite similar to the CartesianVTK format but just with regular
 * block patches. It is also similar to the new PeanoPatchFileFormat which was developed
 * in the same time.
 *
 * Thanks to Roland Haas for support in Munich and Frankfurt when coming up with the file format.
 *
 * @author Sven Köppel
 *
 **/
class exahype::plotters::CarpetHDF5Writer {
  typedef tarch::la::Vector<DIMENSIONS, double> dvec;
  tarch::logging::Log _log;

public:
  // information from the device::init() process
  const int           solverUnknowns; ///< The number of unknowns in the Solver (ie. number of PDEs)
  const int           writtenUnknowns; ///< The number of written out quantities.
  const std::string   basisFilename; ///< The filename prefix as it is common in ExaHyPE plotters
  const int           basisSize; ///< this is _orderPlusOne in ADERDG context and _numberOfCellsPerAxis-2*ghostZones in FV context
  const std::string   select; ///< A selection string for passing further parameters throught the ExaHyPE specification file.

  // set up during construction: Dimensional reduction
  int                 dim; ///< Dimension of the output generated. Do not change this. Setup by constructor.
  CarpetHDF5Slicer   *slicer; ///< Subslice, if present. Otherwise nullptr.
  kernels::index     *writtenCellIdx; ///< Index of a whole cell as in ExaHyPE
  kernels::index     *singleFieldIdx; ///< index of a whole component as in Carpet: Only one value per point
  int                 allFieldsSize; ///< as a service: basisSize^DIM * writtenUnkowns
  int                 singleFieldSize; ///< just basisSize^DIM

  // Things to be counted by this instance
  int                 component; ///< An internal counter of the components (=patches) written out in one plot cycle
  int                 iteration; ///< An internal counter of the number of plot cycle runned. It is kind of global.
  char**              writtenQuantitiesNames; // not const as we check for good names in constructor

  // HDF5 specific data types
  H5::H5File         *single_file; ///< Pointer to a single H5File. nullptr if not allUnknownsInOneFile
  H5::H5File        **seperate_files; ///< List of pointers to H5Files, a single nullptrs if allUnknownsInOneFile.
  H5::DataSpace       patch_space; ///< DataSpaces describing a component/patch: basisSize^D elements.
  H5::DataSpace       dtuple; ///< DataSpace describing a dim-dimensional tuple, ie dim numbers.

  const bool          oneFilePerTimestep; ///< Constant deciding whether to write one file (or series of files) per timestep
  const bool          allUnknownsInOneFile; ///< Constant deciding whether all unknowns should go into a single file or split files instead.

  /**
   * cf. also the documentation in the ADERDG2CarpetHDF5.h
   * 
   * oneFilePerTimestep: You might want to have this to view the result during computation
   *     as HDF5 is very lazy at writing. Note that writing out data this form is not compilant
   *     with most CarpetHDF5 readers (ie. the visit reader). You must join seperate HDF5 files afterwards
   *     manually.
   * 
   * allUnknownsInOneFile: Write different fields in a single H5 combined file. Typically for Cactus
   *     as structure of arrays is to write each unknown in its own file (ie. one file per physical field).
   *
   **/
  CarpetHDF5Writer(const std::string& _filename, int _basisSize, int _solverUnknowns, int _writtenUnknowns, const std::string& _select,
		   char** writtenQuantitiesNames, bool oneFilePerTimestep_=false, bool allUnknownsInOneFile_=false);

  void writeBasicGroup(H5::H5File* file);
  
  /**
   * We use a H5File** here in order to modify a passed pointer H5File* to a single
   * H5File.
   **/
  void openH5File(H5::H5File** file, std::string& filename);

  /**
   * Opens or switchs the currently active H5 file or the list of H5 files.
   * 
   **/
  void openH5(int iteration);
  
  void flushH5File(H5::H5File* file);
  
  void startPlotting(double time);
  void finishPlotting();

  /**
   * This is 2D and 3D, allows several unknowns, named fields and all that.
   * 
   * Possible problem: Local timestepping / each patch *could* have its own time.
   * Then the whole plotting approach of CarpetHDF5 fails and we have to collect
   * cells belonging to the same time somehow. Or we have to keep track of the
   * "iteration" number.
   **/
  void plotPatch(
      const dvec& offsetOfPatch, const dvec& sizeOfPatch, const dvec& dx,
      double* mappedCell, double timeStamp);
  
  void plotPatchForSingleUnknown(
      const dvec& offsetOfPatch, const dvec& sizeOfPatch, const dvec& dx,
      double* mappedCell, double timeStamp,
      int writtenUnknown, H5::H5File* target);

}; // class ADERDG2CarpetHDF5Impl

#include <algorithm>
/**
 * A small auxilliary class to simplify cartesian plotting, ie.
 *  a) slicing on planes parallel to the xy, xz or yz plane
 *  b) slicing on lines parallel to the x, y or z axis.
 * Used only in the CarpetHDF5Writer so far, thus the name.
 * 
 * TODO: Move the implementation of this class to CarpetHDF5Writer.cpp
 **/
struct exahype::plotters::CarpetHDF5Slicer {
	static constexpr int disabled = -1;
	typedef tarch::la::Vector<DIMENSIONS, double> dvec;
	typedef tarch::la::Vector<DIMENSIONS, int> ivec;
	
	const int targetDim; ///< The computed lower dimension. Typically 1 or 2.
	const int baseDim; ///< Actually DIMENSIONS. Typically 2 or 3.
	const dvec req; ///< The requested abscissa in each axis, for instance [NaN,NaN,42] for z=42 and [0,0,NaN] for x=x0, y=y0
	const ivec active; ///< (effective) boolean determining wether this axis is not NaN, for instance [0,0,1] for z=z0 and [1,1,0] for x=x0, y=y0
	ivec activeAxes; ///< A vector (starting from 0) indicating the active axis, for instance [2,-1,-1] for z=z0 and [0,1,-1] for x=x0, y=y0
	ivec runningAxes; ///< A vector indicating the free axis indices, for instance [0,1,-1] for z=z0 and [2,-1,-1] for x=x0, y=y0
	
	CarpetHDF5Slicer(const dvec& _req, const ivec& _active, int _baseDim=DIMENSIONS) : 
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
	
	/// The inverse of active
	int running(int d) const { return active(d) ? 0 : 1; }

	/**
	 * Coarse patch selection criterion, as in all VTK plotters.
	 **/
	bool shallIPlotPatch(const dvec& offsetOfPatch, const dvec& sizeOfPatch) {
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
}; // class CarpetHDF5Slicer

inline std::ostream& operator<<(std::ostream &s,const exahype::plotters::CarpetHDF5Slicer& c) {
  s << "Slicer, Reducing Dimension " << c.baseDim << " to " << c.targetDim << ":\n";
  s << "   req = " << c.req << "\n";
  s << "   active = " << c.active << "\n";
  s << "   activeAxes = " << c.activeAxes << "\n";
  s << "   runningAxes = " << c.runningAxes << "\n";
  return s;
}


#endif /* H5 */
#endif /* _EXAHYPE_PLOTTERS_CARPET_HDF5_WRITER_ */
