#ifndef _EXAHYPE_PLOTTERS_CARPET_HDF5_WRITER_
#define _EXAHYPE_PLOTTERS_CARPET_HDF5_WRITER_

#include "exahype/plotters/Plotter.h"
#include "kernels/KernelUtils.h" // idx::kernels

namespace exahype {
  namespace plotters {
    class CarpetHDF5Writer;
  }
}

#ifdef HDF5 // Only if H5 is present

// HDF5 library, only available if HDF5 is on the path
#include "H5Cpp.h"

class exahype::plotters::CarpetHDF5Writer {
  typedef tarch::la::Vector<DIMENSIONS, double> dvec;
  tarch::logging::Log _log;

public:
  // information from the device::init() process
  const int           solverUnknowns;
  const int           writtenUnknowns;
  const std::string   basisFilename;
  const int           basisSize; ///< this is _orderPlusOne in ADERDG context and _numberOfCellsPerAxis-2*ghostZones in FV context
  const std::string   select;

  // Things to be counted by this instance
  int                 component;
  int                 iteration;
  char**              writtenQuantitiesNames; // not const as we check for good names in constructor
  const kernels::index writtenCellIdx, singleFieldIdx;

  // HDF5 specific data types
  H5::H5File *single_file, **seperate_files;
  H5::DataSpace patch_space, dtuple;

  const bool oneFilePerTimestep, allUnknownsInOneFile;

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

#endif /* H5 */
#endif /* _EXAHYPE_PLOTTERS_CARPET_HDF5_WRITER_ */
