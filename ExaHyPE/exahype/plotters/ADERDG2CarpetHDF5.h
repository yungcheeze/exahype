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
 **/
 
#ifndef _EXAHYPE_PLOTTERS_ADERDG_2_CARPETHDF5_H_
#define _EXAHYPE_PLOTTERS_ADERDG_2_CARPETHDF5_H_

#include "exahype/plotters/Plotter.h"

namespace exahype {
  namespace plotters {
    class ADERDG2CarpetHDF5;
    class ADERDG2CarpetHDF5Writer;
  }
}

/**
 * <h2>Writing CarpetHDF5 files which are compatible to Cactus/EinsteinToolkit</h2>
 * 
 * This plotters produces files similar to the CarpetHDF5 file format. This file format is
 * produced by the https://www.carpetcode.org/ Carpet code for http://cactuscode.org/ Cactus
 * by the http://einsteintoolkit.org/ Einsteintoolkit.
 *
 * By using this plotter, users can use their Cactus postprocessing tools seamless with ExaHyPE
 * output without even knowing. There are also readers included into Visit and Amira for
 * opening the CarpetHDF5 file format.
 *
 * <h3>Current limitation</h3>
 *
 * Currently, we print only 2D. The extension to 3D is trivial.
 * 
 * Currently, we print only one field. Printing several fields in a single file is supported, but
 * not in the same table/mesh. We will extend this.
 * 
 * Currently, we print all timesteps in a single file and one file per MPI rank. This is also what
 * Cactus does, but we're free to split files once they get to large.
 *
 * <h3>How to build this plotter into your release</h3>
 * 
 * By default, this plotter is excluded from compiling. If you enable it in your spec file, the code
 * will stop at startup, throwing a message that this plotter is not supported.
 * 
 * In order to use this plotter, define "HDF5" and provide the HDF5 serial C++ header path and
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
 */
class exahype::plotters::ADERDG2CarpetHDF5 : public exahype::plotters::Plotter::Device {
 public:
  int           fileCounter;
  std::string   filename;
  int           order;
  int           solverUnknowns;
  int           writtenUnknowns;
  std::string   select;
  double        time;
  int           iteration;
  int           component;

  /**
   * Pimpl idiom: In order to avoid any HDF5 dependency all HDF5 logic is hidden inside this
   * class (instance).
   **/
  ADERDG2CarpetHDF5Writer* writer;

  static std::string getIdentifier();

  ADERDG2CarpetHDF5(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
  virtual ~ADERDG2CarpetHDF5();

  virtual void init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, const std::string& select);

  virtual void plotPatch(
        const int cellDescriptionsIndex,
        const int element);
  
  virtual void plotPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp);
  
  virtual void startPlotting( double time );
  virtual void finishPlotting();
  
  // TODO: These ADER interpolating routines should be in some generic library
  void interpolateCartesianPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double *u,
    double *mappedCell,
    double timeStamp
  );
};

#endif/* _EXAHYPE_PLOTTERS_ADERDG_2_CARPETHDF5_H_ */
