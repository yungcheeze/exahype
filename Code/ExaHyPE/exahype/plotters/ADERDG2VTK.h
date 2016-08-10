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
 
#ifndef _EXAHYPE_PLOTTERS_ADERDG_2_VTK_H_
#define _EXAHYPE_PLOTTERS_ADERDG_2_VTK_H_

#include "exahype/plotters/Plotter.h"

#include "tarch/plotter/griddata/blockstructured/PatchWriterUnstructured.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"

namespace exahype {
  namespace plotters {
    class ADERDG2VTK;

    class ADERDG2CartesianVerticesVTKAscii;
    class ADERDG2CartesianVerticesVTKBinary;
    class ADERDG2CartesianCellsVTKAscii;
    class ADERDG2CartesianCellsVTKBinary;

    class ADERDG2LegendreVerticesVTKAscii;
    class ADERDG2LegendreVerticesVTKBinary;
    class ADERDG2LegendreCellsVTKAscii;
    class ADERDG2LegendreCellsVTKBinary;
  }
}

/**
 * Common VTK class. Usually not used directly but through one of the subclasses.
 */
class exahype::plotters::ADERDG2VTK: public exahype::plotters::Plotter::Device {
 private:
  int           _fileCounter;
  const bool    _isBinary;
  const bool    _isCartesian;
  const bool    _plotCells;
  std::string   _filename;
  int           _order;
  int           _solverUnknowns;
  int           _writtenUnknowns;
  std::string   _select;


  tarch::la::Vector<DIMENSIONS, double>  _regionOfInterestLeftBottomFront;
  tarch::la::Vector<DIMENSIONS, double>  _regionOfInterestRightTopBack;

  tarch::plotter::griddata::blockstructured::PatchWriterUnstructured*        _patchWriter;
  tarch::plotter::griddata::blockstructured::PatchWriter::SinglePatchWriter* _gridWriter;

  tarch::plotter::griddata::Writer::VertexDataWriter*  _vertexTimeStampDataWriter;
  tarch::plotter::griddata::Writer::CellDataWriter*    _cellTimeStampDataWriter;
  tarch::plotter::griddata::Writer::VertexDataWriter*  _vertexDataWriter;
  tarch::plotter::griddata::Writer::CellDataWriter*    _cellDataWriter;

 public:
  ADERDG2VTK(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing, bool isBinary, bool isCartesian, bool plotCells);
  virtual ~ADERDG2VTK();

  virtual void init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, const std::string& select);

  virtual void plotPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp);

  virtual void startPlotting( double time );
  virtual void finishPlotting();
};


class exahype::plotters::ADERDG2CartesianVerticesVTKAscii: public exahype::plotters::ADERDG2VTK {
  public:
    static std::string getIdentifier();
    ADERDG2CartesianVerticesVTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
};


class exahype::plotters::ADERDG2CartesianVerticesVTKBinary: public exahype::plotters::ADERDG2VTK {
  public:
    static std::string getIdentifier();
    ADERDG2CartesianVerticesVTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
};


class exahype::plotters::ADERDG2CartesianCellsVTKAscii: public exahype::plotters::ADERDG2VTK {
  public:
    static std::string getIdentifier();
    ADERDG2CartesianCellsVTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
};


class exahype::plotters::ADERDG2CartesianCellsVTKBinary: public exahype::plotters::ADERDG2VTK {
  public:
    static std::string getIdentifier();
    ADERDG2CartesianCellsVTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
};



class exahype::plotters::ADERDG2LegendreVerticesVTKAscii: public exahype::plotters::ADERDG2VTK {
  public:
    static std::string getIdentifier();
    ADERDG2LegendreVerticesVTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
};


class exahype::plotters::ADERDG2LegendreVerticesVTKBinary: public exahype::plotters::ADERDG2VTK {
  public:
    static std::string getIdentifier();
    ADERDG2LegendreVerticesVTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
};


class exahype::plotters::ADERDG2LegendreCellsVTKAscii: public exahype::plotters::ADERDG2VTK {
  public:
    static std::string getIdentifier();
    ADERDG2LegendreCellsVTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
};


class exahype::plotters::ADERDG2LegendreCellsVTKBinary: public exahype::plotters::ADERDG2VTK {
  public:
    static std::string getIdentifier();
    ADERDG2LegendreCellsVTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
};


#endif
