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
 
#include "FiniteVolumes2VTKAscii.h"
#include "ADERDG2VTKAscii.h"
#include "tarch/parallel/Node.h"

// @todo 16/05/03:Dominic Etienne Charreir Plotter depends now on kernels.
// Should thus be placed in kernel module or the solver
// should provide a function that computes solution values
// at equidistant grid points
#include "kernels/DGMatrices.h"
#include "peano/utils/Loop.h"

tarch::logging::Log exahype::plotters::FiniteVolumes2VTKAscii::_log("exahype::plotters::FiniteVolumes2VTKAscii");

exahype::plotters::FiniteVolumes2VTKAscii::FiniteVolumes2VTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
  Device(postProcessing),
  _fileCounter(-1) {
}


std::string exahype::plotters::FiniteVolumes2VTKAscii::getIdentifier() {
  return ADERDG2VTKAscii::getIdentifier();
}


void exahype::plotters::FiniteVolumes2VTKAscii::init(
  const std::string& filename,
  int                numberOfCellsPerAxis,
  int                unknowns,
  int writtenUnknowns,
  const std::string& select
){
  _filename             = filename;
  _numberOfCellsPerAxis = numberOfCellsPerAxis;
  _unknowns             = unknowns;
  _select               = select;
  _patchWriter          = nullptr;

  double x;
  x = Parser::getValueFromPropertyString( select, "left" );
  _regionOfInterestLeftBottomFront(0) = x!=x ? std::numeric_limits<double>::min() : x;
  x = Parser::getValueFromPropertyString( select, "bottom" );
  _regionOfInterestLeftBottomFront(1) = x!=x ? std::numeric_limits<double>::min() : x;
#ifdef Din3
  x = Parser::getValueFromPropertyString( select, "front" );
  _regionOfInterestLeftBottomFront(2) = x!=x ? std::numeric_limits<double>::min() : x;
#endif


  x = Parser::getValueFromPropertyString( select, "right" );
  _regionOfInterestRightTopBack(0) = x!=x ? std::numeric_limits<double>::max() : x;
  x = Parser::getValueFromPropertyString( select, "top" );
  _regionOfInterestRightTopBack(1) = x!=x ? std::numeric_limits<double>::max() : x;
#ifdef Dim3
  x = Parser::getValueFromPropertyString( select, "back" );
  _regionOfInterestRightTopBack(2) = x!=x ? std::numeric_limits<double>::max() : x;
#endif
}


void exahype::plotters::FiniteVolumes2VTKAscii::startPlotting( double time ) {
  _fileCounter++;

  assertion( _patchWriter==nullptr );

  _patchWriter =
      new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
          new tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter());

  _gridWriter = _patchWriter->createSinglePatchWriter();
  _timeStampDataWriter = _patchWriter->createCellDataWriter("time", 1);

  for (int i = 0; i < _unknowns; i++) {
    std::ostringstream identifier;
    identifier << "Q" << i;
    if ( _select.find(identifier.str())!=std::string::npos || _select=="{all}" ) {
      _cellDataWriter.push_back(_patchWriter->createCellDataWriter(identifier.str(), 1));
    }
  }

  assertion( _patchWriter!=nullptr );
  assertion( _gridWriter!=nullptr );
  assertion( _timeStampDataWriter!=nullptr );
}


void exahype::plotters::FiniteVolumes2VTKAscii::finishPlotting() {
  assertion( _patchWriter!=nullptr );
  assertion( _gridWriter!=nullptr );
  assertion( _timeStampDataWriter!=nullptr );

  _gridWriter->close();
  _timeStampDataWriter->close();
  for (auto& p : _cellDataWriter) {
    p->close();
  }

  std::ostringstream snapshotFileName;
  snapshotFileName << _filename
#ifdef Parallel
                   << "-rank-" << tarch::parallel::Node::getInstance().getRank()
#endif
                   << "-" << _fileCounter << ".vtk";

  _patchWriter->writeToFile(snapshotFileName.str());

  for (auto& p : _cellDataWriter) {
    delete p;
  }
  _cellDataWriter.clear();

  delete _timeStampDataWriter;
  delete _gridWriter;
  delete _patchWriter;

  _patchWriter = nullptr;
  _timeStampDataWriter = nullptr;
  _gridWriter = nullptr;
}


exahype::plotters::FiniteVolumes2VTKAscii::~FiniteVolumes2VTKAscii() {
}


void exahype::plotters::FiniteVolumes2VTKAscii::plotPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
  double timeStamp) {
  if (
    tarch::la::allSmaller(_regionOfInterestLeftBottomFront,offsetOfPatch+sizeOfPatch)
    &&
    tarch::la::allGreater(_regionOfInterestRightTopBack,offsetOfPatch)
  ) {
    logDebug("plotPatch(...)","offset of patch: "<<offsetOfPatch
    <<", size of patch: "<<sizeOfPatch
    <<", time stamp: "<<timeStamp);

    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );
    assertion( _timeStampDataWriter!=nullptr );

    int cellIndex = _gridWriter->plotPatch(offsetOfPatch, sizeOfPatch, _numberOfCellsPerAxis).second;

    // @todo 16/05/03:Dominic Etienne Charrier
    // This is depending on the choice of basis/implementation.
    // The equidistant grid projection should therefore be moved into the solver.
    dfor(i,_numberOfCellsPerAxis) {
      _timeStampDataWriter->plotCell(cellIndex, timeStamp);

      int unknownPlotter = 0;
      for (int unknown=0; unknown < _unknowns; unknown++) {
        std::ostringstream identifier;
        identifier << "Q" << unknown;

        if ( _select.find(identifier.str())!=std::string::npos || _select.find("all")!=std::string::npos ) {
          double value = u[peano::utils::dLinearisedWithoutLookup(i,_numberOfCellsPerAxis)];
          assertion3(!std::isnan(value),unknown,_unknowns,_numberOfCellsPerAxis);
          _cellDataWriter[unknownPlotter]->plotCell(cellIndex, value);
          unknownPlotter++;
        }
      }
      cellIndex++;
    }
  }
}
