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
 
#include "exahype/plotters/ADERDG2AsciiVTK.h"
#include "tarch/parallel/Node.h"

// @todo 16/05/03:Dominic Etienne Charreir Plotter depends now on kernels.
// Should thus be placed in kernel module or the solver
// should provide a function that computes solution values
// at equidistant grid points
#include "kernels/DGMatrices.h"
#include "peano/utils/Loop.h"


exahype::plotters::ADERDG2AsciiVTK::ADERDG2AsciiVTK():
  _fileCounter(-1) {
}


std::string exahype::plotters::ADERDG2AsciiVTK::getIdentifier() {
  return "vtk::ascii";
}


void exahype::plotters::ADERDG2AsciiVTK::init(
  const std::string& filename,
  int                order,
  int                unknowns,
  const std::string& select
){
  _filename    = filename;
  _order       = order;
  _unknowns    = unknowns;
  _select      = select;
  _patchWriter = nullptr;




  double x;
  x = getValueFromPropertyString( select, "left" );
  _regionOfInterestLeftBottomFront(0) = x!=x ? std::numeric_limits<double>::min() : x;
  x = getValueFromPropertyString( select, "bottom" );
  _regionOfInterestLeftBottomFront(1) = x!=x ? std::numeric_limits<double>::min() : x;
#ifdef Din3
  x = getValueFromPropertyString( select, "front" );
  _regionOfInterestLeftBottomFront(2) = x!=x ? std::numeric_limits<double>::min() : x;
#endif


  x = getValueFromPropertyString( select, "right" );
  _regionOfInterestRightTopBack(0) = x!=x ? std::numeric_limits<double>::max() : x;
  x = getValueFromPropertyString( select, "top" );
  _regionOfInterestRightTopBack(1) = x!=x ? std::numeric_limits<double>::max() : x;
#ifdef Dim3
  x = getValueFromPropertyString( select, "back" );
  _regionOfInterestRightTopBack(2) = x!=x ? std::numeric_limits<double>::max() : x;
#endif
}


void exahype::plotters::ADERDG2AsciiVTK::startPlotting( double time ) {
  _fileCounter++;

  assertion( _patchWriter==nullptr );

  _patchWriter =
      new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
          new tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter());

  _gridWriter = _patchWriter->createSinglePatchWriter();
  _timeStampDataWriter = _patchWriter->createVertexDataWriter("time", 1);

  for (int i = 0; i < _unknowns; i++) {
    std::ostringstream identifier;
    identifier << "Q" << i;
    if ( _select.find(identifier.str())!=std::string::npos || _select=="{all}" ) {
      _vertexDataWriter.push_back(
        _patchWriter->createVertexDataWriter(identifier.str(), 1));
    }
  }

  assertion( _patchWriter!=nullptr );
  assertion( _gridWriter!=nullptr );
  assertion( _timeStampDataWriter!=nullptr );
}


void exahype::plotters::ADERDG2AsciiVTK::finishPlotting() {
  assertion( _patchWriter!=nullptr );
  assertion( _gridWriter!=nullptr );
  assertion( _timeStampDataWriter!=nullptr );

  _gridWriter->close();
  _timeStampDataWriter->close();
  for (auto& p : _vertexDataWriter) {
    p->close();
  }

  std::ostringstream snapshotFileName;
  snapshotFileName << _filename
#ifdef Parallel
                   << "-rank-" << tarch::parallel::Node::getInstance().getRank()
#endif
                   << "-" << _fileCounter << ".vtk";

  _patchWriter->writeToFile(snapshotFileName.str());

  for (auto& p : _vertexDataWriter) {
    delete p;
  }
  _vertexDataWriter.clear();

  delete _timeStampDataWriter;
  delete _gridWriter;
  delete _patchWriter;

  _patchWriter = nullptr;
  _timeStampDataWriter = nullptr;
  _gridWriter = nullptr;
}


exahype::plotters::ADERDG2AsciiVTK::~ADERDG2AsciiVTK() {
}


void exahype::plotters::ADERDG2AsciiVTK::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {
  if (
    tarch::la::allSmaller(_regionOfInterestLeftBottomFront,offsetOfPatch+sizeOfPatch)
    &&
    tarch::la::allGreater(_regionOfInterestRightTopBack,offsetOfPatch)
  ) {

  assertion( _patchWriter!=nullptr );
  assertion( _gridWriter!=nullptr );
  assertion( _timeStampDataWriter!=nullptr );

  int vertexIndex =
      _gridWriter->plotPatch(offsetOfPatch, sizeOfPatch, _order).first;

// @todo 16/05/03:Dominic Etienne Charrier
// This is the dirty way to perform the projection
// Feel free to modify.
// This is depending on the choice of basis/implementation.
// The equidistant grid projection should therefore be moved into the solver.
  dfor(i,_order+1) {
    _timeStampDataWriter->plotVertex(vertexIndex, timeStamp);

    int unknownPlotter = 0;
    for (int unknown=0; unknown < _unknowns; unknown++) {
      std::ostringstream identifier;
      identifier << "Q" << unknown;

      if ( _select.find(identifier.str())!=std::string::npos || _select.find("all")!=std::string::npos ) {
        double value = 0.0;
        dfor(ii,_order+1) { // Gauss-Legendre node indices
          int iGauss = peano::utils::dLinearisedWithoutLookup(ii,_order + 1);
          value += kernels::equidistantGridProjector1d[_order][ii(1)][i(1)] *
                   kernels::equidistantGridProjector1d[_order][ii(0)][i(0)] *
                   #ifdef Dim3
                   kernels::equidistantGridProjector1d[_order][ii(2)][i(2)] *
                   #endif
                   u[iGauss * _unknowns + unknown];
          assertion3(value == value, offsetOfPatch, sizeOfPatch, iGauss);
        }

        _vertexDataWriter[unknownPlotter]->plotVertex(vertexIndex, value);
        unknownPlotter++;
      }
    }
    vertexIndex++;
  }
  }
}
