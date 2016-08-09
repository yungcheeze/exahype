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
 
#include "ADERDG2VTKBinary.h"
#include "tarch/parallel/Node.h"
// @todo 16/05/03:Dominic Etienne Charreir Plotter depends now on kernels.
// Should thus be placed in kernel module or the solver
// should provide a function that computes solution values
// at equidistant grid points
#include "kernels/DGMatrices.h"
#include "peano/utils/Loop.h"


exahype::plotters::ADERDG2VTKBinary::ADERDG2VTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
  Device(postProcessing),
  _fileCounter(-1) {
}


std::string exahype::plotters::ADERDG2VTKBinary::getIdentifier() {
  return "vtk::Cartesian::binary";
}


void exahype::plotters::ADERDG2VTKBinary::init(
  const std::string& filename,
  int                orderPlusOne,
  int                unknowns,
  int                writtenUnknowns,
  const std::string& select
) {
  _filename          = filename;
  _order             = orderPlusOne-1;
  _solverUnknowns    = unknowns;
  _select            = select;
  _patchWriter       = nullptr;
  _writtenUnknowns   = writtenUnknowns;

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


void exahype::plotters::ADERDG2VTKBinary::startPlotting( double time ) {
  _fileCounter++;

  assertion( _patchWriter==nullptr );

  if (_writtenUnknowns>0) {
    _patchWriter =
        new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
            new tarch::plotter::griddata::unstructured::vtk::
                VTKBinaryFileWriter());

    _gridWriter          = _patchWriter->createSinglePatchWriter();
    _timeStampDataWriter = _patchWriter->createVertexDataWriter("time", 1);
    _vertexDataWriter    = _patchWriter->createVertexDataWriter("Q", _writtenUnknowns);

    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );
    assertion( _timeStampDataWriter!=nullptr );
  }

  _postProcessing->startPlotting( time );
}


void exahype::plotters::ADERDG2VTKBinary::finishPlotting() {
  _postProcessing->finishPlotting();

  if (_writtenUnknowns>0) {
    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );
    assertion( _timeStampDataWriter!=nullptr );

    _gridWriter->close();
    _timeStampDataWriter->close();
    _vertexDataWriter->close();

    std::ostringstream snapshotFileName;
    snapshotFileName << _filename
    #ifdef Parallel
                   << "-rank-" << tarch::parallel::Node::getInstance().getRank()
    #endif
                   << "-" << _fileCounter << ".vtk";

    // See issue #47 for discussion whether to quit program on failure:
    // _patchWriter should raise/throw the C++ Exception or return something in case
    // of failure.
    _patchWriter->writeToFile(snapshotFileName.str());

    delete _vertexDataWriter;
    delete _timeStampDataWriter;
    delete _gridWriter;
    delete _patchWriter;

    _vertexDataWriter    = nullptr;
    _patchWriter         = nullptr;
    _timeStampDataWriter = nullptr;
    _gridWriter          = nullptr;
  }
}


exahype::plotters::ADERDG2VTKBinary::~ADERDG2VTKBinary() {
}


void exahype::plotters::ADERDG2VTKBinary::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {
  if (
    tarch::la::allSmaller(_regionOfInterestLeftBottomFront,offsetOfPatch+sizeOfPatch)
    &&
    tarch::la::allGreater(_regionOfInterestRightTopBack,offsetOfPatch)
  ) {
    assertion( _writtenUnknowns==0 || _patchWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _gridWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _timeStampDataWriter!=nullptr );

    int vertexIndex = _writtenUnknowns==0 ? -1 : _gridWriter->plotPatch(offsetOfPatch, sizeOfPatch, _order).first;

    double* interpoland = new double[_solverUnknowns];
    double* value       = _writtenUnknowns==0 ? nullptr : new double[_writtenUnknowns];

    // @todo 16/05/03:Dominic Etienne Charrier
    // This is depending on the choice of basis/implementation.
    // The equidistant grid projection should therefore be moved into the solver.
    dfor(i,_order+1) {
      if (_writtenUnknowns>0) {
        _timeStampDataWriter->plotVertex(vertexIndex, timeStamp);
      }

      for (int unknown=0; unknown < _solverUnknowns; unknown++) {
        interpoland[unknown] = 0.0;
        dfor(ii,_order+1) { // Gauss-Legendre node indices
          int iGauss = peano::utils::dLinearisedWithoutLookup(ii,_order + 1);
          interpoland[unknown] += kernels::equidistantGridProjector1d[_order][ii(1)][i(1)] *
                   kernels::equidistantGridProjector1d[_order][ii(0)][i(0)] *
                   #ifdef Dim3
                   kernels::equidistantGridProjector1d[_order][ii(2)][i(2)] *
                   #endif
                   u[iGauss * _solverUnknowns + unknown];
          assertion3(interpoland[unknown] == interpoland[unknown], offsetOfPatch, sizeOfPatch, iGauss);
        }
      }

      assertion(sizeOfPatch(0)==sizeOfPatch(1));
      _postProcessing->mapQuantities(
        offsetOfPatch,
        sizeOfPatch,
        offsetOfPatch + i.convertScalar<double>()* (sizeOfPatch(0)/(_order)),
        interpoland,
        value,
        timeStamp
      );

      if (_writtenUnknowns>0) {
        _vertexDataWriter->plotVertex(vertexIndex, value, _writtenUnknowns );
      }

      vertexIndex++;
    }

    delete[] interpoland;
    delete[] value;
  }
}
