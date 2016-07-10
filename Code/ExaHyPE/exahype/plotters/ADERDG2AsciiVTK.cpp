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

int exahype::plotters::ADERDG2AsciiVTK::FileCounter(0);

exahype::plotters::ADERDG2AsciiVTK::ADERDG2AsciiVTK() {
}


std::string exahype::plotters::ADERDG2AsciiVTK::getIdentifier() {
  return "vtk::ascii";
}


void exahype::plotters::ADERDG2AsciiVTK::init(
  const std::string& filename,
  int                order,
  int                unknowns,
  const std::string& select) {
  _filename    = filename;
  _order       = order;
  _unknowns    = unknowns;
  _select      = select;
  _patchWriter =
      new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
          new tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter());

  _gridWriter = _patchWriter->createSinglePatchWriter();
  _timeStampDataWriter = _patchWriter->createVertexDataWriter("time", 1);

  for (int i = 0; i < unknowns; i++) {
    std::ostringstream identifier;
    identifier << "Q" << i;
    if ( _select.find(identifier.str())!=std::string::npos || _select=="{all}" ) {
      _vertexDataWriter.push_back(
        _patchWriter->createVertexDataWriter(identifier.str(), 1));
    }
  }
}

exahype::plotters::ADERDG2AsciiVTK::~ADERDG2AsciiVTK() {
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
                   << "-" << FileCounter << ".vtk";

  _patchWriter->writeToFile(snapshotFileName.str());

  FileCounter++;

  for (auto& p : _vertexDataWriter) {
    delete p;
  }
  delete _timeStampDataWriter;
  delete _gridWriter;
  delete _patchWriter;
}

void exahype::plotters::ADERDG2AsciiVTK::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {
  int vertexIndex =
      _gridWriter->plotPatch(offsetOfPatch, sizeOfPatch, _order).first;

// @todo 16/05/03:Dominic Etienne Charrier
// This is the dirty way to perform the projection
// Feel free to modify.
// This is depending on the choice of basis/implementation.
// The equidistant grid projection should therefore be moved into the solver.
#if DIMENSIONS == 2
  for (int j = 0; j < _order + 1; j++) {  // regular grid node indices
    for (int i = 0; i < _order + 1; i++) {
      _timeStampDataWriter->plotVertex(vertexIndex, timeStamp);

      int unknownPlotter = 0;
      for (int unknown=0; unknown < _unknowns; unknown++) {
        std::ostringstream identifier;
        identifier << "Q" << unknown;

        if ( _select.find(identifier.str())!=std::string::npos || _select=="{all}" ) {
          double value = 0;
          for (int jj = 0; jj < _order + 1;
               jj++) {  // Gauss-Legendre node indices
            for (int ii = 0; ii < _order + 1; ii++) {
              int iGauss = jj * (_order + 1) + ii;
              value += kernels::equidistantGridProjector1d[_order][ii][i] *
                       kernels::equidistantGridProjector1d[_order][jj][j] *
                       u[iGauss * _unknowns + unknown];
              assertion3(value == value, offsetOfPatch, sizeOfPatch, iGauss);
            }
          }

          _vertexDataWriter[unknownPlotter]->plotVertex(vertexIndex, value);
          unknownPlotter++;
        }
      }
      vertexIndex++;
    }
  }
#else
  // @todo 16/05/03:Dominic Etienne Charrier: THIS IS UNTESTED CODE!!!
  // This is the dirty way to perform the projection
  // Feel free to modify.
  // This is depending on the choice of basis/implementation.
  // The equidistant grid projection should therefore be moved into the solver.
  for (int k = 0; k < _order + 1; k++) {  // regular grid node indices
    for (int j = 0; j < _order + 1; j++) {
      for (int i = 0; i < _order + 1; i++) {
        _timeStampDataWriter->plotVertex(vertexIndex, timeStamp);

        int unknownPlotter = 0;
        for (int unknown=0; unknown < _unknowns; unknown++) {
          std::ostringstream identifier;
          identifier << "Q" << unknown;

          if ( _select.find(identifier.str())!=std::string::npos || _select=="{all}" ) {
            double value = 0;
            for (int jj = 0; jj < _order + 1;
                 jj++) {  // Gauss-Legendre node indices
              for (int ii = 0; ii < _order + 1; ii++) {
                int iGauss = jj * (_order + 1) + ii;
                value += kernels::equidistantGridProjector1d[_order][ii][i] *
                         kernels::equidistantGridProjector1d[_order][jj][j] *
                         u[iGauss * _unknowns + unknown];
                assertion3(value == value, offsetOfPatch, sizeOfPatch, iGauss);
              }
            }

            _vertexDataWriter[unknownPlotter]->plotVertex(vertexIndex, value);
            unknownPlotter++;
          }
        }
        vertexIndex++;
      }
    }
  }
#endif
}
