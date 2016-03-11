#include "exahype/plotters/ADERDG2AsciiVTK.h"
#include "tarch/parallel/Node.h"

int exahype::plotters::ADERDG2AsciiVTK::FileCounter(0);

exahype::plotters::ADERDG2AsciiVTK::ADERDG2AsciiVTK(const std::string& filename,
                                                    int order, int unknowns)
    : _filename(filename), _order(order), _unknowns(unknowns) {
  _patchWriter =
      new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
          new tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter());

  _gridWriter = _patchWriter->createSinglePatchWriter();
  _timeStampDataWriter = _patchWriter->createVertexDataWriter("time", 1);

  /// @todo Remove
  assertionEquals(order, 3);

  for (int i = 0; i < unknowns; i++) {
    std::ostringstream identifier;
    identifier << "Q" << i;
    _vertexDataWriter.push_back(
        _patchWriter->createVertexDataWriter(identifier.str(), 1));
  }
}

exahype::plotters::ADERDG2AsciiVTK::~ADERDG2AsciiVTK() {
  _gridWriter->close();
  _timeStampDataWriter->close();
  for (
      std::vector<tarch::plotter::griddata::Writer::VertexDataWriter*>::iterator
          p = _vertexDataWriter.begin();
      p != _vertexDataWriter.end(); p++) {
    (*p)->close();
  }

  std::ostringstream snapshotFileName;
  snapshotFileName << _filename
#ifdef Parallel
                   << "-rank-" << tarch::parallel::Node::getInstance().getRank()
#endif
                   << "-" << FileCounter << ".vtk";

  _patchWriter->writeToFile(snapshotFileName.str());

  FileCounter++;

  for (
      std::vector<tarch::plotter::griddata::Writer::VertexDataWriter*>::iterator
          p = _vertexDataWriter.begin();
      p != _vertexDataWriter.end(); p++) {
    delete *p;
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

  for (int i = 0; i < tarch::la::aPowI(DIMENSIONS, _order + 1); i++) {
    _timeStampDataWriter->plotVertex(vertexIndex, timeStamp);
    int unknown = 0;
    for (std::vector<
             tarch::plotter::griddata::Writer::VertexDataWriter*>::iterator p =
             _vertexDataWriter.begin();
         p != _vertexDataWriter.end(); p++) {
      (*p)->plotVertex(vertexIndex, u[i * _unknowns + unknown]);
      unknown++;
    }

    vertexIndex++;
  }
}
