#ifndef _EXAHYPE_PLOTTERS_ADERDG_2_BINARY_VTK_H_
#define _EXAHYPE_PLOTTERS_ADERDG_2_BINARY_VTK_H_


#include "exahype/plotters/Plotter.h"


#include "tarch/plotter/griddata/blockstructured/PatchWriterUnstructured.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTKBinaryFileWriter.h"



namespace exahype {
  namespace plotters {
    class ADERDG2BinaryVTK;
  }
}


class exahype::plotters::ADERDG2BinaryVTK: public exahype::plotters::Plotter::Device {
  private:
    static int FileCounter;

    const std::string _filename;
    const int         _order;
    const int         _unknowns;

    tarch::plotter::griddata::blockstructured::PatchWriterUnstructured*           _patchWriter;
    tarch::plotter::griddata::blockstructured::PatchWriter::SinglePatchWriter*    _gridWriter;

    tarch::plotter::griddata::Writer::VertexDataWriter*                           _timeStampDataWriter;
    std::vector< tarch::plotter::griddata::Writer::VertexDataWriter* >            _vertexDataWriter;
  public:
    ADERDG2BinaryVTK(const std::string& filename, int order, int unknowns);
    virtual ~ADERDG2BinaryVTK();

    virtual void plotPatch(
      const tarch::la::Vector<DIMENSIONS,double>&  offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS,double>&  sizeOfPatch,
      double* u,
      double timeStamp
    );
};


#endif
