#include "tarch/plotter/griddata/unstructured/"

#ifndef EXAHYPE_ADERDG_ADERDGIO_H_
#define EXAHYPE_ADERDG_ADERDGIO_H_

namespace exahype {
  namespace aderdg {
    namespace io {
      /**
       * @brief Exports the DG solution to VTK.
       *
       * @todo: DEC: This is not the final interface yet.
       *             What we want to output is user specific too.
       *
       * @param[in] luh       An array of size basisSize**2 containing the solution DoF.
       * @param[in] center    An array of size dim containing the cell center.
       * @param[in] dx        An array of size dim containing the extent of the cell.
       * @param[in] nvar      The number of conserved variables.
       * @param[in] basisSize The size of the 1D basis.
       */
      template <int dim>
      double exportToVTK(
          tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter*                 _vtkWriter,
          tarch::plotter::griddata::unstructured::UnstructuredGridWriter::VertexWriter*   _vertexWriter,
          tarch::plotter::griddata::unstructured::UnstructuredGridWriter::CellWriter*     _cellWriter,
          tarch::plotter::griddata::Writer::VertexDataWriter*                             _vertexValueWriter,
          const double * const luh,
          const double * const center,
          const double * const dx,
          const int nvar,
          const int basisSize
      );
    } // namespace io
  }  // namespace aderdg
}  // namespace exahype

// 2D specialisations
template <>
double exahype::aderdg::io::exportToVTK<2>(
    tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter*                 _vtkWriter,
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter::VertexWriter*   _vertexWriter,
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter::CellWriter*     _cellWriter,
    tarch::plotter::griddata::Writer::VertexDataWriter*                             _vertexValueWriter,
    const double * const luh,
    const double * const center,
    const double * const dx,
    const int nvar,
    const int basisSize
);

// 3D specialisations
template <>
double exahype::aderdg::io::exportToVTK<3>(
    tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter*                 _vtkWriter,
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter::VertexWriter*   _vertexWriter,
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter::CellWriter*     _cellWriter,
    tarch::plotter::griddata::Writer::VertexDataWriter*                             _vertexValueWriter,
    const double * const luh,
    const double * const center,
    const double * const dx,
    const int nvar,
    const int basisSize
);

// No includes since there are no partial template function specialisations allowed by the standard
// and full specialisations belong in a source file.

#endif /* EXAHYPE_ADERDG_ADERDGIO_H_ */
