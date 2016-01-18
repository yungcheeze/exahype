#include "exahype/aderdg/ADERDGIO.h"

#include "string.h"

#include "kernels/geometry/ElementMapping.h"
#include "kernels/quad/GaussLegendre.h"

#include "kernels/aderdg/DGMatrices.h"

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
) {

}

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
) {
  // basis

  int    indexMapping      [EXAHYPE_ORDER+1][EXAHYPE_ORDER+1];
  double uniformCoordinates[DIMENSIONS     ][EXAHYPE_NBASIS_POWER_DIMENSIONS];
  double uniformPartition  [TWO_POWER_D    ][EXAHYPE_NBASIS_POWER_DIMENSIONS];
  double uniformDoF        [EXAHYPE_NVARS][EXAHYPE_NBASIS_POWER_DIMENSIONS];

  // define sub nodes
  for (int ii=0; ii < basisSize; ii++) {
    for (int jj=0; jj < basisSize; jj++) {
      const int uniformNodeIndex = ii + basisSize * jj;

      indexMapping[ii][jj] = uniformNodeIndex;
      uniformCoordinates[0][uniformNodeIndex] = (double) ii / (double) EXAHYPE_ORDER;
      uniformCoordinates[1][uniformNodeIndex] = (double) jj / (double) EXAHYPE_ORDER;
    }
  }

  // define sub quadrangles/hexahedrons
  for (int ii=0; ii < basisSize; ii++) {
    for (int jj=0; jj < basisSize; jj++) {
      const int uniformNodeIndex = ii + basisSize * jj;

      uniformPartition[0][uniformNodeIndex] = indexMapping[ii  ][jj  ];
      uniformPartition[1][uniformNodeIndex] = indexMapping[ii+1][jj  ];
      uniformPartition[2][uniformNodeIndex] = indexMapping[ii+1][jj+1];
      uniformPartition[3][uniformNodeIndex] = indexMapping[ii  ][jj+1];
    }
  }
  // Map Gauss-Legendre nodes to equidistant subgrid coordinates
  //          std::memset((double *) &subData[0],0,sizeof(double) * EXAHYPE_NVARS * EXAHYPE_NBASIS_POWER_DIMENSIONS);
  for (int ii=0; ii<basisSize; ii++) { // mem zero
    for (int jj=0; jj<basisSize; jj++) {
      const int uniformNodeIndex = ii + basisSize * jj;

      for (int ivar=0; ivar < nvar; ivar++) {
        uniformDoF[ivar][uniformNodeIndex] = 0;
      }
    }
  }

  for (int ii=0; ii<basisSize; ii++) { // project on subgrid coordinates
    for (int jj=0; jj<basisSize; jj++) {
      const int uniformNodeIndex = ii + basisSize * jj;

      for (int mm=0; mm<basisSize; mm++) { // project on subgrid coordinates
        for (int nn=0; nn<basisSize; nn++) {
          const int nodeIndex     = mm + basisSize * nn;
          const int dofStartIndex = nodeIndex * nvar;

          for (int ivar=0; ivar < nvar; ivar++) {
            uniformDoF[ivar][uniformNodeIndex] += luh[dofStartIndex+ivar] * aderdg::subOutputMatrix[nodeIndex][uniformNodeIndex];
          }
        }
      }
    }
  }

  // helpers
  double x,y;

  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    for (int jj=0; jj<basisSize; jj++) {
      // location and index of nodal degrees of freedom
      const int nodeIndex     = ii + basisSize * jj;
      const int dofStartIndex = nodeIndex * nvar;

      const double r = uniformCoordinates[0][nodeIndex];
      const double s = uniformCoordinates[1][nodeIndex];

      // ADERDG VTK export kernel
      geometry::mapping2d(center[0],center[1],dx[0],dx[0],r,s,&x,&y);
      tarch::la::Vector<DIMENSIONS,double> currentVertexPosition(x,y);

      // input basisSize,nvar,vertexWriter*,cellWriter*,VertexValueWriter*

      const int vtkNodeIndex = _vertexWriter->plotVertex(currentVertexPosition);
      _cellWriter->plotPoint(vtkNodeIndex);

      for (int ivar=4; ivar < 5; ivar++) {
        //                const double dofValue = u[dofStartIndex+ivar];
        const double dofValue = uniformDoF[ivar][nodeIndex];
        _vertexValueWriter->plotVertex(vtkNodeIndex,dofValue);
      }
    }
  }
}
