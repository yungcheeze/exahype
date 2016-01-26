/*
#include "exahype/aderdg/ADERDG.h"

#include "kernels/quad/GaussLegendre.h"
#include "kernels/geometry/ElementMapping.h"

#include "kernels/initialvalues/PDEInitialValues.h"

// 3D specialisation
template <>
void exahype::aderdg::initialValues<3>(
    double * luh,
    const double * const center,
    const double * const dx,
    const int nvar,
    const int basisSize
) {
  // todo insert your code here
}

// 2D specialisation
template <>
void exahype::aderdg::initialValues<2>(
    double * luh,
    const double * const center,
    const double * const dx,
    const int nvar,
    const int basisSize
) {
  double x,y;

  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    for (int jj=0; jj<basisSize; jj++) {
      // location and index of nodal degrees of freedom
      const int nodeIndex = ii + basisSize * jj;

      // const double * const nodeIndex = { ii, jj }
      // exahype::aderdg::ADERDGInitialValue<2>(luh,nvar,value)

      const double qr = exahype::quad::gaussLegendreNodes[ii];
      const double qs = exahype::quad::gaussLegendreNodes[jj];
      exahype::geometry::mapping2d(center[0],center[1],dx[0],dx[1],qr,qs,&x,&y);

      // read and set the initial values
      const int dofStartIndex  = nodeIndex * nvar;
      exahype::pde::PDEInitialValue2d(x,y,nvar,&luh[dofStartIndex]);
    }
  }
}
*/
