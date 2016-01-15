#include "exahype/aderdg/ADERDG.h"

#include "kernels/quad/GaussLegendre.h"
#include "kernels/geometry/ElementMapping.h"

// 3D specialisation
template <>
void exahype::aderdg::initialValues<3>(
    double * luh,
    double * center,
    const int nvar,
    const int basisSize
) {
  // todo insert your code here
}

// 2D specialisation
template <>
void exahype::aderdg::initialValues<2>(
    double * luh,
    double * center,
    const int nvar,
    const int basisSize
) {
  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    for (int jj=0; jj<basisSize; jj++) {
      // location and index of nodal degrees of freedom
      const int nodeIndex = ii + basisSize * jj;

      // const double * const nodeIndex = { ii, jj }
      // exahype::aderdg::ADERDGInitialValue<2>(luh,nvar,value)

      const double qr = exahype::quad::gaussLegendreNodes[ii];
      const double qs = exahype::quad::gaussLegendreNodes[jj];
      exahype::geometry::mapping2d(center(0),center(1),dx,dy,qr,qs,&x,&y);

      // read initial condition
      exahype::problem::PDEInitialValue2d(x,y,nvar,value);

      // set the DoF
      const int dofStartIndex  = nodeIndex * nvar;

      for (int ivar=0; ivar < nvar; ivar++) {
        luh[dofStartIndex+ivar] = value[ivar];
      }
    }
  }
}
