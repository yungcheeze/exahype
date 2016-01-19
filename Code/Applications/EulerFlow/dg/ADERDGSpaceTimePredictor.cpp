#include "EulerFlow/dg/ADERDG.h"

#include "tarch/la/ScalarOperations.h"

#include "EulerFlow/quad/GaussLegendre.h"

#include "EulerFlow/dg/DGMatrices.h"
#include "EulerFlow/dg/Constants.h"
#include "EulerFlow/problem/Problem.h"

#include "stdlib.h"

#include "string.h"

// explicit specialisations
template <>
void exahype::dg::spaceTimePredictor<3>(
    double * restrict /*local*/ lQi,
    double * restrict /*local*/ lFi,
    const double * restrict const /*in*/ luh,
    double * restrict /*out*/ lQhi,
    double * restrict /*out*/ lFhi,
    double * restrict /*out*/ lQhbnd,
    double * restrict /*out*/ lFhbnd,
    double * restrict /*local*/ rhs0,
    double * restrict /*local*/ rhs,
    double * restrict /*local*/ tmp,
    const double * restrict const /*in*/ dx,
    const double /*in*/ dt
) {
  constexpr int dim         = DIMENSIONS;         // 3
  constexpr int dimTimesTwo = (2*DIMENSIONS);     // 6
  constexpr int nvar        = EXAHYPE_NVARS;
  constexpr int basisSize   = (EXAHYPE_ORDER+1);

  // todo insert your code here
}


template <>
void exahype::dg::spaceTimePredictor<2>(
    double * restrict /*local*/ lQi,
    double * restrict /*local*/ lFi,
    const double * restrict const /*in*/ luh,
    double * restrict /*out*/ lQhi,
    double * restrict /*out*/ lFhi,
    double * restrict /*out*/ lQhbnd,
    double * restrict /*out*/ lFhbnd,
    double * restrict /*local*/ rhs0,
    double * restrict /*local*/ rhs,
    double * restrict /*local*/ tmp,
    const double * restrict const /*in*/ dx,
    const double /*in*/ dt
) {
  constexpr int dim         = DIMENSIONS;         // 2
  constexpr int dimTimesTwo = (2*DIMENSIONS);     // 4
  constexpr int basisSize   = (EXAHYPE_ORDER+1);
  constexpr int nvar        = EXAHYPE_NVARS;

  // helper variables
  constexpr int numberOfSpaceTimeDof  = nvar * power(basisSize,dim+1);
  constexpr int numberOfDof           = nvar * power(basisSize, dim);

  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    for (int jj=0; jj<basisSize; jj++) {
      for (int ll=0; ll<basisSize; ll++) { // loop over dof
        // location and index of nodal degrees of freedom
        const int nodeIndex          = ii + basisSize * jj;
        const int spaceTimeNodeIndex = ii + basisSize * jj  + basisSize * basisSize * ll;

        const int dofStartIndex           = nodeIndex * nvar;
        const int spaceTimeDofStartIndex  = spaceTimeNodeIndex * nvar;

        for (int ivar=0; ivar < nvar; ivar++) {
          // Trivial initial guess (can be significantly improved)
          lQi[spaceTimeDofStartIndex+ivar] = luh[dofStartIndex+ivar];

          // Compute the contribution of the initial condition uh to the time update. I prefer to compute it once
          // and store it in rhs0, but if you think it is faster, you can also recompute this contribution
          // inside the Picard loop (DO iter = 1, N+1)
          rhs0[spaceTimeDofStartIndex+ivar] =
              quad::gaussLegendreWeights[ii] *
              quad::gaussLegendreWeights[jj] *
              dg::F0[ll] *
              luh[dofStartIndex+ivar];
        }
      }
    }
  }
  // Above seems to work!

  double*    Q;
  double*    f;
  double*    g;

  //double* dqdt = (double*) std::malloc(nvar * basisSize * sizeof(double)); // todo this is just for debugging; in general, do not use mallocs

  // Discrete Picard iterations. This set of nested loops should (theoretically) be a dream for vectorization, since they are rather independent...
  for (int iter=1; iter < basisSize+1; iter++) {

    // Compute the fluxes (once these fluxes are available, the subsequent operations are independent from each other)
    for (int ll=0; ll<basisSize; ll++) { // loop over dof (time)
      for (int ii=0; ii<basisSize; ii++) { // loop over dof
        for (int jj=0; jj<basisSize; jj++) {
          const int nodeIndex          = ii + basisSize * jj;
          const int spaceTimeNodeIndex = nodeIndex  + basisSize * basisSize * ll;

          const int spaceTimeDofStartIndex     = spaceTimeNodeIndex * nvar;
          const int spaceTimeFluxDofStartIndex = spaceTimeDofStartIndex * dim;

          Q = &lQi [spaceTimeDofStartIndex        ];
          f = &lFi[spaceTimeFluxDofStartIndex     ];
          g = &lFi[spaceTimeFluxDofStartIndex+nvar];
          exahype::problem::PDEFlux(Q,nvar,f,g);
        }
      }
      // Above seems okay!

      // x direction (independent from the y and z derivatives)
      // Kxi : basisSize * basisSize
      // lFh : nvar * basisSize

      // Compute the "derivatives" (contributions of the stiffness matrix)
      // x direction (independent from the y and z derivatives)
      for (int ii=0; ii<basisSize; ii++) { // loop over dof
        for (int jj=0; jj<basisSize; jj++) {
          const int nodeIndex              = ii + basisSize * jj;
          const int spaceTimeNodeIndex     = nodeIndex  + basisSize * basisSize * ll;
          const int spaceTimeDofStartIndex = spaceTimeNodeIndex * nvar;

          double weight = quad::gaussLegendreWeights[ll] *
              quad::gaussLegendreWeights[jj];

// COMPUTE SPATIAL DERIVATIVES FOR TESTING PURPOSES
//for (int ll=0; ll < basisSize; ll++) { // set tmp = 0
//  for(int ivar=0; ivar < nvar; ivar++) {
//    tmp[ivar + nvar*ll] = 0.;
//  }
//}
//
//for(int mm=0; mm < basisSize; mm++) {
//  const int mmNodeIndex        = mm + basisSize * jj;
//  const int mmSpaceTimeNodeIndex         = mmNodeIndex  + basisSize * basisSize * ll;
//  const int mmSpaceTimeDofStartIndex     = mmSpaceTimeNodeIndex * nvar;
//  const int mmSpaceTimeFluxDofStartIndex = mmSpaceTimeDofStartIndex * dim;
//
//  Q = &(lQi [mmSpaceTimeDofStartIndex]);
//  for(int ivar=0; ivar < nvar; ivar++) {
//    tmp[ivar] += 1./dxPatch * dg::dudx[ii][mm] * Q[ivar];
//  }
//}
          for(int ivar=0; ivar < nvar; ivar++) {
            rhs[spaceTimeDofStartIndex+ivar] = rhs0[spaceTimeDofStartIndex+ivar];
          }

          for(int mm=0; mm < basisSize; mm++) {
            const int mmNodeIndex                  = mm + basisSize * jj;
            const int mmSpaceTimeNodeIndex         = mmNodeIndex  + basisSize * basisSize * ll;
            const int mmSpaceTimeDofStartIndex     = mmSpaceTimeNodeIndex * nvar;
            const int mmSpaceTimeFluxDofStartIndex = mmSpaceTimeDofStartIndex * dim;

            f = &lFi[mmSpaceTimeFluxDofStartIndex];

            for(int ivar=0; ivar < nvar; ivar++) {
              rhs[spaceTimeDofStartIndex+ivar]
                  -= weight * dt/dx[0] * dg::Kxi[mm][ii] * f[ivar];
            }
          }
        }
      }
      // Above seems okay!

      // Compute the "derivatives" (contributions of the stiffness matrix)
      // y direction (independent from the x and z derivatives)
      for (int ii=0; ii<basisSize; ii++) { // loop over dof
        for (int jj=0; jj<basisSize; jj++) {
          const int nodeIndex              = ii + basisSize * jj;
          const int spaceTimeNodeIndex     = nodeIndex  + basisSize * basisSize * ll;
          const int spaceTimeDofStartIndex = spaceTimeNodeIndex * nvar;

          double weight = quad::gaussLegendreWeights[ll] *
              quad::gaussLegendreWeights[ii];

          for(int mm=0; mm < basisSize; mm++) {
            const int mmNodeIndex                  = ii + basisSize * mm;
            const int mmSpaceTimeNodeIndex         = mmNodeIndex  + basisSize * basisSize * ll;
            const int mmSpaceTimeDofStartIndex     = mmSpaceTimeNodeIndex * nvar;
            const int mmSpaceTimeFluxDofStartIndex = mmSpaceTimeDofStartIndex * dim;

            g = &lFi[mmSpaceTimeFluxDofStartIndex+nvar];

            for(int ivar=0; ivar < nvar; ivar++) {
              rhs[spaceTimeDofStartIndex+ivar]
                  -= weight * dt/dx[1] * dg::Kxi[mm][jj] * g[ivar];
            }
          }
        }
      }
    } // end of time dof loop

    // Above seems okay!

    for (int ii=0; ii<basisSize; ii++) {  // loop over dof
      for (int jj=0; jj<basisSize; jj++) {
        const int nodeIndex = ii + basisSize * jj;

        double iWeight = 1./(quad::gaussLegendreWeights[ii] * quad::gaussLegendreWeights[jj]);

        for (int ll=0; ll < basisSize; ll++) { // set tmp = 0
          for(int ivar=0; ivar < nvar; ivar++) {
            tmp[ivar + nvar*ll] = 0.;
          }
        }

        for (int ll=0; ll<basisSize; ll++) { // loop over dof

          for(int nn=0; nn < basisSize; nn++) {
            const int nnSpaceTimeNodeIndex     = nodeIndex  + basisSize * basisSize * nn;
            const int nnSpaceTimeDofStartIndex = nnSpaceTimeNodeIndex * nvar;

            for(int ivar=0; ivar < nvar; ivar++) {
              tmp[ivar + nvar*ll] += iWeight * dg::iK1[ll][nn] * rhs[nnSpaceTimeDofStartIndex+ivar];
            }
          }
        }

        for (int ll=0; ll<basisSize; ll++) { // loop over dof
          const int spaceTimeNodeIndex     = nodeIndex  + basisSize * basisSize * ll;
          const int spaceTimeDofStartIndex = spaceTimeNodeIndex * nvar;

          for(int ivar=0; ivar < nvar; ivar++) {
            lQi[spaceTimeDofStartIndex+ivar] = tmp[ivar + nvar*ll];
          }
        }

// UNCOMMENT FOR DEBUGGING PURPOSES
//        // dqdt
//        for (int ll=0; ll < basisSize; ll++) { // set tmp = 0
//          for(int ivar=0; ivar < nvar; ivar++) {
//            dqdt[ivar + nvar*ll] = 0.;
//          }
//        }
//
//        for (int ll=0; ll<basisSize; ll++) { // loop over dof
//
//          for(int ivar=0; ivar < nvar; ivar++) {
//
//            for(int nn=0; nn < basisSize; nn++) {
//              const int nnSpaceTimeNodeIndex         = nodeIndex  + basisSize * basisSize * nn;
//              const int nnSpaceTimeDofStartIndex     = nnSpaceTimeNodeIndex * nvar;
//
//              dqdt[ivar + nvar*ll] += 1./dt * dg::dudx[ll][nn] *
//                  lQi[nnSpaceTimeDofStartIndex+ivar];
//            }
//          }
//        }
      }
    }
  } // end of Picard iteration

  /////////////////////////////////////////////////
  // Post processing of the predictor:
  // Immediately compute the time-averaged space-time polynomials
  /////////////////////////////////////////////////
  constexpr int numberOfFluxDof  = numberOfDof * dim;

  memset((double *) lQhi,0,sizeof(double) * numberOfDof);
  memset((double *) lFhi,0,sizeof(double) * numberOfFluxDof);

  // memory layout of lFhi:
  // lFhi = [ lFhi_x | lFhi_y ] ordered as
  // (a) lFhi_x[nDOF_y][nDOF_x][nVar]
  // (b) lFhi_y[nDOF_y][nDOF_x][nVar]
  //
  double * lFhi_x = &lFhi[0];
  double * lFhi_y = &lFhi[numberOfDof];

  for (int ii=0; ii<basisSize; ii++) { // loop over dof
    for (int jj=0; jj<basisSize; jj++) {
      const int nodeIndex         = ii + basisSize * jj;
      const int dofStartIndex     = nodeIndex * nvar;
      //const int fluxDofStartIndex = dim * dofStartIndex;

      for (int ll=0; ll<basisSize; ll++) { // loop over dof
        const int spaceTimeNodeIndex         = nodeIndex  + basisSize * basisSize * ll;
        const int spaceTimeDofStartIndex     = spaceTimeNodeIndex * nvar;
        const int spaceTimeFluxDofStartIndex = spaceTimeDofStartIndex * dim;

        Q = &lQi[spaceTimeDofStartIndex];

        f = &lFi[spaceTimeFluxDofStartIndex     ];
        g = &lFi[spaceTimeFluxDofStartIndex+nvar];

        double weight = quad::gaussLegendreWeights[ll];

        //double * temp = &(lQhi[dofStartIndex]);
        for(int ivar=0; ivar < nvar; ivar++) {
          lQhi[dofStartIndex+ivar] += weight * Q[ivar];

          lFhi_x[dofStartIndex+ivar] += weight * f[ivar];
          lFhi_y[dofStartIndex+ivar] += weight * g[ivar];
        }
      }
    }
  }

  // TODO
  // Ideas for improvement
  // (a) reorder lFhi_y, lFhi_z do have better access mechanism
  //     in SurfaceIntegral??? e.g. lFhi_y(nvar, i,j,k) becomes lFhi_y(nvar,j,i,k)
  //     and then the usual access mechanism


  /////////////////////////////////////////////////
  // Compute the bounday-extrapolated values for Q and F*n
  /////////////////////////////////////////////////
  constexpr int numberOfFaceDof = nvar * power(basisSize, dim-1);

  memset((double *) &lQhbnd[0],0,sizeof(double) * numberOfFaceDof * dimTimesTwo);
  memset((double *) &lFhbnd[0],0,sizeof(double) * numberOfFaceDof * dimTimesTwo);

  //
  // x-direction: face 0 (left) and face 1 (right)
  //
  for (int jj=0; jj<basisSize; jj++) {
    const int nodeIndex      = jj;
    const int dofStartIndexL = EXAHYPE_FACE_LEFT  * numberOfFaceDof + nodeIndex * nvar;
    const int dofStartIndexR = EXAHYPE_FACE_RIGHT * numberOfFaceDof + nodeIndex * nvar;

    for (int mm=0; mm<basisSize; mm++) { // loop over dof
      const int mmNodeIndex         = mm  + basisSize * jj;
      const int mmDofStartIndex     = mmNodeIndex * nvar;

      Q = &lQhi[mmDofStartIndex];

      for(int ivar=0; ivar < nvar; ivar++) {
        lQhbnd[dofStartIndexL+ivar] += dg::FLCoeff[mm] * Q[ivar];
        lQhbnd[dofStartIndexR+ivar] += dg::FRCoeff[mm] * Q[ivar];

        lFhbnd[dofStartIndexL+ivar] += dg::FLCoeff[mm] * lFhi_x[mmDofStartIndex+ivar]; // lFhi_x * FLCoeff
        lFhbnd[dofStartIndexR+ivar] += dg::FRCoeff[mm] * lFhi_x[mmDofStartIndex+ivar]; // lFhi_x * FRCoeff
      }
    }
    continue;
  }

  //
  // y-direction: face 2 (left) and face 3 (right)
  //
  for (int ii=0; ii<basisSize; ii++) {
    const int nodeIndex      = ii;
    const int dofStartIndexL = EXAHYPE_FACE_FRONT * numberOfFaceDof + nodeIndex * nvar;
    const int dofStartIndexR = EXAHYPE_FACE_BACK  * numberOfFaceDof + nodeIndex * nvar;

    for (int mm=0; mm<basisSize; mm++) {
      const int mmNodeIndex         = ii  + basisSize * mm;
      const int mmDofStartIndex     = mmNodeIndex * nvar;

      Q = &lQhi [mmDofStartIndex         ];

      for(int ivar=0; ivar < nvar; ivar++) {
        lQhbnd[dofStartIndexL+ivar] += dg::FLCoeff[mm] * Q[ivar];
        lQhbnd[dofStartIndexR+ivar] += dg::FRCoeff[mm] * Q[ivar];

        lFhbnd[dofStartIndexL+ivar] += dg::FLCoeff[mm] * lFhi_y[mmDofStartIndex+ivar]; // lFhi_y * FLCoeff
        lFhbnd[dofStartIndexR+ivar] += dg::FRCoeff[mm] * lFhi_y[mmDofStartIndex+ivar]; // lFhi_y * FRCoeff
      }
    }
    continue;
  }

  // clean up
  //  std::free(dqdt);
}
