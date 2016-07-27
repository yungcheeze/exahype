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

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>

#include "exahype/tests/kernels/c/FinitevolumesMusclTest.h"
#include "kernels/finitevolumes/muscl/c/3d/solutionUpdate.cpph"

#if DIMENSIONS == 3

namespace exahype {
namespace tests {
namespace c {

namespace {
// Will be set inside testSolutionUpdate
static double a;  // linear advection coefficient (x direction)
static double b;  // linear advection coefficient (y direction)
static double c;  // linear advection coefficient (z direction)
}  // namespace

void FinitevolumesMusclTest::testFlux(const double *const Q, double **F) {
  double *f = F[0];
  double *g = F[1];
  double *h = F[2];

  f[0] = a * Q[0];
  f[1] = a * Q[1];
  f[2] = a * Q[1];
  f[3] = a * Q[1];
  f[4] = a * Q[1];

  g[0] = b * Q[0];
  g[1] = b * Q[1];
  g[2] = b * Q[1];
  g[3] = b * Q[1];
  g[4] = b * Q[1];

  h[0] = c * Q[0];
  h[1] = c * Q[1];
  h[2] = c * Q[1];
  h[3] = c * Q[1];
  h[4] = c * Q[1];
}

void FinitevolumesMusclTest::testEigenvalues(const double *const Q,
                                             const int normalNonZeroIndex,
                                             double *lambda) {
  switch (normalNonZeroIndex) {
    case 0:
      std::fill(lambda, lambda + 5, a);
      break;
    case 1:
      std::fill(lambda, lambda + 5, b);
      break;
    case 2:
      std::fill(lambda, lambda + 5, c);
      break;
    default:
      assert(false);
  }
}

void FinitevolumesMusclTest::testSolutionUpdate() {
  // linear advection x
  {
    a = 1.23;
    const double dt = 0.234;
    const double cfl = 1.0;
    const double dx_scalar = a * dt / cfl;
    const tarch::la::Vector<DIMENSIONS, double> dx(dx_scalar, dx_scalar,
                                                   dx_scalar);

    const int basisSize = 4;  // 4 points per dimension in cell
    const int basisSize2 = basisSize * basisSize;
    const int basisSize3 = basisSize2 * basisSize;
    const int numberOfVariables = 5;
    double *luh[3 * 3 * 3];  // 9 cells
    for (int i = 0; i < 3 * 3 * 3; i++) {
      luh[i] = new double[basisSize3 * numberOfVariables];
    }

    // initialize (lower half -1.0, upper half +1.0)
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          if (k == 0) {
            // set complete cell to -1.0
            std::fill(
                luh[i * 3 * 3 + j * 3 + k],
                luh[i * 3 * 3 + j * 3 + k] + basisSize3 * numberOfVariables,
                -1.0);
          } else if (k == 2) {
            // set complete cell to +1.0
            std::fill(
                luh[i * 3 * 3 + j * 3 + k],
                luh[i * 3 * 3 + j * 3 + k] + basisSize3 * numberOfVariables,
                +1.0);
          } else {  // k == 1 (middle)
            for (int ii = 0; ii < basisSize; ii++) {
              for (int jj = 0; jj < basisSize; jj++) {
                for (int kk = 0; kk < basisSize; kk++) {
                  double value;
                  if (kk < 1) {
                    value = -1.0;
                  } else {
                    value = +1.0;
                  }
                  for (int l = 0; l < numberOfVariables; l++) {
                    luh[i * 3 * 3 + j * 3 + k]
                       [ii * basisSize2 * numberOfVariables +
                        jj * basisSize * numberOfVariables +
                        kk * numberOfVariables + l] = value;
                  }
                }
              }
            }
          }
        }
      }
    }

    std::cout << std::setprecision(5);
    for (int i = 0; i < basisSize; i++) {
      for (int j = 0; j < basisSize; j++) {
        for (int k = 0; k < basisSize; k++) {
          std::cout << i << " " << j << " " << k << ":";
          for (int l = 0; l < numberOfVariables; l++) {
            std::cout << luh[1 * 3 * 3 + 1 * 3 + 1]
                            [i * basisSize2 * numberOfVariables +
                             j * basisSize * numberOfVariables +
                             k * numberOfVariables + l]
                      << ", ";
          }
          std::cout << std::endl;
        }
      }
    }

    std::cout << "---------------------------------------" << std::endl;

    // do time step
    kernels::finitevolumes::muscl::c::solutionUpdate<testFlux, testEigenvalues>(
        luh, dx, dt, numberOfVariables, basisSize);

    for (int i = 0; i < basisSize; i++) {
      for (int j = 0; j < basisSize; j++) {
        for (int k = 0; k < basisSize; k++) {
          std::cout << i << " " << j << " " << k << ":";
          for (int l = 0; l < numberOfVariables; l++) {
            std::cout << luh[1 * 3 * 3 + 1 * 3 + 1]
                            [i * basisSize2 * numberOfVariables +
                             j * basisSize * numberOfVariables +
                             k * numberOfVariables + l]
                      << ", ";
          }
          std::cout << std::endl;
        }
      }
    }

    // check

    for (int i = 0; i < 3 * 3 * 3; i++) {
      delete[] luh[i];
    }
  }

  // linear advection y
  {

  }

  // linear advection z
  {}

  assert("Test done" && false);
}

}  // namepsace c
}  // namespace tests
}  // namespace exahype

#endif  // DIMENSIONS == 3
