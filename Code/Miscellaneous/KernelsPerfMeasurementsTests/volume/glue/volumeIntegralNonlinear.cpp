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
 

#include <cstring>
extern double** Kxi;

void volumeIntegralNonlinear(double* lduh, const double* const lFhi,
                             const double* dx,
                             const int numberOfVariables, const int basisSize) {
  const int basisSize2 = basisSize * basisSize;
  const int basisSize3 = basisSize2 * basisSize;

  // Initialize the update DOF
  std::memset(lduh, 0, basisSize3 * numberOfVariables * sizeof(double));

  // x-direction
  const int x_offset = 0 * basisSize3 * numberOfVariables;
  for (int i = 0; i < basisSize; i++) {
    const int i_offset = i * basisSize2 * numberOfVariables;

    for (int j = 0; j < basisSize; j++) {
      const int j_offset = j * basisSize * numberOfVariables;
      const double weight = weights2[i] *
                            weights2[j];
      const double updateSize = weight / dx[0];

      // Fortran: lduh(l, k, j, i) += lFhi_x(l, m, j, i) * Kxi(m, k)
      // Matrix product: (l, m) * (m, k) = (l, k)
      for (int k = 0; k < basisSize; k++) {
        const int k_offset = k * numberOfVariables;

        for (int l = 0; l < numberOfVariables; l++) {
          for (int m = 0; m < basisSize; m++) {
            lduh[i_offset + j_offset + k_offset + l] +=
                Kxi[k][m] *
                lFhi[x_offset + i_offset + j_offset + m * numberOfVariables +
                     l] *
                updateSize;
          }
        }
      }
    }
  }

  // y-direction
  const int y_offset = 1 * basisSize3 * numberOfVariables;
  for (int i = 0; i < basisSize; i++) {
    const int i_offset = i * basisSize2 * numberOfVariables;

    for (int j = 0; j < basisSize; j++) {
      const int j_offset = j * basisSize * numberOfVariables;
      const double weight = weights2[i] *
                            weights2[j];
      const double updateSize = weight / dx[1];

      // Fortran: lduh(l, j, k, i) += lFhi_y(l,m,j,i) * Kxi(m, k)
      // Matrix product: (l, m) * (m, k) = (l, k)
      for (int k = 0; k < basisSize; k++) {
        for (int l = 0; l < numberOfVariables; l++) {
          for (int m = 0; m < basisSize; m++) {
            lduh[i_offset + k * basisSize * numberOfVariables +
                 j * numberOfVariables + l] +=
                Kxi[k][m] *
                lFhi[y_offset + i_offset + j_offset + m * numberOfVariables +
                     l] *
                updateSize;
          }
        }
      }
    }
  }

  // z-direction
  const int z_offset = 2 * basisSize3 * numberOfVariables;
  for (int i = 0; i < basisSize; i++) {
    const int i_offset = i * basisSize2 * numberOfVariables;

    for (int j = 0; j < basisSize; j++) {
      const int j_offset = j * basisSize * numberOfVariables;
      const double weight = weights2[i] *
                            weights2[j];
      const double updateSize = weight / dx[2];

      // Fortran: lduh(l, j, i, k) += lFhi_z(l, m, j, i) * Kxi(m, k)
      // Matrix product (l, m) * (m, k) = (l, k)
      for (int k = 0; k < basisSize; k++) {
        for (int l = 0; l < numberOfVariables; l++) {
          for (int m = 0; m < basisSize; m++) {
            lduh[k * basisSize2 * numberOfVariables +
                 i * basisSize * numberOfVariables + j * numberOfVariables +
                 l] += Kxi[k][m] *
                       lFhi[z_offset + i_offset + j_offset +
                            m * numberOfVariables + l] *
                       updateSize;
          }
        }
      }
    }
  }
}

