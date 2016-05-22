#include "kernels/aderdg/generic/Kernels.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

#include "string.h"

#if DIMENSIONS == 2
void singleLevelFaceUnknownsProlongation(
    double* lQhbndFine,
    const double* lQhbndCoarse,
    const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize) {
  for (int m1 = 0; m1 < basisSize; ++m1) {
    for (int ivar = 0; ivar < numberOfVariables; ++ivar) {
      const int mNodeIndex     = m1;
      const int mDofStartIndex = mNodeIndex * numberOfVariables;

      for (int n1 = 0; n1 < basisSize; ++n1) {
        const int nNodeIndex = n1;
        const int nDofStartIndex = nNodeIndex * numberOfVariables;

        lQhbndFine[mDofStartIndex+ivar] +=
            lQhbndCoarse[nDofStartIndex + ivar] *
            kernels::fineGridProjector1d[basisSize-1][subfaceIndex[0]][n1][m1]; // q_n * Pi_nm
      }
    }
  }
}

void kernels::aderdg::generic::c::faceUnknownsProlongation(double* lQhbndFine,
                              double* lFhbndFine,
                              const double* lQhbndCoarse,
                              const double* lFhbndCoarse,
                              const int coarseGridLevel,
                              const int fineGridLevel,
                              const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex,
                              const int numberOfVariables, const int basisSize){
  const int levelDelta = fineGridLevel - coarseGridLevel;

  double * lQhbndFineTemp = new double[basisSize*numberOfVariables];
  double * lFhbndFineTemp = new double[basisSize*numberOfVariables];

  double * workPointerQhbnd1 = 0;
  double * workPointerFhbnd1 = 0;

  double * workPointerQhbnd2 = 0;
  double * workPointerFhbnd2 = 0;

  // This ensures that the last workPointerQhbnd1 points to lQhbndFine.
  // Analogous is done for workPointerFhbnd1.
  if (levelDelta % 2 == 0) {
    workPointerQhbnd1 = lQhbndFineTemp;
    workPointerFhbnd1 = lFhbndFineTemp;
  } else {
    workPointerQhbnd1 = lQhbndFine;
    workPointerFhbnd1 = lFhbndFine;
  }

  tarch::la::Vector<DIMENSIONS-1,int> subfaceIndexPrevious (subfaceIndex);
  tarch::la::Vector<DIMENSIONS-1,int> subfaceIndexCurrent;

  tarch::la::Vector<DIMENSIONS-1,int> subintervalIndex;
  // This loop step by step decodes subfaceIndex[0] into a tertiary basis
  // starting with the highest significance 3^(levelDelta-1).
  // Per step, the digit corresponding to the current significance then determines
  // the subinterval for the one level prolongation (which is performed
  // in the same step).
  for (int l = 1; l < levelDelta+1; ++l) {
    const int significance = tarch::la::aPowI(levelDelta-l,3);
    subfaceIndexCurrent[0] = subfaceIndexPrevious[0] % significance;
    subintervalIndex[0]    = (subfaceIndexPrevious[0] - subfaceIndexCurrent[0])/significance;
    assertion(subintervalIndex[0] < 3);

    // Zero the values of 'workPointer1' since we compute a sum for each entry.
    memset(workPointerQhbnd1, 0, sizeof(double)*numberOfVariables*basisSize);
    memset(workPointerFhbnd1, 0, sizeof(double)*numberOfVariables*basisSize);

    // Apply the one level prolongation operator.
    // Use the coarse grid unknowns in the first iteration
    // and use workPointerQhbnd2,workPointerFhbnd2 in the following.
    if (l==1) {
      singleLevelFaceUnknownsProlongation(
          workPointerQhbnd1,
          lQhbndCoarse,
          subintervalIndex,
          numberOfVariables,
          basisSize);

      singleLevelFaceUnknownsProlongation(
          workPointerFhbnd1,
          lFhbndCoarse,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    } else {
      singleLevelFaceUnknownsProlongation(
          workPointerQhbnd1,
          workPointerQhbnd2,
          subintervalIndex,
          numberOfVariables,
          basisSize);

      singleLevelFaceUnknownsProlongation(
          workPointerFhbnd1,
          workPointerFhbnd2,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    }

    // Prepare next iteration.
    subfaceIndexPrevious = subfaceIndexCurrent;

    workPointerQhbnd2 = workPointerQhbnd1;
    workPointerFhbnd2 = workPointerFhbnd1;
    // One check is sufficient since we change both work pointers
    // for lQhbnd and lFhbnd at the same time.
    if (workPointerQhbnd1 == lQhbndFineTemp) {
      workPointerQhbnd1 = lQhbndFine;
      workPointerFhbnd1 = lFhbndFine;
    } else {
      workPointerQhbnd1 = lQhbndFineTemp;
      workPointerFhbnd1 = lFhbndFineTemp;
    }
  }

  // Clean up.
  delete [] lQhbndFineTemp;
  delete [] lFhbndFineTemp;
}

void singleLevelFaceUnknownsRestriction(
    double* lQhbndCoarse,
    const double* lQhbndFine,
    const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize) {
  for (int m1 = 0; m1 < basisSize; ++m1) {
    for (int ivar = 0; ivar < numberOfVariables; ++ivar) {
      const int mNodeIndex     = m1;
      const int mDofStartIndex = mNodeIndex * numberOfVariables;

      for (int n1 = 0; n1 < basisSize; ++n1) {
        const int nNodeIndex = n1;
        const int nDofStartIndex = nNodeIndex * numberOfVariables;

        // Pi_mn * q_n
        // Note that the indices are interchanged
        // in expression 'fineGridProjector1d[order][digit][m][n]'
        // in comparsion to the one level prolongation operator.
        lQhbndCoarse[mDofStartIndex+ivar] +=
            kernels::gaussLegendreWeights[basisSize-1][n1] *
            kernels::fineGridProjector1d[basisSize-1][subfaceIndex[0]][m1][n1] *
            lQhbndFine[nDofStartIndex + ivar] /
            kernels::gaussLegendreWeights[basisSize-1][m1] / 3.0;
      }
    }
  }
}

void accumulate(
    double * inOutArray,
    const double * inArray,
    const int length) {
  for (int i = 0; i < length; ++i) {
    inOutArray[i] += inArray[i];
  }
}

void kernels::aderdg::generic::c::faceUnknownsRestriction(double* lQhbndCoarse,
                             double* lFhbndCoarse,
                             const double* lQhbndFine,
                             const double* lFhbndFine,
                             const int coarseGridLevel,
                             const int fineGridLevel,
                             const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex,
                             const int numberOfVariables, const int basisSize){
  const int levelDelta     = fineGridLevel - coarseGridLevel;

  double * lQhbndCoarseTemp1 = new double[basisSize*numberOfVariables];
  double * lFhbndCoarseTemp1 = new double[basisSize*numberOfVariables];
  double * lQhbndCoarseTemp2 = new double[basisSize*numberOfVariables];
  double * lFhbndCoarseTemp2 = new double[basisSize*numberOfVariables];

  double * workPointerQhbnd1 = 0;
  double * workPointerQhbnd2 = 0;
  double * workPointerFhbnd1 = 0;
  double * workPointerFhbnd2 = 0;

  // This ensures that the last workPointerQhbnd1 points to lQhbndFine.
  // Analogous is done for workPointerFhbnd1.
  workPointerQhbnd1 = lQhbndCoarseTemp1;
  workPointerFhbnd1 = lFhbndCoarseTemp1;

  tarch::la::Vector<DIMENSIONS-1,int> subfaceIndexCurrent(subfaceIndex);

  tarch::la::Vector<DIMENSIONS-1,int> subintervalIndex;
  // This loop step by step decodes subfaceIndex[0] into a tertiary basis
  // starting with the lowest significance 3^0 (in contrast to the the prolongation).
  // Per step, the digit corresponding to the current significance then determines
  // the subinterval for the one level prolongation (which is performed
  // in the same step).
  for (int l = 1; l < levelDelta+1; ++l) {
    subintervalIndex[0]    = subfaceIndexCurrent[0] % 3;
    subfaceIndexCurrent[0] = (subfaceIndexCurrent[0] - subintervalIndex[0])/3;
    assertion(subintervalIndex[0] < 3);

    // Zero the values of workPointerQhbnd1 and workPointerFhbnd1 since we compute a sum for each entry.
    memset(workPointerQhbnd1, 0, sizeof(double)*numberOfVariables*basisSize);
    memset(workPointerFhbnd1, 0, sizeof(double)*numberOfVariables*basisSize);

    // Apply the one level restriction operator.
    if (l==1) {
      singleLevelFaceUnknownsRestriction(
          workPointerQhbnd1,
          lQhbndFine,
          subintervalIndex,
          numberOfVariables,
          basisSize);

      singleLevelFaceUnknownsRestriction(
          workPointerFhbnd1,
          lFhbndFine,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    } else {
      singleLevelFaceUnknownsRestriction(
          workPointerQhbnd1,
          workPointerQhbnd2,
          subintervalIndex,
          numberOfVariables,
          basisSize);

      singleLevelFaceUnknownsRestriction(
          workPointerFhbnd1,
          workPointerFhbnd2,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    }

    // Prepare next iteration.
    workPointerQhbnd2 = workPointerQhbnd1;
    workPointerFhbnd2 = workPointerFhbnd1;
    // One check is sufficient since we change both work pointers
    // for lQhbnd and lFhbnd at the same time.
    if (workPointerQhbnd1 == lQhbndCoarseTemp1) {
      workPointerQhbnd1 = lQhbndCoarseTemp2;
      workPointerFhbnd1 = lFhbndCoarseTemp2;
    } else {
      workPointerQhbnd1 = lQhbndCoarseTemp1;
      workPointerFhbnd1 = lFhbndCoarseTemp1;
    }
  }

  // Add to coarse grid degrees of freedom:
  accumulate(lQhbndCoarse,workPointerQhbnd2,numberOfVariables*basisSize);
  accumulate(lFhbndCoarse,workPointerFhbnd2,numberOfVariables*basisSize);

  // Clean up.
  delete [] lQhbndCoarseTemp1;
  delete [] lFhbndCoarseTemp1;
  delete [] lQhbndCoarseTemp2;
  delete [] lFhbndCoarseTemp2;
}

void singleLevelVolumeUnknownsProlongation(
    double* luhFine,
    const double* luhCoarse,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize) {
  for (int m2 = 0; m2 < basisSize; ++m2) {
    for (int m1 = 0; m1 < basisSize; ++m1) {
      for (int ivar = 0; ivar < numberOfVariables; ivar++) {
        const int mNodeIndex     = basisSize*m2 + m1;
        const int mDofStartIndex = mNodeIndex * numberOfVariables;

        for (int n2 = 0; n2 < basisSize; ++n2) {
          for (int n1 = 0; n1 < basisSize; ++n1) {
            const int nNodeIndex     = basisSize*n2 + n1;
            const int nDofStartIndex = nNodeIndex * numberOfVariables;

            luhFine[mDofStartIndex+ivar] +=
                luhCoarse[nDofStartIndex + ivar] *
                kernels::fineGridProjector1d[basisSize-1][subcellIndex[0]][n1][m1] *
                kernels::fineGridProjector1d[basisSize-1][subcellIndex[1]][n2][m2]; // \sum q_n * P_nm
          }
        }
      }
    }
  }
}

void kernels::aderdg::generic::c::volumeUnknownsProlongation(double* luhFine,
                                const double* luhCoarse,
                                const int coarseGridLevel,
                                const int fineGridLevel,
                                const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
                                const int numberOfVariables, const int basisSize){
  const int levelDelta = fineGridLevel - coarseGridLevel;

  double * luhFineTemp = new double[basisSize*basisSize*numberOfVariables];

  double * workPointerUh1 = 0;
  double * workPointerUh2 = 0;

  // This ensures that the last workPointeruh1 points to luhFine.
  // Analogous is done for workPointerFhbnd1.
  if (levelDelta % 2 == 0) {
    workPointerUh1 = luhFineTemp;
  } else {
    workPointerUh1 = luhFine;
  }

  tarch::la::Vector<DIMENSIONS,int> subcellIndexPrevious (subcellIndex);
  tarch::la::Vector<DIMENSIONS,int> subcellIndexCurrent;

  tarch::la::Vector<DIMENSIONS,int> subintervalIndex;
  // This loop step by step decodes subfaceIndex[0] into a tertiary basis
  // starting with the highest significance 3^(levelDelta-1).
  // Per step, the digit corresponding to the current significance then determines
  // the subinterval for the one level prolongation (which is performed
  // in the same step).
  for (int l = 1; l < levelDelta+1; ++l) {
    const int significance = tarch::la::aPowI(levelDelta-l,3);
    subcellIndexCurrent[0] = subcellIndexPrevious[0] % significance;
    subcellIndexCurrent[1] = subcellIndexPrevious[1] % significance;
    subintervalIndex[0]    = (subcellIndexPrevious[0] - subcellIndexCurrent[0])/significance;
    subintervalIndex[1]    = (subcellIndexPrevious[1] - subcellIndexCurrent[1])/significance;
    assertion(subintervalIndex[0] < 3);
    assertion(subintervalIndex[1] < 3);

    // Zero the values of 'workPointer1' since we compute a sum for each entry.
    memset(workPointerUh1, 0, basisSize*basisSize*numberOfVariables*sizeof(double));

    // Apply the one level prolongation operator.
    // Use the coarse grid unknowns in the first iteration
    // and use workPointeruh2,workPointerFhbnd2 in the following.
    if (l==1) {
      singleLevelVolumeUnknownsProlongation(
          workPointerUh1,
          luhCoarse,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    } else {
      singleLevelVolumeUnknownsProlongation(
          workPointerUh1,
          workPointerUh2,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    }

    // Prepare next iteration.
    subcellIndexPrevious = subcellIndexCurrent;

    workPointerUh2 = workPointerUh1;
    // One check is sufficient since we change both work pointers
    // for luh and lFhbnd at the same time.
    if (workPointerUh1 == luhFineTemp) {
      workPointerUh1 = luhFine;
    } else {
      workPointerUh1 = luhFineTemp;
    }
  }

  // Clean up.
  delete [] luhFineTemp;
}

void singleLevelVolumeUnknownsRestriction(
    double* luhCoarse,
    const double* luhFine,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize) {
  for (int m2 = 0; m2 < basisSize; ++m2) {
    for (int m1 = 0; m1 < basisSize; ++m1) {
      for (int ivar = 0; ivar < numberOfVariables; ivar++) {
        const int mNodeIndex     = basisSize*m2 + m1;
        const int mDofStartIndex = mNodeIndex * numberOfVariables;

        for (int n2 = 0; n2 < basisSize; ++n2) {
          for (int n1 = 0; n1 < basisSize; ++n1) {
            const int nNodeIndex     = basisSize*n2 + n1;
            const int nDofStartIndex = nNodeIndex * numberOfVariables;

            // P_mn * q_n
            // Note that the indices are interchanged
            // in expression 'fineGridProjector1d[order][digit][m][n]'
            // in comparsion to the one level prolongation operator.
            luhCoarse[mDofStartIndex+ivar] +=
                kernels::gaussLegendreWeights[basisSize-1][n1] *
                kernels::gaussLegendreWeights[basisSize-1][n2] *
                kernels::fineGridProjector1d[basisSize-1][subcellIndex[0]][m1][n1] *
                kernels::fineGridProjector1d[basisSize-1][subcellIndex[1]][m2][n2] *
                luhFine[nDofStartIndex + ivar] /
                kernels::gaussLegendreWeights[basisSize-1][m1] /
                kernels::gaussLegendreWeights[basisSize-1][m2] / 9.0; // 9=3^DIMENSIONS
          }
        }
      }
    }
  }
}

void kernels::aderdg::generic::c::volumeUnknownsRestriction(double* luhCoarse,
                               const double* luhFine,
                               const int coarseGridLevel,
                               const int fineGridLevel,
                               const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
                               const int numberOfVariables, const int basisSize){
  const int levelDelta     = fineGridLevel - coarseGridLevel;

  double * luhCoarseTemp1 = new double[basisSize*basisSize*numberOfVariables];
  double * luhCoarseTemp2 = new double[basisSize*basisSize*numberOfVariables];

  double * workPointerUh1 = 0;
  double * workPointerUh2 = 0;

  // This ensures that the last workPointeruh1 points to luhFine.
  // Analogous is done for workPointerFhbnd1.
  workPointerUh1 = luhCoarseTemp1;

  tarch::la::Vector<DIMENSIONS,int> subcellIndexCurrent(subcellIndex);

  tarch::la::Vector<DIMENSIONS,int> subintervalIndex;
  // This loop step by step decodes subcellIndex[0] into a tertiary basis
  // starting with the lowest significance 3^0 (in contrast to the prolongation decoding).
  // Per step, the digit corresponding to the current significance then determines
  // the subinterval for the one level prolongation (which is performed
  // in the same step).
  for (int l = 1; l < levelDelta+1; ++l) {
    subintervalIndex[0]    = subcellIndexCurrent[0] % 3;
    subintervalIndex[1]    = subcellIndexCurrent[1] % 3;
    subcellIndexCurrent[0] = (subcellIndexCurrent[0] - subintervalIndex[0])/3;
    subcellIndexCurrent[1] = (subcellIndexCurrent[1] - subintervalIndex[1])/3;
    assertion(subintervalIndex[0] < 3);
    assertion(subintervalIndex[1] < 3);

    // Zero the values of workPointeruh1 since we compute a sum for each entry.
    memset(workPointerUh1, 0, basisSize*basisSize*numberOfVariables*sizeof(double));

    // Apply the one level restriction operator.
    if (l==1) {
      singleLevelVolumeUnknownsRestriction(
          workPointerUh1,
          luhFine,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    } else {
      singleLevelVolumeUnknownsRestriction(
          workPointerUh1,
          workPointerUh2,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    }

    // Prepare next iteration.
    workPointerUh2 = workPointerUh1;
    if (workPointerUh1 == luhCoarseTemp1) {
      workPointerUh1 = luhCoarseTemp2;
    } else {
      workPointerUh1 = luhCoarseTemp1;
    }
  }

  // Add to coarse grid degrees of freedom:
  accumulate(luhCoarse,workPointerUh2,basisSize*basisSize*numberOfVariables);

  // Clean up.
  delete [] luhCoarseTemp1;
  delete [] luhCoarseTemp2;
}
#endif
