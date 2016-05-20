#include "kernels/aderdg/generic/Kernels.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

#include "string.h"

#if DIMENSIONS == 2
void singleLevelFaceUnknownsProlongation(
    double* lQhbndFine,
    double* lFhbndFine,
    const double* lQhbndCoarse,
    const double* lFhbndCoarse,
    const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize) {
  for (int mm = 0; mm < basisSize; mm++) {
    for (int ivar = 0; ivar < numberOfVariables; ivar++) {
      const int mmNodeIndex     = mm;
      const int mmDofStartIndex = mmNodeIndex * numberOfVariables;

      for (int ii = 0; ii < basisSize; ii++) {
        const int iiNodeIndex = ii;
        const int iiDofStartIndex = iiNodeIndex * numberOfVariables;

        lQhbndFine[mmDofStartIndex+ivar]
                          += lQhbndCoarse[iiDofStartIndex + ivar] *
                          kernels::fineGridProjector1d[basisSize-1][subfaceIndex[0]][ii][mm]; // q_i * Pi_im

        lFhbndFine[mmDofStartIndex+ivar]
                          += lFhbndCoarse[iiDofStartIndex + ivar] *
                          kernels::fineGridProjector1d[basisSize-1][subfaceIndex[0]][ii][mm]; // f_i * Pi_im
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
          workPointerQhbnd1,workPointerFhbnd1,
          lQhbndCoarse,lFhbndCoarse,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    } else {
      singleLevelFaceUnknownsProlongation(
          workPointerQhbnd1,workPointerFhbnd1,
          workPointerQhbnd2,workPointerFhbnd2,
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

//  memset(lQhbndFine, 0, sizeof(double)*numberOfVariables*basisSize);
//  memset(lFhbndFine, 0, sizeof(double)*numberOfVariables*basisSize);

  assertionNumericalEquals(lQhbndFine[3],0);
//  assertionNumericalEquals(lFhbndFine[0],1);
}

void singleLevelFaceUnknownsRestriction(
    double* lQhbndCoarse,
    double* lFhbndCoarse,
    const double* lQhbndFine,
    const double* lFhbndFine,
    const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize) {
  for (int mm = 0; mm < basisSize; mm++) {
    for (int ivar = 0; ivar < numberOfVariables; ivar++) {
      const int mmNodeIndex     = mm;
      const int mmDofStartIndex = mmNodeIndex * numberOfVariables;

      for (int ii = 0; ii < basisSize; ii++) {
        const int iiNodeIndex = ii;
        const int iiDofStartIndex = iiNodeIndex * numberOfVariables;
        lQhbndCoarse[mmDofStartIndex+ivar]
                          += kernels::gaussLegendreWeights[basisSize-1][ii] *
                             kernels::fineGridProjector1d[basisSize-1][subfaceIndex[0]][mm][ii] *
                             lQhbndFine[iiDofStartIndex + ivar] /
                             kernels::gaussLegendreWeights[basisSize-1][mm] / 3; // Pi_mi * q_i

        lFhbndCoarse[mmDofStartIndex+ivar]
                          += kernels::gaussLegendreWeights[basisSize-1][ii] *
                             kernels::fineGridProjector1d[basisSize-1][subfaceIndex[0]][mm][ii] *
                             lFhbndFine[iiDofStartIndex + ivar] /
                             kernels::gaussLegendreWeights[basisSize-1][mm] / 3; // Pi_mi * f_i
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
    // Note that the indices are interchanged
    // in expression 'fineGridProjector1d[order][digit][mm][ii]'
    // in comparsion to the one level prolongation operator.
    if (l==1) {
      singleLevelFaceUnknownsRestriction(
          workPointerQhbnd1,workPointerFhbnd1,
          lQhbndFine,lFhbndFine,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    } else {
      singleLevelFaceUnknownsRestriction(
          workPointerQhbnd1,workPointerFhbnd1,
          workPointerQhbnd2,workPointerFhbnd2,
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
//  memset(workPointerQhbnd2, 0, sizeof(double)*numberOfVariables*basisSize);
//  memset(workPointerFhbnd2, 0, sizeof(double)*numberOfVariables*basisSize);

  accumulate(lQhbndCoarse,workPointerQhbnd2,numberOfVariables*basisSize);
  accumulate(lFhbndCoarse,workPointerFhbnd2,numberOfVariables*basisSize);

//  assertionNumericalEquals(lQhbndCoarse[0],1);
//  assertionNumericalEquals(workPointerQhbnd2[0]*gaussLegendreWeights[basisSize-1][1],1./3.);
//  assertionNumericalEquals(lFhbndCoarse[0],0);

  // Clean up.
  delete [] lQhbndCoarseTemp1;
  delete [] lFhbndCoarseTemp1;
  delete [] lQhbndCoarseTemp2;
  delete [] lFhbndCoarseTemp2;
}

void kernels::aderdg::generic::c::volumeUnknownsProlongation(double* luhFine,
                                const double* luhCoarse,
                                const int coarseGridLevel,
                                const int fineGridLevel,
                                const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
                                const int numberOfVariables, const int basisSize){
  // @todo Please implement!
}

void kernels::aderdg::generic::c::volumeUnknownsRestriction(double* luhCoarse,
                               const double* luhFine,
                               const int coarseGridLevel,
                               const int fineGridLevel,
                               const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
                               const int numberOfVariables, const int basisSize){
  // @todo Please implement!
}
#endif
