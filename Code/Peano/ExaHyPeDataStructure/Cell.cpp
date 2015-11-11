#include "ExaHyPeDataStructure/Cell.h"

#include "ExaHyPeDataStructure/PdeInfo.h"

exahype::Cell::Cell():
Base() {
  // do nothing
}


exahype::Cell::Cell(const Base::DoNotCallStandardConstructor& value):
      Base(value) {
  // do nothing
}

exahype::Cell::Cell(const Base::PersistentCell& argument):
      Base(argument) {
  // do nothing
}

// ! Start: Required for the multiscalelinkedcell toolbox
int
exahype::Cell::getPatchIndex() const {
  return _cellData.getPatchIndex(0);
}
int
exahype::Cell::getPatchIndex(int elementIndex) const {
  return _cellData.getPatchIndex(elementIndex);
}

int power(const int base,int const exponent) {
  int power=1;
  for (int i=0; i<exponent; i++) {
    power *= base;
  }
  return power;
}

void
exahype::Cell::initCellInComputeTree(const int level,
    const tarch::la::Vector<DIMENSIONS,double> offset,
    const tarch::la::Vector<DIMENSIONS,double> size) {
  const int indexOfFirstPatch = PatchDescriptionHeap::getInstance().createData(NUM_PDE);
  assertion( indexOfFirstPatch >= 0 );

  for (int pde = 0; pde < NUM_PDE; pde++) {
    int basisSize            = exahype::order[pde]+1;

    int numberOfSpaceTimeDof = exahype::numberOfVariables[pde] * power(basisSize,DIMENSIONS+1);
    int numberOfDof          = exahype::numberOfVariables[pde] * power(basisSize,DIMENSIONS  );
    int numberOfDofOnFaces   = exahype::numberOfVariables[pde] * power(basisSize,DIMENSIONS-1) * DIMENSIONS_TIMES_TWO;

    _cellData.setPatchIndex(pde,indexOfFirstPatch+pde);
    records::PatchDescription patchDescriptionForPde =
        PatchDescriptionHeap::getInstance().getData(indexOfFirstPatch)[pde];

    // Allocate space-time DoF
    patchDescriptionForPde.setSpaceTimePredictorDof   (DataHeap::getInstance().createData(numberOfSpaceTimeDof,numberOfSpaceTimeDof));
    patchDescriptionForPde.setSpaceTimeVolumeFluxesDof(DataHeap::getInstance().createData(numberOfSpaceTimeDof,numberOfSpaceTimeDof));

    // Allocate volume DoF
    patchDescriptionForPde.setDof                     (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
    patchDescriptionForPde.setUpdateDof               (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
    patchDescriptionForPde.setPredictorDof            (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
    patchDescriptionForPde.setVolumeFluxesDof         (DataHeap::getInstance().createData(numberOfDof,numberOfDof));

    // Allocate face DoF
    patchDescriptionForPde.setExtrapolatedPredictorDof(DataHeap::getInstance().createData(numberOfDofOnFaces,numberOfDofOnFaces));
    patchDescriptionForPde.setNormalFluxesDof         (DataHeap::getInstance().createData(numberOfDofOnFaces,numberOfDofOnFaces));
    patchDescriptionForPde.setFluctuationsDof         (DataHeap::getInstance().createData(numberOfDofOnFaces,numberOfDofOnFaces));

    // Pass geometry information to the patch description
    patchDescriptionForPde.setLevel(level);
    patchDescriptionForPde.setOffset(offset);
    patchDescriptionForPde.setSize(size);
  }
}
// ! End: Required for the multiscalelinkedcell toolbox

