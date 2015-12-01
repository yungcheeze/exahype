#include "EulerFlow3d/Cell.h"

#include "EulerFlow3d/multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "tarch/la/ScalarOperations.h"

exahype::Cell::Cell():
Base() {
  // do nothing
}


exahype::Cell::Cell(const Base::DoNotCallStandardConstructor& value):
          Base(value) {
  _cellData.setCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

exahype::Cell::Cell(const Base::PersistentCell& argument):
          Base(argument) {
  // do nothing
}

// ! Begin of code for multiscalelinkedcell toolbox.
int
exahype::Cell::getCellDescriptionsIndex() const {
  return _cellData.getCellDescriptionsIndex();
}

void
exahype::Cell::initCellWithDefaultValues() {
  _cellData.setCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}
// ! End of code for multiscalelinkedcell toolbox.

// ! Begin of code for DG method/multiscalelinkedcell toolbox
void
exahype::Cell::initCellInComputeTree(const int level,
                                     const tarch::la::Vector<DIMENSIONS,double> size,
                                     const int numberOfPDEs,
                                     const int* order,
                                     const int* numberOfVariables) {
  const int indexOfCellDescriptions = CellDescriptionHeap::getInstance().createData(numberOfPDEs);
  assertion( indexOfCellDescriptions >= 0 );
  _cellData.setCellDescriptionsIndex(indexOfCellDescriptions);

  for (int pde = 0; pde < numberOfPDEs; pde++) {
    int basisSize            = exahype::order[pde]+1;

    //  @todo: consider to precompute these values/to use a lookup table at this point
    int numberOfSpaceTimeDof = EXAHYPE_PATCH_SIZE_TOTAL                            * numberOfVariables[pde]          * tarch::la::aPowI(DIMENSIONS+1,basisSize);
    int numberOfDof          = EXAHYPE_PATCH_SIZE_TOTAL                            * exahype::numberOfVariables[pde] * tarch::la::aPowI(DIMENSIONS  ,basisSize);
    int numberOfDofOnFaces   = EXAHYPE_PATCH_SIZE_TOTAL_TIMES_DIMENSIONS_TIMES_TWO * exahype::numberOfVariables[pde] * tarch::la::aPowI(DIMENSIONS-1,basisSize);

    records::CellDescription& cellDescriptionForPde =
        CellDescriptionHeap::getInstance().getData(indexOfCellDescriptions)[pde];

    // Allocate space-time DoF
    cellDescriptionForPde.setSpaceTimePredictor (DataHeap::getInstance().createData(numberOfSpaceTimeDof,numberOfSpaceTimeDof));
    cellDescriptionForPde.setSpaceTimeVolumeFlux(DataHeap::getInstance().createData(DIMENSIONS*numberOfSpaceTimeDof,DIMENSIONS*numberOfSpaceTimeDof));

    // Allocate volume DoF
    cellDescriptionForPde.setSolution   (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
    cellDescriptionForPde.setUpdate     (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
    cellDescriptionForPde.setPredictor  (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
    cellDescriptionForPde.setVolumeFlux (DataHeap::getInstance().createData(DIMENSIONS*numberOfDof,DIMENSIONS*numberOfDof));

    // Allocate face DoF
    cellDescriptionForPde.setExtrapolatedPredictor(DataHeap::getInstance().createData(numberOfDofOnFaces,numberOfDofOnFaces));
    cellDescriptionForPde.setFluctuation          (DataHeap::getInstance().createData(numberOfDofOnFaces,numberOfDofOnFaces));

    // Pass geometry information to the cellDescription description
    cellDescriptionForPde.setLevel (level);
    cellDescriptionForPde.setSize  (size);
  }
}
// ! End of code for DG method/multiscalelinkedcell toolbox.


