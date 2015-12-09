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
                                     const int order,
                                     const int numberOfVariables) {
  const int indexOfCellDescriptions = CellDescriptionHeap::getInstance().createData(numberOfPDEs);
  assertion( indexOfCellDescriptions >= 0 );
  _cellData.setCellDescriptionsIndex(indexOfCellDescriptions);

  for (int pde = 0; pde < numberOfPDEs; pde++) {
    int basisSize = EXAHYPE_ORDER+1;

    //  @todo: consider to precompute these values/to use a lookup table at this point
    int numberOfSpaceTimeDof           = EXAHYPE_NVARS * tarch::la::aPowI(DIMENSIONS+1,basisSize);
    int numberOfSpaceTimeVolumeFluxDof = DIMENSIONS*numberOfSpaceTimeDof;

    int numberOfDof            = EXAHYPE_NVARS * tarch::la::aPowI(DIMENSIONS  ,basisSize);
    int numberOfVolumeFluxDof  = DIMENSIONS * numberOfDof;
    int numberOfDofOnFace     = DIMENSIONS_TIMES_TWO * EXAHYPE_NVARS * tarch::la::aPowI(DIMENSIONS-1,basisSize);

    records::CellDescription& cellDescriptionForPde =
        CellDescriptionHeap::getInstance().getData(indexOfCellDescriptions)[pde];

    for (int i=0; i<EXAHYPE_PATCH_SIZE_X+2; i++) { // loop over patches
      for (int j=0; j<EXAHYPE_PATCH_SIZE_Y+2; j++) {
        const int patchIndex = i + (EXAHYPE_PATCH_SIZE_X+2) * j;

        // Initialize patch space-time DoF to invalid value
        cellDescriptionForPde.setSpaceTimePredictor (patchIndex,multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
        cellDescriptionForPde.setSpaceTimeVolumeFlux(patchIndex,multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);

        // Initialize patch volume DoF to invalid value
        cellDescriptionForPde.setSolution   (patchIndex,multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
        cellDescriptionForPde.setUpdate     (patchIndex,multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
        cellDescriptionForPde.setPredictor  (patchIndex,multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
        cellDescriptionForPde.setVolumeFlux (patchIndex,multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);

        // Allocate space-time and volume DoF only for non-ghost patch cells
        if (i > 0 || i < EXAHYPE_PATCH_SIZE_X+1) {
          if (j > 0 || j < EXAHYPE_PATCH_SIZE_Y+1) {
            // Allocate space-time DoF
            cellDescriptionForPde.setSpaceTimePredictor (patchIndex,DataHeap::getInstance().createData(numberOfSpaceTimeDof,numberOfSpaceTimeDof));
            cellDescriptionForPde.setSpaceTimeVolumeFlux(patchIndex,DataHeap::getInstance().createData(numberOfSpaceTimeVolumeFluxDof,numberOfSpaceTimeVolumeFluxDof));

            // Allocate volume DoF
            cellDescriptionForPde.setSolution   (patchIndex,DataHeap::getInstance().createData(numberOfDof,numberOfDof));
            cellDescriptionForPde.setUpdate     (patchIndex,DataHeap::getInstance().createData(numberOfDof,numberOfDof));
            cellDescriptionForPde.setPredictor  (patchIndex,DataHeap::getInstance().createData(numberOfDof,numberOfDof));
            cellDescriptionForPde.setVolumeFlux (patchIndex,DataHeap::getInstance().createData(numberOfVolumeFluxDof,numberOfVolumeFluxDof));
          }
        }

        // Allocate face DoF
        cellDescriptionForPde.setExtrapolatedPredictor(patchIndex,DataHeap::getInstance().createData(numberOfDofOnFace,numberOfDofOnFace));
        cellDescriptionForPde.setFluctuation          (patchIndex,DataHeap::getInstance().createData(numberOfDofOnFace,numberOfDofOnFace));
      }
    }

    // Pass geometry information to the cellDescription description
    cellDescriptionForPde.setLevel (level);
    cellDescriptionForPde.setSize  (size);
  }
}
// ! End of code for DG method/multiscalelinkedcell toolbox.


