#include "exahype/Cell.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "tarch/la/ScalarOperations.h"

#include "exahype/Constants.h"

exahype::Cell::Cell():
Base() {
  // do nothing
}


exahype::Cell::Cell(const Base::DoNotCallStandardConstructor& value):
                  Base(value) {
  _cellData.setADERDGCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

exahype::Cell::Cell(const Base::PersistentCell& argument):
                  Base(argument) {
  // do nothing
}

// ! Begin of code for multiscalelinkedcell toolbox.
int
exahype::Cell::getADERDGCellDescriptionsIndex() const {
  return _cellData.getADERDGCellDescriptionsIndex();
}

void
exahype::Cell::initCellWithDefaultValues() {
  _cellData.setADERDGCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}
// ! End of code for multiscalelinkedcell toolbox.

// ! Begin of code for DG method/multiscalelinkedcell toolbox
void
exahype::Cell::initCellInComputeTree(const int level,
                                     const tarch::la::Vector<DIMENSIONS,double> size,
                                     const int numberOfPDEs,
                                     const int order,
                                     const int numberOfVariables) {
  const int indexOfADERDGCellDescriptions = ADERDGADERDGCellDescriptionHeap::getInstance().createData(numberOfPDEs);
  assertion( indexOfADERDGCellDescriptions >= 0 );
  _cellData.setADERDGCellDescriptionsIndex(indexOfADERDGCellDescriptions);

  for (int pde = 0; pde < numberOfPDEs; pde++) {
    int basisSize = EXAHYPE_ORDER+1;

    //  @todo: consider to precompute these values/to use a lookup table at this point
    int numberOfSpaceTimeDof           = EXAHYPE_NVARS * tarch::la::aPowI(DIMENSIONS+1,basisSize);
    int numberOfSpaceTimeVolumeFluxDof = DIMENSIONS*numberOfSpaceTimeDof;

    int numberOfDof            = EXAHYPE_NVARS * tarch::la::aPowI(DIMENSIONS  ,basisSize);
    int numberOfVolumeFluxDof  = DIMENSIONS * numberOfDof;
    int numberOfDofOnFace      = DIMENSIONS_TIMES_TWO * EXAHYPE_NVARS * tarch::la::aPowI(DIMENSIONS-1,basisSize);

    records::ADERDGCellDescription& cellDescriptionForPde =
        ADERDGADERDGCellDescriptionHeap::getInstance().getData(indexOfADERDGCellDescriptions)[pde];

    // Allocate space-time DoF
    cellDescriptionForPde.setSpaceTimePredictor (DataHeap::getInstance().createData(numberOfSpaceTimeDof,numberOfSpaceTimeDof));
    cellDescriptionForPde.setSpaceTimeVolumeFlux(DataHeap::getInstance().createData(numberOfSpaceTimeVolumeFluxDof,numberOfSpaceTimeVolumeFluxDof));

    // Allocate volume DoF
    cellDescriptionForPde.setSolution   (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
    cellDescriptionForPde.setUpdate     (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
    cellDescriptionForPde.setPredictor  (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
    cellDescriptionForPde.setVolumeFlux (DataHeap::getInstance().createData(numberOfVolumeFluxDof,numberOfVolumeFluxDof));

    // Allocate face DoF
    cellDescriptionForPde.setExtrapolatedPredictor(DataHeap::getInstance().createData(numberOfDofOnFace,numberOfDofOnFace));
    cellDescriptionForPde.setFluctuation          (DataHeap::getInstance().createData(numberOfDofOnFace,numberOfDofOnFace));

    // Pass geometry information to the cellDescription description
    cellDescriptionForPde.setLevel (level);
    cellDescriptionForPde.setSize  (size);
  }
}
// ! End of code for DG method/multiscalelinkedcell toolbox.


