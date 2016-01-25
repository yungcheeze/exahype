#include "exahype/Cell.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "tarch/la/ScalarOperations.h"


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
// ! End of code for multiscalelinkedcell toolbox.






// ! Begin of code for DG method


void exahype::Cell::init(
  const int                                    level,
  const tarch::la::Vector<DIMENSIONS,double>&  size,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre
) {
  // @todo Tobias Weinzierl
  // New Workflow
  // Ask the solver whether there is an ADERDG cell here at this position (we might
  // need the offset as well). If yes, take the corresponding Solver Description and ask
  // it to create the right Cell Description.

  // @todo Tobias Weinzierl
  // Delegate to solver-specific code fragments

//  _cellData.setADERDGCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);

  /*
  // ! Begin of code for multiscalelinkedcell toolbox and DG method
  if (!fineGridCell.isRefined()) {      // We only want to initialize ADERDGCellDescriptions on the initial fine grid
    logDebug("enterCell(...)","initialising ADERDGCellDescription: " << "fine grid level: " << fineGridVerticesEnumerator.getLevel() << ", fine grid position of cell: " << fineGridPositionOfCell);
*/



/*
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
    cellDescriptionForPde.setTimeStamp(0.0);
  }
*/
}
// ! End of code for DG method/multiscalelinkedcell toolbox.


