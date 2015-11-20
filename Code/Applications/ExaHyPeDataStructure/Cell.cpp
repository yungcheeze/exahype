#include "ExaHyPeDataStructure/Cell.h"

#include "ExaHyPeDataStructure/Constants.h"

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
exahype::Cell::getCellDescriptionIndex() const {
  return _cellData.getCellDescriptionIndex(0);
}
int
exahype::Cell::getCellDescriptionIndex(int elementIndex) const {
  return _cellData.getCellDescriptionIndex(elementIndex);
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
  const int indexOfFirstCellDescription = CellDescriptionHeap::getInstance().createData(NUM_PDE);
  assertion( indexOfFirstCellDescription >= 0 );

  for (int pde = 0; pde < NUM_PDE; pde++) {
    int basisSize            = exahype::order[pde]+1;

    int numberOfSpaceTimeDof = exahype::numberOfVariables[pde] * power(basisSize,DIMENSIONS+1);
    int numberOfDof          = exahype::numberOfVariables[pde] * power(basisSize,DIMENSIONS  );
    int numberOfDofOnFace    = exahype::numberOfVariables[pde] * power(basisSize,DIMENSIONS-1);

    _cellData.setCellDescriptionIndex(pde,indexOfFirstCellDescription+pde);
    records::CellDescription cellDescriptionForPde =
        CellDescriptionHeap::getInstance().getData(indexOfFirstCellDescription)[pde];

    // Allocate space-time DoF
    cellDescriptionForPde.setSpaceTimePredictor (DataHeap::getInstance().createData(numberOfSpaceTimeDof,numberOfSpaceTimeDof));
    cellDescriptionForPde.setSpaceTimeVolumeFlux(DataHeap::getInstance().createData(DIMENSIONS*numberOfSpaceTimeDof,DIMENSIONS*numberOfSpaceTimeDof));

    // Allocate volume DoF
    cellDescriptionForPde.setSolution   (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
    cellDescriptionForPde.setUpdate     (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
    cellDescriptionForPde.setPredictor  (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
    cellDescriptionForPde.setVolumeFlux (DataHeap::getInstance().createData(DIMENSIONS*numberOfDof,DIMENSIONS*numberOfDof));

    // Allocate face DoF
    for (int face=0; face < 6; face++) {
      cellDescriptionForPde.setExtrapolatedPredictor(face,DataHeap::getInstance().createData(numberOfDofOnFace,numberOfDofOnFace));
      cellDescriptionForPde.setFluctuation          (face,DataHeap::getInstance().createData(numberOfDofOnFace,numberOfDofOnFace));
    }

    // Pass geometry information to the cellDescription description
    cellDescriptionForPde.setLevel(level);
    cellDescriptionForPde.setOffset(offset);
    cellDescriptionForPde.setSize(size);
  }
}
// ! End: Required for the multiscalelinkedcell toolbox

