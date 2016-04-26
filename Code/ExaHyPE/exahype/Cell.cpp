#include "exahype/Cell.h"
#include "exahype/State.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "tarch/la/ScalarOperations.h"

#include "kernels/KernelCalls.h"

#include "exahype/solvers/Solver.h"

#include "exahype/records/ADERDGCellDescription.h"

tarch::logging::Log exahype::Cell::_log("exahype::Cell");

exahype::Cell::Cell() : Base() {
  _cellData.setADERDGCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

exahype::Cell::Cell(const Base::DoNotCallStandardConstructor& value)
: Base(value) {
  _cellData.setADERDGCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

exahype::Cell::Cell(const Base::PersistentCell& argument) : Base(argument) {
  // do nothing
}

int exahype::Cell::getADERDGCellDescriptionsIndex() const {
  return _cellData.getADERDGCellDescriptionsIndex();
}

void exahype::Cell::addNewCellDescription(
    const int                                            solverNumber,
    const exahype::records::ADERDGCellDescription::Type  cellType,
    const exahype::records::ADERDGCellDescription::RefinementEvent refinementEvent,
    const exahype::records::ADERDGCellDescription::RefinementEvent augmentationEvent,
    const int                                            level,
    const int                                            parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>&            fineGridPositionOfCell,
    const tarch::la::Vector<DIMENSIONS, double>&         size,
    const tarch::la::Vector<DIMENSIONS, double>&         cellCentre
) {
  if (!ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      _cellData.getADERDGCellDescriptionsIndex())) {
    const int ADERDGCellDescriptionIndex =
        ADERDGCellDescriptionHeap::getInstance().createData(0, 0);
    _cellData.setADERDGCellDescriptionsIndex(ADERDGCellDescriptionIndex);
  }

  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
  switch (solver->getType()) {
  case exahype::solvers::Solver::ADER_DG: {
    exahype::records::ADERDGCellDescription newCellDescription;
    newCellDescription.setSolverNumber(solverNumber);

    // Default AMR settings
    newCellDescription.setType(cellType);
    newCellDescription.setParentIndex(parentIndex);
    newCellDescription.setLevel(level);
    newCellDescription.setFineGridPositionOfCell(fineGridPositionOfCell);
    newCellDescription.setParent(false);
    newCellDescription.setHasNeighboursOfTypeCell(false);
    newCellDescription.setRefinementEvent(refinementEvent);
    newCellDescription.setAugmentationEvent(augmentationEvent);

    std::bitset<DIMENSIONS_TIMES_TWO> riemannSolvePerformed; // default construction: no bit set
    newCellDescription.setRiemannSolvePerformed(riemannSolvePerformed);

    // Pass geometry information to the cellDescription description
    newCellDescription.setSize(size);
    newCellDescription.setOffset(cellCentre);

    // Default field data indices
    newCellDescription.setSpaceTimePredictor(-1);
    newCellDescription.setSpaceTimeVolumeFlux(-1);
    newCellDescription.setPredictor(-1);
    newCellDescription.setSolution(-1);
    newCellDescription.setUpdate(-1);
    newCellDescription.setVolumeFlux(-1);
    newCellDescription.setExtrapolatedPredictor(-1);
    newCellDescription.setFluctuation(-1);

    ADERDGCellDescriptionHeap::getInstance()
    .getData(_cellData.getADERDGCellDescriptionsIndex())
    .push_back(newCellDescription);

  } break;
  default: {
    logDebug("addNewCellDescription(...)", "could not add a cell descriptions for this solver. cell="
        << toString() << ", level=" << level
        << ", size=" << size
        << ",offset=" << cellCentre);
  } break;
  }
}

void exahype::Cell::initialiseCellDescription(const int solverNumber) {
  // Ensure that -1 value is invalid index (cf. addNewCellDescription)
  assertion(!DataHeap::getInstance().isValidIndex(-1));

  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      _cellData.getADERDGCellDescriptionsIndex()),
      toString());

  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
  switch (solver->getType()) {
  case exahype::solvers::Solver::ADER_DG: {
    for (std::vector<exahype::records::ADERDGCellDescription>::
        iterator p = ADERDGCellDescriptionHeap::getInstance().getData(
            _cellData.getADERDGCellDescriptionsIndex()).begin();
        p != ADERDGCellDescriptionHeap::getInstance().getData(
            _cellData.getADERDGCellDescriptionsIndex()).end();
        p++) { // Loop over coarse grid cell descriptions.

      if (solverNumber==p->getSolverNumber()) {
        if (p->getType()==exahype::records::ADERDGCellDescription::Cell) {
          if (!DataHeap::getInstance().isValidIndex(
              p->getSolution())) {
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getSpaceTimePredictor()));
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getSpaceTimeVolumeFlux()));
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getPredictor()));
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getUpdate()));
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getVolumeFlux()));

            p->setPredictorTimeStamp(
                solver->getMinPredictorTimeStamp());

            const int spaceTimeUnknownsPerCell =
                solver->getSpaceTimeUnknownsPerCell();
            const int spaceTimeFluxUnknownsPerCell =
                solver->getSpaceTimeFluxUnknownsPerCell();
            const int unknownsPerCell = solver->getUnknownsPerCell();
            const int fluxUnknownsPerCell = solver->getFluxUnknownsPerCell();

            // Allocate space-time DoF
            p->setSpaceTimePredictor(
                DataHeap::getInstance().createData(spaceTimeUnknownsPerCell,
                    spaceTimeUnknownsPerCell));
            p->setSpaceTimeVolumeFlux(
                DataHeap::getInstance().createData(
                    spaceTimeFluxUnknownsPerCell,
                    spaceTimeFluxUnknownsPerCell));

            // Allocate volume DoF
            p->setSolution(DataHeap::getInstance().createData(
                unknownsPerCell, unknownsPerCell));
            p->setUpdate(DataHeap::getInstance().createData(
                unknownsPerCell, unknownsPerCell));
            p->setPredictor(DataHeap::getInstance().createData(
                unknownsPerCell, unknownsPerCell));
            p->setVolumeFlux(DataHeap::getInstance().createData(
                fluxUnknownsPerCell, fluxUnknownsPerCell));
          }
        }

        // Allocate face DoF
        if (p->getType()==exahype::records::ADERDGCellDescription::Cell
            ||
            p->getType()==exahype::records::ADERDGCellDescription::Shell
            ||
            (p->getType()==exahype::records::ADERDGCellDescription::VirtualShell
                &&
                p->getHasNeighboursOfTypeCell())) {

          if (!DataHeap::getInstance().isValidIndex(
              p->getExtrapolatedPredictor())) {
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getFluctuation()));

            const int unknownsPerCellBoundary =
                solver->getUnknownsPerCellBoundary();

            p->setExtrapolatedPredictor(
                DataHeap::getInstance().createData(unknownsPerCellBoundary,
                    unknownsPerCellBoundary));
            p->setFluctuation(DataHeap::getInstance().createData(
                unknownsPerCellBoundary, unknownsPerCellBoundary));
          }
        }
      }
    }
  } break;
  default: {
    logDebug("initialiseCellDescription(...)", "solver is not associated with any cell descriptions of this cell. cell="
        << toString() << ", level=" << level
        << ", size=" << size
        << ",offset=" << cellCentre);
  } break;
  }
}

void exahype::Cell::cleanCellDescription(const int solverNumber) {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      _cellData.getADERDGCellDescriptionsIndex()),
      toString());

  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
  switch (solver->getType()) {
  case exahype::solvers::Solver::ADER_DG: {
    for (std::vector<exahype::records::ADERDGCellDescription>::
        iterator p = ADERDGCellDescriptionHeap::getInstance().getData(
            _cellData.getADERDGCellDescriptionsIndex()).begin();
        p != ADERDGCellDescriptionHeap::getInstance().getData(
            _cellData.getADERDGCellDescriptionsIndex()).end();
        p++) { // Loop over coarse grid cell descriptions.

      if (solverNumber==p->getSolverNumber()) {
        assertion1(p->getType()==exahype::records::ADERDGCellDescription::Cell
            ||
            p->getType()==exahype::records::ADERDGCellDescription::Shell
            ||
            p->getType()==exahype::records::ADERDGCellDescription::VirtualShell,
            toString());

        if (p->getType()==exahype::records::ADERDGCellDescription::Shell
            ||
            p->getType()==exahype::records::ADERDGCellDescription::VirtualShell) {
          if (DataHeap::getInstance().isValidIndex(
              p->getSolution())) {
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getSpaceTimePredictor()));
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getSpaceTimeVolumeFlux()));
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getPredictor()));
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getUpdate()));
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getVolumeFlux()));

            DataHeap::getInstance().deleteData(
                p->getSpaceTimePredictor());
            DataHeap::getInstance().deleteData(
                p->getSpaceTimeVolumeFlux());
            DataHeap::getInstance().deleteData(
                p->getSolution());
            DataHeap::getInstance().deleteData(
                p->getUpdate());
            DataHeap::getInstance().deleteData(
                p->getUpdate());
          }
        }

        if(p->getType()==exahype::records::ADERDGCellDescription::VirtualShell
            &&
            !p->getHasNeighboursOfTypeCell()) {
          if (!DataHeap::getInstance().isValidIndex(
              p->getExtrapolatedPredictor())) {
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getFluctuation()));

            DataHeap::getInstance().deleteData(
                p->getExtrapolatedPredictor());
            DataHeap::getInstance().deleteData(
                p->getFluctuation());
          }
        }
      }
    }
  } break;
  default:
    logDebug("cleanCellDescription(...)", "solver is not associated with any cell descriptions of this cell. cell="
        << toString() << ", level=" << level
        << ", size=" << size
        << ",offset=" << cellCentre);
    break;
  }
}

exahype::Cell::SubcellPosition exahype::Cell::getSubcellPositionOfVirtualShell(
    const exahype::records::ADERDGCellDescription& pChild) const {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      _cellData.getADERDGCellDescriptionsIndex()),
      toString());
  assertion1(
      solvers::RegisteredSolvers[pChild.getSolverNumber()]->getType()==
          exahype::solvers::Solver::ADER_DG,toString());
  assertion1(
      pChild->getType()==exahype::records::ADERDGCellDescription::VirtualShell,
      toString());

  exahype::Cell::SubcellPosition subcellPosition;
  // initialisation
  subcellPosition.parentIndex = pChild.getParentIndex();
  exahype::records::ADERDGCellDescription* pParent;
  assertion(ADERDGCellDescriptionHeap::getInstance().
            isValidIndex(subcellPosition.parentIndex));
  for (std::vector<exahype::records::ADERDGCellDescription>::
      iterator p = ADERDGCellDescriptionHeap::getInstance().getData(
          subcellPosition.parentIndex).begin();
      p != ADERDGCellDescriptionHeap::getInstance().getData(
          subcellPosition.parentIndex).end();
      p++) { // Loop over cell descriptions

    if (p->getSolverNumber()==pChild.getSolverNumber()) {
      pParent = &(*p);
      subcellPosition.parentIndex = p->getParentIndex();
      assertion(ADERDGCellDescriptionHeap::getInstance().
          isValidIndex(subcellPosition.parentIndex));
    }
  }

  // recursion
  int i = 1;
  while (
      (!pParent->getHasNeighboursOfTypeCell()
       &&
       pParent->getType()==exahype::records::ADERDGCellDescription::VirtualShell)
      ||
      pParent->getType()!=exahype::records::ADERDGCellDescription::Cell
  ) {
    for (std::vector<exahype::records::ADERDGCellDescription>::
        iterator p = ADERDGCellDescriptionHeap::getInstance().getData(
            subcellPosition.parentIndex).begin();
        p != ADERDGCellDescriptionHeap::getInstance().getData(
            subcellPosition.parentIndex).end();
        p++) { // Loop over cell descriptions
      if (p->getSolverNumber()==pChild.getSolverNumber()) {
        pParent = &(*p);
        subcellPosition.parentIndex = p->getParentIndex();
        assertion(ADERDGCellDescriptionHeap::getInstance().
            isValidIndex(subcellPosition.parentIndex));
      }
    }
  }
  assertion((p->getHasNeighboursOfTypeCell()
      &&
      p->getType()==exahype::records::ADERDGCellDescription::VirtualShell)
      ||
      p->getType()==exahype::records::ADERDGCellDescription::Cell);

  // compute subcell index
  const int levelDelta = pChild.getLevel() - pParent->getLevel();
  double scaling = tarch::la::aPowI(levelDelta,3);
  for (int xi=0; xi < DIMENSIONS; xi++) {
    subcellPosition.subcellIndex[xi] = tarch::la::round(
        scaling * pChild.getOffset(xi) - pParent->getOffset(xi)));
  }

  return subcellPosition;
}

//exahype::Cell::SubcellPosition exahype::Cell::getSubcellPositionOfVirtualShell(
//    const int solverNumber,
//    const int parentIndex,
//    const tarch::la::Vector<DIMENSIONS,int>& fineGridPositionOfCell) const {
//
//  assertion1(static_cast<unsigned int>(solverNumber) <
//      ADERDGCellDescriptionHeap::getInstance().getData(
//          _cellData.getADERDGCellDescriptionsIndex()).size(),
//          toString());
//  exahype::Cell::SubcellPosition subcellPosition;
//
//  // initialise recursion
//  for (int d=0; d < DIMENSIONS; d++) {
//    subcellPosition.subcellIndex[d] = fineGridPositionOfCell[d];
//  }
//
//  subcellPosition.parentIndex = parentIndex;
//  assertion(ADERDGCellDescriptionHeap::getInstance().
//      isValidIndex(subcellPosition.parentIndex));
//  exahype::records::ADERDGCellDescription* p = &(ADERDGCellDescriptionHeap::getInstance().
//      getData(subcellPosition.parentIndex)[solverNumber]);
//
//  // recursion
//  int i = 1;
//  while (
//      (!p->getHasNeighboursOfTypeCell() && p->getType()==exahype::records::ADERDGCellDescription::VirtualShell)
//      ||
//      p->getType()!=exahype::records::ADERDGCellDescription::Cell
//  ) {
//    for (int d=0; d < DIMENSIONS; d++) {
//      subcellPosition.subcellIndex[d] +=
//          tarch::la::aPowI(i,3) * p->getFineGridPositionOfCell(d);
//    }
//
//    subcellPosition.parentIndex = p->getParentIndex();
//    assertion(ADERDGCellDescriptionHeap::getInstance().
//        isValidIndex(subcellPosition.parentIndex));
//    p = &(ADERDGCellDescriptionHeap::getInstance().
//        getData(subcellPosition.parentIndex)[solverNumber]);
//    i++;
//  }
//
//  assertion((p->getHasNeighboursOfTypeCell() &&
//      p->getType()==exahype::records::ADERDGCellDescription::VirtualShell)
//      ||
//      p->getType()==exahype::records::ADERDGCellDescription::Cell)
//  return subcellPosition;
//}

//void exahype::Cell::init(
//    const int level, const tarch::la::Vector<DIMENSIONS, double>& size,
//    const tarch::la::Vector<DIMENSIONS, double>& cellCentre) {
//  assertion1(!ADERDGCellDescriptionHeap::getInstance().isValidIndex(
//      _cellData.getADERDGCellDescriptionsIndex()),
//      toString());
//  const int ADERDGCellDescriptionIndex =
//      ADERDGCellDescriptionHeap::getInstance().createData(0, 0);
//  _cellData.setADERDGCellDescriptionsIndex(ADERDGCellDescriptionIndex);
//
//  int solverNumber = 0;
//  for (std::vector<exahype::solvers::Solver*>::const_iterator p =
//      solvers::RegisteredSolvers.begin();
//      p != solvers::RegisteredSolvers.end(); p++) {
//    exahype::solvers::Solver* solver = *p;
//
//    // Has to be +1 here
//    if (level == solver->getMinimumTreeDepth() + 1) {
//      switch (solver->getType()) {
//      case exahype::solvers::Solver::ADER_DG: {
//        exahype::records::ADERDGCellDescription newCellDescription;
//
//        // Pass geometry information to the cellDescription description
//        newCellDescription.setLevel(level);
//        newCellDescription.setSize(size);
//        newCellDescription.setOffset(cellCentre);
//
//        // todo 16/02/25:Dominic Etienne Charrier:
//        // This should be set to max_double.
//        // Move time step synchronisation in time step synchronisation
//        // mapping.
//        newCellDescription.setPredictorTimeStamp(
//            solver->getMinPredictorTimeStamp());
//        newCellDescription.setSolverNumber(solverNumber);
//
//        const int spaceTimeUnknownsPerCell =
//            solver->getSpaceTimeUnknownsPerCell();
//        const int SpaceTimeFluxUnknownsPerCell =
//            solver->getSpaceTimeFluxUnknownsPerCell();
//        const int unknownsPerCell = solver->getUnknownsPerCell();
//        const int fluxUnknownsPerCell = solver->getFluxUnknownsPerCell();
//        const int unknownsPerCellBoundary =
//            solver->getUnknownsPerCellBoundary();
//
//        // Allocate space-time DoF
//        newCellDescription.setSpaceTimePredictor(
//            DataHeap::getInstance().createData(spaceTimeUnknownsPerCell,
//                spaceTimeUnknownsPerCell));
//        newCellDescription.setSpaceTimeVolumeFlux(
//            DataHeap::getInstance().createData(SpaceTimeFluxUnknownsPerCell,
//                SpaceTimeFluxUnknownsPerCell));
//
//        // Allocate volume DoF
//        newCellDescription.setSolution(DataHeap::getInstance().createData(
//            unknownsPerCell, unknownsPerCell));
//        newCellDescription.setUpdate(DataHeap::getInstance().createData(
//            unknownsPerCell, unknownsPerCell));
//        newCellDescription.setPredictor(DataHeap::getInstance().createData(
//            unknownsPerCell, unknownsPerCell));
//        newCellDescription.setVolumeFlux(DataHeap::getInstance().createData(
//            fluxUnknownsPerCell, fluxUnknownsPerCell));
//
//        // Allocate face DoF
//        newCellDescription.setExtrapolatedPredictor(
//            DataHeap::getInstance().createData(unknownsPerCellBoundary,
//                unknownsPerCellBoundary));
//        newCellDescription.setFluctuation(DataHeap::getInstance().createData(
//            unknownsPerCellBoundary, unknownsPerCellBoundary));
//
//        ADERDGCellDescriptionHeap::getInstance()
//        .getData(_cellData.getADERDGCellDescriptionsIndex())
//        .push_back(newCellDescription);
//      } break;
//      }
//    } else {
//      logDebug("init(...)", "cell is not associated with any solver. cell="
//          << toString() << ", level=" << level
//          << ", size=" << size
//          << ",offset=" << cellCentre);
//    }
//    solverNumber++;
//  }
//}
