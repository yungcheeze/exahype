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


void exahype::Cell::setupMetaData() {
  assertion1( !ADERDGCellDescriptionHeap::getInstance().isValidIndex( _cellData.getADERDGCellDescriptionsIndex()), toString() );

  const int ADERDGCellDescriptionIndex =
    ADERDGCellDescriptionHeap::getInstance().createData(0, 0);
  _cellData.setADERDGCellDescriptionsIndex(ADERDGCellDescriptionIndex);
}


void exahype::Cell::shutdownMetaData() {
  assertion1( ADERDGCellDescriptionHeap::getInstance().isValidIndex( _cellData.getADERDGCellDescriptionsIndex() ), toString() );

  ADERDGCellDescriptionHeap::getInstance().deleteData( _cellData.getADERDGCellDescriptionsIndex() );
}


void exahype::Cell::addNewCellDescription(
    const int solverNumber,
    const exahype::records::ADERDGCellDescription::Type cellType,
    const exahype::records::ADERDGCellDescription::RefinementEvent
        refinementEvent,
    const int level,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const tarch::la::Vector<DIMENSIONS, double>& size,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre
) {
  assertion1( ADERDGCellDescriptionHeap::getInstance().isValidIndex( _cellData.getADERDGCellDescriptionsIndex()), toString() );

  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
  switch (solver->getType()) {
    case exahype::solvers::Solver::ADER_DG: {
      exahype::records::ADERDGCellDescription newCellDescription;
      newCellDescription.setSolverNumber(solverNumber);

      // Default AMR settings
      newCellDescription.setType(cellType);
      newCellDescription.setParentIndex(parentIndex);
      newCellDescription.setLevel(level);
      newCellDescription.setRefinementEvent(refinementEvent);

      std::bitset<DIMENSIONS_TIMES_TWO>
          riemannSolvePerformed;  // default construction: no bit set
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
      logDebug("addNewCellDescription(...)",
               "could not add a cell descriptions for this solver. cell="
                   << toString() << ", level=" << level << ", size=" << size
                   << ",offset=" << cellCentre);
    } break;
  }
}

// todo this method is deprecated about to change.
// Will be replaced by
// allocateCellDescriptionFields(std::vector<...>::iterator)
// similar to cleanCellDescriptionFields(std::vector<...>::iterator)
void exahype::Cell::initialiseCellDescription(const int solverNumber) {
  // Ensure that -1 value is invalid index (cf. addNewCellDescription)
  assertion(!DataHeap::getInstance().isValidIndex(-1));

  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 _cellData.getADERDGCellDescriptionsIndex()),
             toString());

  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
  switch (solver->getType()) {
    case exahype::solvers::Solver::ADER_DG: {
      for (std::vector<exahype::records::ADERDGCellDescription>::iterator p =
               ADERDGCellDescriptionHeap::getInstance()
                   .getData(_cellData.getADERDGCellDescriptionsIndex())
                   .begin();
           p !=
           ADERDGCellDescriptionHeap::getInstance()
               .getData(_cellData.getADERDGCellDescriptionsIndex())
               .end();
           p++) {  // Loop over coarse grid cell descriptions.

        if (solverNumber == p->getSolverNumber()) {
          if (p->getType() == exahype::records::ADERDGCellDescription::Cell) {
            if (!DataHeap::getInstance().isValidIndex(p->getSolution())) {
              assertion(!DataHeap::getInstance().isValidIndex(
                  p->getSpaceTimePredictor()));
              assertion(!DataHeap::getInstance().isValidIndex(
                  p->getSpaceTimeVolumeFlux()));
              assertion(
                  !DataHeap::getInstance().isValidIndex(p->getPredictor()));
              assertion(!DataHeap::getInstance().isValidIndex(p->getUpdate()));
              assertion(
                  !DataHeap::getInstance().isValidIndex(p->getVolumeFlux()));

              p->setPredictorTimeStamp(solver->getMinPredictorTimeStamp());

              const int spaceTimeUnknownsPerCell =
                  solver->getSpaceTimeUnknownsPerCell();
              const int spaceTimeFluxUnknownsPerCell =
                  solver->getSpaceTimeFluxUnknownsPerCell();
              const int unknownsPerCell = solver->getUnknownsPerCell();
              const int fluxUnknownsPerCell = solver->getFluxUnknownsPerCell();

              // Allocate space-time DoF
              p->setSpaceTimePredictor(DataHeap::getInstance().createData(
                  spaceTimeUnknownsPerCell, spaceTimeUnknownsPerCell));
              p->setSpaceTimeVolumeFlux(DataHeap::getInstance().createData(
                  spaceTimeFluxUnknownsPerCell, spaceTimeFluxUnknownsPerCell));

              // Allocate volume DoF
              p->setSolution(DataHeap::getInstance().createData(
                  unknownsPerCell, unknownsPerCell));
              p->setUpdate(DataHeap::getInstance().createData(unknownsPerCell,
                                                              unknownsPerCell));
              p->setPredictor(DataHeap::getInstance().createData(
                  unknownsPerCell, unknownsPerCell));
              p->setVolumeFlux(DataHeap::getInstance().createData(
                  fluxUnknownsPerCell, fluxUnknownsPerCell));
            }
          }

          // Allocate face DoF
          if (p->getType() == exahype::records::ADERDGCellDescription::Cell ||
              p->getType() ==
                  exahype::records::ADERDGCellDescription::Ancestor ||
              p->getType() ==
                  exahype::records::ADERDGCellDescription::Descendant) {
            if (!DataHeap::getInstance().isValidIndex(
                    p->getExtrapolatedPredictor())) {
              assertion(
                  !DataHeap::getInstance().isValidIndex(p->getFluctuation()));

              const int unknownsPerCellBoundary =
                  solver->getUnknownsPerCellBoundary();

              p->setExtrapolatedPredictor(DataHeap::getInstance().createData(
                  unknownsPerCellBoundary, unknownsPerCellBoundary));
              p->setFluctuation(DataHeap::getInstance().createData(
                  unknownsPerCellBoundary, unknownsPerCellBoundary));
            }
          }
        }
      }
    } break;
    default: {
      logDebug("initialiseCellDescription(...)",
               "solver is not associated with any cell descriptions of this "
               "cell. cell="
                   << toString() << ", level=" << level << ", size=" << size
                   << ",offset=" << cellCentre);
    } break;
  }
}

std::vector<exahype::records::ADERDGCellDescription>::iterator
exahype::Cell::deleteCellDescription(
    std::vector<exahype::records::ADERDGCellDescription>::iterator p) {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 _cellData.getADERDGCellDescriptionsIndex()),
             toString());

  return ADERDGCellDescriptionHeap::getInstance()
      .getData(_cellData.getADERDGCellDescriptionsIndex())
      .erase(p);
}

// todo Will be renamed to
// deallocateUnusedCellDescriptionFields(std::vector<...>::iterator)
void exahype::Cell::cleanCellDescription(
    std::vector<exahype::records::ADERDGCellDescription>::iterator p) {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 _cellData.getADERDGCellDescriptionsIndex()),
             toString());

  const solvers::Solver* solver =
      solvers::RegisteredSolvers[p->getSolverNumber()];
  switch (solver->getType()) {
    case exahype::solvers::Solver::ADER_DG: {
      switch (p->getType()) {
        case exahype::records::ADERDGCellDescription::Erased:
        case exahype::records::ADERDGCellDescription::EmptyAncestor:
        case exahype::records::ADERDGCellDescription::EmptyDescendant:
        case exahype::records::ADERDGCellDescription::Ancestor:
        case exahype::records::ADERDGCellDescription::Descendant:
          if (DataHeap::getInstance().isValidIndex(p->getSolution())) {
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getSpaceTimePredictor()));
            assertion(!DataHeap::getInstance().isValidIndex(
                p->getSpaceTimeVolumeFlux()));
            assertion(!DataHeap::getInstance().isValidIndex(p->getPredictor()));
            assertion(!DataHeap::getInstance().isValidIndex(p->getUpdate()));
            assertion(
                !DataHeap::getInstance().isValidIndex(p->getVolumeFlux()));

            DataHeap::getInstance().deleteData(p->getSpaceTimePredictor());
            DataHeap::getInstance().deleteData(p->getSpaceTimeVolumeFlux());
            DataHeap::getInstance().deleteData(p->getSolution());
            DataHeap::getInstance().deleteData(p->getUpdate());
            DataHeap::getInstance().deleteData(p->getUpdate());
          }
          break;
        default:
          break;
      }
      switch (p->getType()) {
        case exahype::records::ADERDGCellDescription::Erased:
        case exahype::records::ADERDGCellDescription::EmptyAncestor:
        case exahype::records::ADERDGCellDescription::EmptyDescendant:
          if (!DataHeap::getInstance().isValidIndex(
                  p->getExtrapolatedPredictor())) {
            assertion(
                !DataHeap::getInstance().isValidIndex(p->getFluctuation()));

            DataHeap::getInstance().deleteData(p->getExtrapolatedPredictor());
            DataHeap::getInstance().deleteData(p->getFluctuation());
          }
          break;
        default:
          break;
      }
      break;
      default:
        logDebug("cleanCellDescription(...)",
                 "solver is not associated with any cell descriptions of this "
                 "cell. cell="
                     << toString() << ", level=" << level << ", size=" << size
                     << ",offset=" << cellCentre);
        break;
    }
  }
}

exahype::Cell::SubcellPosition
exahype::Cell::computeSubcellPositionOfCellOrAncestor(
    const exahype::records::ADERDGCellDescription& pChild) const {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 _cellData.getADERDGCellDescriptionsIndex()),
             toString());
  assertion1(solvers::RegisteredSolvers[pChild.getSolverNumber()]->getType() ==
                 exahype::solvers::Solver::ADER_DG,
             toString());
  assertion1(pChild.getType() == exahype::records::ADERDGCellDescription::Cell,
             toString());

  exahype::Cell::SubcellPosition subcellPosition;
  // initialisation
  subcellPosition.parentIndex = pChild.getParentIndex();
  exahype::records::ADERDGCellDescription* pParent = 0;
  for (std::vector<exahype::records::ADERDGCellDescription>::iterator p =
           ADERDGCellDescriptionHeap::getInstance()
               .getData(subcellPosition.parentIndex)
               .begin();
       p !=
       ADERDGCellDescriptionHeap::getInstance()
           .getData(subcellPosition.parentIndex)
           .end();
       p++) {  // Loop over cell descriptions

    if (p->getSolverNumber() == pChild.getSolverNumber()) {
      pParent = &(*p);
      assertion(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
          subcellPosition.parentIndex));
    }
  }

  if (pParent != 0) {
    // recursion
    while (pParent->getType() ==
           exahype::records::ADERDGCellDescription::EmptyAncestor) {
      for (std::vector<exahype::records::ADERDGCellDescription>::iterator p =
               ADERDGCellDescriptionHeap::getInstance()
                   .getData(subcellPosition.parentIndex)
                   .begin();
           p !=
           ADERDGCellDescriptionHeap::getInstance()
               .getData(subcellPosition.parentIndex)
               .end();
           p++) {  // Loop over cell descriptions
        if (p->getSolverNumber() == pChild.getSolverNumber()) {
          pParent = &(*p);
          subcellPosition.parentIndex = p->getParentIndex();
          assertion(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
              subcellPosition.parentIndex));
        }
      }
    }
    assertion(pParent->getType() ==
              exahype::records::ADERDGCellDescription::Ancestor);

    // compute subcell index
    const int levelDelta = pChild.getLevel() - pParent->getLevel();
    double scaling = tarch::la::aPowI(levelDelta, 3);
    for (int xi = 0; xi < DIMENSIONS; xi++) {
      subcellPosition.subcellIndex[xi] = tarch::la::round(
          scaling * pChild.getOffset(xi) - pParent->getOffset(xi));
    }
  }

  return subcellPosition;
}

exahype::Cell::SubcellPosition
exahype::Cell::computeSubcellPositionOfDescendant(
    const exahype::records::ADERDGCellDescription& pChild) const {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 _cellData.getADERDGCellDescriptionsIndex()),
             toString());
  assertion1(solvers::RegisteredSolvers[pChild.getSolverNumber()]->getType() ==
                 exahype::solvers::Solver::ADER_DG,
             toString());
  assertion1(
      pChild.getType() == exahype::records::ADERDGCellDescription::Descendant,
      toString());

  exahype::Cell::SubcellPosition subcellPosition;
  // initialisation
  subcellPosition.parentIndex = pChild.getParentIndex();
  exahype::records::ADERDGCellDescription* pParent = 0;
  assertion(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      subcellPosition.parentIndex));
  for (std::vector<exahype::records::ADERDGCellDescription>::iterator p =
           ADERDGCellDescriptionHeap::getInstance()
               .getData(subcellPosition.parentIndex)
               .begin();
       p !=
       ADERDGCellDescriptionHeap::getInstance()
           .getData(subcellPosition.parentIndex)
           .end();
       p++) {  // Loop over cell descriptions

    if (p->getSolverNumber() == pChild.getSolverNumber()) {
      pParent = &(*p);
      subcellPosition.parentIndex = p->getParentIndex();
      assertion(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
          subcellPosition.parentIndex));
    }
  }

  if (pParent != 0) {
    // recursion
    while (pParent->getType() ==
           exahype::records::ADERDGCellDescription::EmptyDescendant) {
      for (std::vector<exahype::records::ADERDGCellDescription>::iterator p =
               ADERDGCellDescriptionHeap::getInstance()
                   .getData(subcellPosition.parentIndex)
                   .begin();
           p !=
           ADERDGCellDescriptionHeap::getInstance()
               .getData(subcellPosition.parentIndex)
               .end();
           p++) {  // Loop over cell descriptions
        if (p->getSolverNumber() == pChild.getSolverNumber()) {
          pParent = &(*p);
          subcellPosition.parentIndex = p->getParentIndex();
          assertion(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
              subcellPosition.parentIndex));
        }
      }
    }
    assertion(pParent->getType() ==
                  exahype::records::ADERDGCellDescription::Descendant ||
              pParent->getType() ==
                  exahype::records::ADERDGCellDescription::Cell);

    // compute subcell index
    const int levelDelta = pChild.getLevel() - pParent->getLevel();
    double scaling = tarch::la::aPowI(levelDelta, 3);
    for (int xi = 0; xi < DIMENSIONS; xi++) {
      subcellPosition.subcellIndex[xi] = tarch::la::round(
          scaling * pChild.getOffset(xi) - pParent->getOffset(xi));
    }
  } else {
    std::cerr << "exahype::Cell::computeSubcellPositionOfDescendant: parent of "
                 "descendant could not be found!"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  return subcellPosition;
}
