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
    const int solverNumber,
    const exahype::Cell::CellDescriptionType cellType,
    const int level,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const tarch::la::Vector<DIMENSIONS, double>& size,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre) {
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

      // AMR information
      newCellDescription.setType(cellType);
      newCellDescription.setParentIndex(parentIndex);
      newCellDescription.setLevel(level);
      newCellDescription.setFineGridPositionOfCell(fineGridPositionOfCell);
      newCellDescription.setParent(false);
      newCellDescription.setHasNeighboursOfTypeCell(false);
      newCellDescription.setRefinementNecessary(false);
      newCellDescription.setVirtualRefinementNecessary(false);

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
    }
  }
}

void exahype::Cell::initialiseCellDescription(const int solverNumber) {
  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
  switch (solver->getType()) {
    case exahype::solvers::Solver::ADER_DG: {
      // @todo ask Tobias
      // Ensure that default data field index initialisation value is invalid.
      assertion(!DataHeap::getInstance().isValidIndex(-1));

      // calling method should check if index is valid
      assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
          _cellData.getADERDGCellDescriptionsIndex()),
          toString());
      assertion1(static_cast<unsigned int>(solverNumber) <
          ADERDGCellDescriptionHeap::getInstance().getData(
          _cellData.getADERDGCellDescriptionsIndex()).size(),
          toString());

      records::ADERDGCellDescription& cellDescription =
          ADERDGCellDescriptionHeap::getInstance().getData(
              _cellData.getADERDGCellDescriptionsIndex())[solverNumber];

      if (cellDescription.getType()==exahype::Cell::RealCell) {
        if (!DataHeap::getInstance().isValidIndex(
            cellDescription.getSolution())) {
          assertion(!DataHeap::getInstance().isValidIndex(
              cellDescription.getSpaceTimePredictor()));
          assertion(!DataHeap::getInstance().isValidIndex(
              cellDescription.getSpaceTimeVolumeFlux()));
          assertion(!DataHeap::getInstance().isValidIndex(
              cellDescription.getPredictor()));
          assertion(!DataHeap::getInstance().isValidIndex(
              cellDescription.getUpdate()));
          assertion(!DataHeap::getInstance().isValidIndex(
              cellDescription.getVolumeFlux()));

          cellDescription.setPredictorTimeStamp(
              solver->getMinPredictorTimeStamp());

          const int spaceTimeUnknownsPerCell =
              solver->getSpaceTimeUnknownsPerCell();
          const int spaceTimeFluxUnknownsPerCell =
              solver->getSpaceTimeFluxUnknownsPerCell();
          const int unknownsPerCell = solver->getUnknownsPerCell();
          const int fluxUnknownsPerCell = solver->getFluxUnknownsPerCell();

          // Allocate space-time DoF
          cellDescription.setSpaceTimePredictor(
              DataHeap::getInstance().createData(spaceTimeUnknownsPerCell,
                                                 spaceTimeUnknownsPerCell));
          cellDescription.setSpaceTimeVolumeFlux(
              DataHeap::getInstance().createData(
                  spaceTimeFluxUnknownsPerCell,
                  spaceTimeFluxUnknownsPerCell));

          // Allocate volume DoF
          cellDescription.setSolution(DataHeap::getInstance().createData(
              unknownsPerCell, unknownsPerCell));
          cellDescription.setUpdate(DataHeap::getInstance().createData(
              unknownsPerCell, unknownsPerCell));
          cellDescription.setPredictor(DataHeap::getInstance().createData(
              unknownsPerCell, unknownsPerCell));
          cellDescription.setVolumeFlux(DataHeap::getInstance().createData(
              fluxUnknownsPerCell, fluxUnknownsPerCell));
        }
      }

      // Allocate face DoF
      if (cellDescription.getType()==RealCell
          ||
          cellDescription.getType()==RealShell
          ||
          (cellDescription.getType()==VirtualShell
          &&
          cellDescription.getHasNeighboursOfTypeCell())) {

        if (!DataHeap::getInstance().isValidIndex(
            cellDescription.getExtrapolatedPredictor())) {
          assertion(!DataHeap::getInstance().isValidIndex(
              cellDescription.getFluctuation()));

          const int unknownsPerCellBoundary =
              solver->getUnknownsPerCellBoundary();

          cellDescription.setExtrapolatedPredictor(
              DataHeap::getInstance().createData(unknownsPerCellBoundary,
                                                 unknownsPerCellBoundary));
          cellDescription.setFluctuation(DataHeap::getInstance().createData(
              unknownsPerCellBoundary, unknownsPerCellBoundary));
        }
      }
    } break;
  }
}

void exahype::Cell::cleanCellDescription(const int solverNumber) {
  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
  switch (solver->getType()) {
    case exahype::solvers::Solver::ADER_DG: {
      assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
          _cellData.getADERDGCellDescriptionsIndex()),
          toString());

      records::ADERDGCellDescription& cellDescription =
          ADERDGCellDescriptionHeap::getInstance().getData(
              this->getADERDGCellDescriptionsIndex())[solverNumber];

      assertion1(cellDescription.getType()==RealCell
                 ||
                 cellDescription.getType()==RealShell
                 ||
                 cellDescription.getType()==VirtualShell,
                toString());

      if (cellDescription.getType()==RealShell
          ||
          cellDescription.getType()==VirtualShell) {
        if (DataHeap::getInstance().isValidIndex(
            cellDescription.getSolution())) {
          assertion(!DataHeap::getInstance().isValidIndex(
              cellDescription.getSpaceTimePredictor()));
          assertion(!DataHeap::getInstance().isValidIndex(
              cellDescription.getSpaceTimeVolumeFlux()));
          assertion(!DataHeap::getInstance().isValidIndex(
              cellDescription.getPredictor()));
          assertion(!DataHeap::getInstance().isValidIndex(
              cellDescription.getUpdate()));
          assertion(!DataHeap::getInstance().isValidIndex(
              cellDescription.getVolumeFlux()));

          DataHeap::getInstance().deleteData(
              cellDescription.getSpaceTimePredictor());
          DataHeap::getInstance().deleteData(
              cellDescription.getSpaceTimeVolumeFlux());
          DataHeap::getInstance().deleteData(
              cellDescription.getSolution());
          DataHeap::getInstance().deleteData(
              cellDescription.getUpdate());
          DataHeap::getInstance().deleteData(
              cellDescription.getUpdate());
        }
      }

      if(cellDescription.getType()==VirtualShell
         &&
         !cellDescription.getHasNeighboursOfTypeCell()) {
        if (!DataHeap::getInstance().isValidIndex(
            cellDescription.getExtrapolatedPredictor())) {
          assertion(!DataHeap::getInstance().isValidIndex(
              cellDescription.getFluctuation()));

          DataHeap::getInstance().deleteData(
              cellDescription.getExtrapolatedPredictor());
          DataHeap::getInstance().deleteData(
              cellDescription.getFluctuation());
        }
      }
    }
  }
}

void exahype::Cell::init(
    const int level, const tarch::la::Vector<DIMENSIONS, double>& size,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre) {
  assertion1(!ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      _cellData.getADERDGCellDescriptionsIndex()),
      toString());
  const int ADERDGCellDescriptionIndex =
      ADERDGCellDescriptionHeap::getInstance().createData(0, 0);
  _cellData.setADERDGCellDescriptionsIndex(ADERDGCellDescriptionIndex);

  int solverNumber = 0;
  for (std::vector<exahype::solvers::Solver*>::const_iterator p =
      solvers::RegisteredSolvers.begin();
      p != solvers::RegisteredSolvers.end(); p++) {
    exahype::solvers::Solver* solver = *p;

    // Has to be +1 here
    if (level == solver->getMinimumTreeDepth() + 1) {
      switch (solver->getType()) {
        case exahype::solvers::Solver::ADER_DG: {
          exahype::records::ADERDGCellDescription newCellDescription;

          // Pass geometry information to the cellDescription description
          newCellDescription.setLevel(level);
          newCellDescription.setSize(size);
          newCellDescription.setOffset(cellCentre);

          // todo 16/02/25:Dominic Etienne Charrier:
          // This should be set to max_double.
          // Move time step synchronisation in time step synchronisation
          // mapping.
          newCellDescription.setPredictorTimeStamp(
              solver->getMinPredictorTimeStamp());
          newCellDescription.setSolverNumber(solverNumber);

          const int spaceTimeUnknownsPerCell =
              solver->getSpaceTimeUnknownsPerCell();
          const int SpaceTimeFluxUnknownsPerCell =
              solver->getSpaceTimeFluxUnknownsPerCell();
          const int unknownsPerCell = solver->getUnknownsPerCell();
          const int fluxUnknownsPerCell = solver->getFluxUnknownsPerCell();
          const int unknownsPerCellBoundary =
              solver->getUnknownsPerCellBoundary();

          // Allocate space-time DoF
          newCellDescription.setSpaceTimePredictor(
              DataHeap::getInstance().createData(spaceTimeUnknownsPerCell,
                                                 spaceTimeUnknownsPerCell));
          newCellDescription.setSpaceTimeVolumeFlux(
              DataHeap::getInstance().createData(SpaceTimeFluxUnknownsPerCell,
                                                 SpaceTimeFluxUnknownsPerCell));

          // Allocate volume DoF
          newCellDescription.setSolution(DataHeap::getInstance().createData(
              unknownsPerCell, unknownsPerCell));
          newCellDescription.setUpdate(DataHeap::getInstance().createData(
              unknownsPerCell, unknownsPerCell));
          newCellDescription.setPredictor(DataHeap::getInstance().createData(
              unknownsPerCell, unknownsPerCell));
          newCellDescription.setVolumeFlux(DataHeap::getInstance().createData(
              fluxUnknownsPerCell, fluxUnknownsPerCell));

          // Allocate face DoF
          newCellDescription.setExtrapolatedPredictor(
              DataHeap::getInstance().createData(unknownsPerCellBoundary,
                                                 unknownsPerCellBoundary));
          newCellDescription.setFluctuation(DataHeap::getInstance().createData(
              unknownsPerCellBoundary, unknownsPerCellBoundary));

          ADERDGCellDescriptionHeap::getInstance()
          .getData(_cellData.getADERDGCellDescriptionsIndex())
          .push_back(newCellDescription);
        } break;
      }
    } else {
      logDebug("init(...)", "cell is not associated with any solver. cell="
               << toString() << ", level=" << level
               << ", size=" << size
               << ",offset=" << cellCentre);
    }
    solverNumber++;
  }
}
