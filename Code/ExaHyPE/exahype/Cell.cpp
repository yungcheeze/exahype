#include "exahype/Cell.h"
#include "exahype/State.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "tarch/la/ScalarOperations.h"


#include "kernels/KernelCalls.h"

#include "exahype/solvers/Solver.h"

#include "exahype/records/ADERDGCellDescription.h"


tarch::logging::Log exahype::Cell::_log( "exahype::Cell" );


exahype::Cell::Cell():
      Base() {
  _cellData.setADERDGCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}


exahype::Cell::Cell(const Base::DoNotCallStandardConstructor& value):
      Base(value) {
  _cellData.setADERDGCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

exahype::Cell::Cell(const Base::PersistentCell& argument):
                      Base(argument) {
  // do nothing
}


int exahype::Cell::getADERDGCellDescriptionsIndex() const {
  return _cellData.getADERDGCellDescriptionsIndex();
}


void exahype::Cell::init(
  const int                                    level,
  const tarch::la::Vector<DIMENSIONS,double>&  size,
  const tarch::la::Vector<DIMENSIONS,double>&  cellOffset
) {
  assertion1( !ADERDGCellDescriptionHeap::getInstance().isValidIndex(_cellData.getADERDGCellDescriptionsIndex()), toString() );
  const int ADERDGCellDescriptionIndex = ADERDGCellDescriptionHeap::getInstance().createData(0,0);
  _cellData.setADERDGCellDescriptionsIndex( ADERDGCellDescriptionIndex );

  int solverNumber = 0;
  for (
    std::vector<exahype::solvers::Solver*>::const_iterator p = solvers::RegisteredSolvers.begin();
    p != solvers::RegisteredSolvers.end();
    p++
  ) {
    exahype::solvers::Solver* solver = *p;

    // Has to be +1 here
    if (level==solver->getMinimumTreeDepth()+1) {
      switch (solver->getType()) {
        case exahype::solvers::Solver::ADER_DG:
        {
          exahype::records::ADERDGCellDescription newCellDescription;

          // Pass geometry information to the cellDescription description
          newCellDescription.setLevel (level);
          newCellDescription.setSize  (size);
          newCellDescription.setOffset(cellOffset);

          // todo 16/02/25:Dominic Etienne Charrier:
          // This should be set to max_double.
          // Move time step synchronisation in time step synchronisation mapping.
          newCellDescription.setPredictorTimeStamp( solver->getMinPredictorTimeStamp() );
          newCellDescription.setSolverNumber ( solverNumber);

          const int spaceTimeUnknownsPerCell     = solver->getSpaceTimeUnknownsPerCell();
          const int SpaceTimeFluxUnknownsPerCell = solver->getSpaceTimeFluxUnknownsPerCell();
          const int unknownsPerCell              = solver->getUnknownsPerCell();
          const int fluxUnknownsPerCell          = solver->getFluxUnknownsPerCell();
          const int unknownsPerCellBoundary      = solver->getUnknownsPerCellBoundary();

          // Allocate space-time DoF
          newCellDescription.setSpaceTimePredictor (DataHeap::getInstance().createData(spaceTimeUnknownsPerCell,spaceTimeUnknownsPerCell));
          newCellDescription.setSpaceTimeVolumeFlux(DataHeap::getInstance().createData(SpaceTimeFluxUnknownsPerCell,SpaceTimeFluxUnknownsPerCell));

          // Allocate volume DoF
          newCellDescription.setSolution   (DataHeap::getInstance().createData(unknownsPerCell,unknownsPerCell));
          newCellDescription.setUpdate     (DataHeap::getInstance().createData(unknownsPerCell,unknownsPerCell));
          newCellDescription.setPredictor  (DataHeap::getInstance().createData(unknownsPerCell,unknownsPerCell));
          newCellDescription.setVolumeFlux (DataHeap::getInstance().createData(fluxUnknownsPerCell,fluxUnknownsPerCell));

          // Allocate face DoF
          newCellDescription.setExtrapolatedPredictor(DataHeap::getInstance().createData(unknownsPerCellBoundary,unknownsPerCellBoundary));
          newCellDescription.setFluctuation          (DataHeap::getInstance().createData(unknownsPerCellBoundary,unknownsPerCellBoundary));

          ADERDGCellDescriptionHeap::getInstance().getData( _cellData.getADERDGCellDescriptionsIndex() ).push_back( newCellDescription );
        }
        break;
      }
    }
    else {
      logDebug( "init(...)","cell is not associated with any solver. cell=" << toString() << ", level=" << level << ", size=" << size << ",offset=" << cellOffset );
    }
    solverNumber++;
  }
}
