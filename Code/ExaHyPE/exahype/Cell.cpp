#include "exahype/Cell.h"

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
  _cellData.setADERDGCellDescriptionsIndex( ADERDGCellDescriptionHeap::getInstance().createData(0,0) );

  for ( int solverNumber = 0; solverNumber < static_cast<int>( exahype::solvers::RegisteredSolvers.size() ); solverNumber++) {
    if (level==exahype::solvers::RegisteredSolvers[solverNumber]->getMinimumTreeDepth()) {
      logDebug( "init(...)","initialising cell description: " << "fine grid level: " << fineGridVerticesEnumerator.getLevel() << ", fine grid position of cell: " << fineGridPositionOfCell);

      switch (exahype::solvers::RegisteredSolvers[solverNumber]->getType()) {
        case exahype::solvers::Solver::ADER_DG:
          {
            exahype::records::ADERDGCellDescription newCellDescription;

            // Pass geometry information to the cellDescription description
            newCellDescription.setLevel (level);
            newCellDescription.setSize  (size);
            newCellDescription.setTimeStamp(0.0);
            newCellDescription.setOffset(cellOffset);
            newCellDescription.setSolverNumber(solverNumber);

            int numberOfSpaceTimeDof           = exahype::solvers::RegisteredSolvers[solverNumber]->getNumberOfVariables() * tarch::la::aPowI(DIMENSIONS+1,exahype::solvers::RegisteredSolvers[solverNumber]->getNodesPerCoordinateAxis());
            int numberOfSpaceTimeVolumeFluxDof = DIMENSIONS*numberOfSpaceTimeDof;

            int numberOfDof            = exahype::solvers::RegisteredSolvers[solverNumber]->getNumberOfVariables() * tarch::la::aPowI(DIMENSIONS,exahype::solvers::RegisteredSolvers[solverNumber]->getNodesPerCoordinateAxis());
            int numberOfVolumeFluxDof  = DIMENSIONS * numberOfDof;
            int numberOfDofOnFace      = DIMENSIONS_TIMES_TWO * exahype::solvers::RegisteredSolvers[solverNumber]->getNumberOfVariables() * tarch::la::aPowI(DIMENSIONS-1,exahype::solvers::RegisteredSolvers[solverNumber]->getNodesPerCoordinateAxis());

            // Allocate space-time DoF
            newCellDescription.setSpaceTimePredictor (DataHeap::getInstance().createData(numberOfSpaceTimeDof,numberOfSpaceTimeDof));
            newCellDescription.setSpaceTimeVolumeFlux(DataHeap::getInstance().createData(numberOfSpaceTimeVolumeFluxDof,numberOfSpaceTimeVolumeFluxDof));

            // Allocate volume DoF
            newCellDescription.setSolution   (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
            newCellDescription.setUpdate     (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
            newCellDescription.setPredictor  (DataHeap::getInstance().createData(numberOfDof,numberOfDof));
            newCellDescription.setVolumeFlux (DataHeap::getInstance().createData(numberOfVolumeFluxDof,numberOfVolumeFluxDof));

            // Allocate face DoF
            newCellDescription.setExtrapolatedPredictor(DataHeap::getInstance().createData(numberOfDofOnFace,numberOfDofOnFace));
            newCellDescription.setFluctuation          (DataHeap::getInstance().createData(numberOfDofOnFace,numberOfDofOnFace));

            ADERDGCellDescriptionHeap::getInstance().getData( _cellData.getADERDGCellDescriptionsIndex() ).push_back( newCellDescription );
          }
          break;
      }
    }
  }
}
