/*
 * ADERDGAPosterioriSubcellLimiter.cpp
 *
 *  Created on: 6 Oct 2016
 *      Author: dominic
 */

#include "ADERDGAPosterioriSubcellLimiter.h"

#include "tarch/Assertions.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

exahype::solvers::ADERDGAPosterioriSubcellLimiter::ADERDGAPosterioriSubcellLimiter(int aderdgSolverNumber,int finiteVolumesSolverNumber) :
            CellWiseCoupling(0.0,0.0),
            _isActive(false),
            _aderdgSolverNumber(aderdgSolverNumber),
            _finiteVolumesSolverNumber(finiteVolumesSolverNumber)
{
  assertion2(_aderdgSolverNumber >= 0 && _aderdgSolverNumber<static_cast<int>(exahype::solvers::RegisteredSolvers.size()),
      _aderdgSolverNumber,exahype::solvers::RegisteredSolvers.size());
  assertion2(_finiteVolumesSolverNumber >= 0 && _finiteVolumesSolverNumber<static_cast<int>(exahype::solvers::RegisteredSolvers.size()),
      _finiteVolumesSolverNumber,exahype::solvers::RegisteredSolvers.size());
  assertion1(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]->getType()==exahype::solvers::Solver::Type::ADER_DG,
      exahype::solvers::Solver::toString(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]->getType()));
  assertion1(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]->getType()==exahype::solvers::Solver::Type::FiniteVolumes,
        exahype::solvers::Solver::toString(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]->getType()));
}

void exahype::solvers::ADERDGAPosterioriSubcellLimiter::coupleSolversBeforeSolutionUpdate(const int cellDescriptionsIndex) {
  _isActive=false;

  exahype::solvers::ADERDGSolver* aderdgSolver =
      static_cast<exahype::solvers::ADERDGSolver*>(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]);
  exahype::solvers::FiniteVolumesSolver* finiteVolumesSolver =
        static_cast<exahype::solvers::FiniteVolumesSolver*>(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]);
  int aderdgElement        = aderdgSolver->tryGetElement(cellDescriptionsIndex,_aderdgSolverNumber);
  int finiteVolumesElement = finiteVolumesSolver->tryGetElement(cellDescriptionsIndex,_finiteVolumesSolverNumber);

  if (aderdgElement!=exahype::solvers::Solver::NotFound) {
    auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
        cellDescriptionsIndex)[aderdgElement];
    // TODO(Dominic):
    // Later, we can add a new FV solver cell description here if there is none
    // registered on the cell description.
    // The finiteVolumesElement!= assertion will then be replaced by an if.
    assertion2(finiteVolumesElement!=exahype::solvers::Solver::NotFound,cellDescriptionsIndex,_aderdgSolverNumber);
    auto& finiteVolumesCellDescription = exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(
        cellDescriptionsIndex)[finiteVolumesElement];

    aderdgCellDescription.setSkipSolutionUpdate(false);
    finiteVolumesCellDescription.setSkipSolutionUpdate(true);

    //    TODO(JM):
    //    In the first method, you need to check if the aderdgSolver needs limiting (based on min max values etc.).
    bool cellIsTroubled = false;
    if (cellIsTroubled) {
      // TODO(JM):
      // In the first method, you need to check if the aderdgSolver needs limiting (based on min max values etc.).
      // If so project the aderdg data onto the fv subgrid.
      // You get all you need from the cell descriptions

      _isActive=true;
      aderdgCellDescription.setSkipSolutionUpdate(true);
      finiteVolumesCellDescription.setSkipSolutionUpdate(false);
    }
  }
}

void exahype::solvers::ADERDGAPosterioriSubcellLimiter::coupleSolversAfterSolutionUpdate(const int cellDescriptionsIndex) {
  if (_isActive) {
    exahype::solvers::ADERDGSolver* aderdgSolver =
        static_cast<exahype::solvers::ADERDGSolver*>(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]);
    exahype::solvers::FiniteVolumesSolver* finiteVolumesSolver =
          static_cast<exahype::solvers::FiniteVolumesSolver*>(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]);
    int aderdgElement        = aderdgSolver->tryGetElement(cellDescriptionsIndex,_aderdgSolverNumber);
    int finiteVolumesElement = finiteVolumesSolver->tryGetElement(cellDescriptionsIndex,_finiteVolumesSolverNumber);
    assertion2(aderdgElement       !=exahype::solvers::Solver::NotFound,cellDescriptionsIndex,_aderdgSolverNumber);
    assertion2(finiteVolumesElement!=exahype::solvers::Solver::NotFound,cellDescriptionsIndex,_finiteVolumesSolverNumber);

    if (aderdgElement!=exahype::solvers::Solver::NotFound) {
      auto& aderdgCellDescription        = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
          cellDescriptionsIndex)[aderdgElement];
      auto& finiteVolumesCellDescription = exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(
          cellDescriptionsIndex)[finiteVolumesElement];

      // TODO(JM):
      // In the second method,  you project the fv solution onto the aderdg solution again.
      // You get all you need from the cell descriptions
    }
  }
}
