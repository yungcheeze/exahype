/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#include "exahype/mappings/SpaceTimePredictor.h"
#include "exahype/solvers/ADERDGSolver.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/utils/Globals.h"

#include "tarch/multicore/Loop.h"

#include "tarch/multicore/Lock.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include <algorithm>

peano::CommunicationSpecification
exahype::mappings::SpaceTimePredictor::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}

peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

// The remainder specs all are nop
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::SpaceTimePredictor::_log(
    "exahype::mappings::SpaceTimePredictor");

#ifdef Parallel
int exahype::mappings::SpaceTimePredictor::_mpiTag =
    tarch::parallel::Node::reserveFreeTag(
        "exahype::mappings::SpaceTimePredictor");
#endif


tarch::multicore::BooleanSemaphore exahype::mappings::SpaceTimePredictor::_semaphoreForRestriction;

#if defined(Debug)
int exahype::mappings::SpaceTimePredictor::
_parentOfCellOrAncestorNotFound = 0;
int exahype::mappings::SpaceTimePredictor::
_parentOfCellOrAncestorFound    = 0;
int exahype::mappings::SpaceTimePredictor::
_parentOfDescendantFound        = 0;
#endif

void exahype::mappings::SpaceTimePredictor::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (fineGridCell.isInitialised()) {
    const int numberOfADERDGCellDescriptions = static_cast<int>(
        ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).size());

    // please use a different UserDefined per mapping/event
    const peano::datatraversal::autotuning::MethodTrace methodTrace =
        peano::datatraversal::autotuning::UserDefined7;
    const int grainSize =
        peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
            numberOfADERDGCellDescriptions, methodTrace);
    pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
      auto& pFine = fineGridCell.getADERDGCellDescription(i);
      // TODO(Dominic): Might need to only allow a certain
      // number of refinement events for Dyn. AMR.
      // In static AMR, we assume always that the refinement event is None.

      // Reset helper variables
      for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
        pFine.setRiemannSolvePerformed(faceIndex,false);

        // TODO(Dominic): Add to docu.
        #ifdef Parallel
        int listingsOfRemoteRank = countListingsOfRemoteRankAtFace(faceIndex,fineGridVertices,fineGridVerticesEnumerator);
        if (listingsOfRemoteRank==0) {
          listingsOfRemoteRank = TWO_POWER_D;
        }
        pFine.setFaceDataExchangeCounter(faceIndex,listingsOfRemoteRank);
        #endif
      }

      exahype::Cell::SubcellPosition subcellPosition;
      switch (pFine.getType()) {
      case exahype::records::ADERDGCellDescription::Cell:
        assertion1(pFine.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,pFine.toString());

        fineGridCell.validateNoNansInADERDGSolver(i,fineGridCell,fineGridVerticesEnumerator,"exahype::mappings::SpaceTimePredictor::enterCell[pre]");
        computeSpaceTimePredictorAndVolumeIntegral(pFine);
        fineGridCell.validateNoNansInADERDGSolver(i,fineGridCell,fineGridVerticesEnumerator,"exahype::mappings::SpaceTimePredictor::enterCell[pre]");
        break;
      case exahype::records::ADERDGCellDescription::Ancestor:
        assertion1(pFine.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,pFine.toString());
        prepareAncestor(pFine);
        break;
      case exahype::records::ADERDGCellDescription::Descendant:
        assertion1(pFine.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,pFine.toString());
        assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(pFine.getParentIndex()),pFine.toString());
        subcellPosition = fineGridCell.computeSubcellPositionOfDescendant(pFine);
        prolongateADERDGFaceData(pFine, subcellPosition.parentIndex,
            subcellPosition.subcellIndex);
        break;
      default:
        break;
      }
    endpfor

    peano::datatraversal::autotuning::Oracle::getInstance()
        .parallelSectionHasTerminated(methodTrace);
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::SpaceTimePredictor::computeSpaceTimePredictorAndVolumeIntegral(
    exahype::records::ADERDGCellDescription& p) {
  exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
              exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);

  assertion1(p.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getSpaceTimePredictor()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getSpaceTimeVolumeFlux()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getSolution()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getUpdate()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getPredictor()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getVolumeFlux()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getExtrapolatedPredictor()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getFluctuation()),p.toString());

  // space-time DoF (basisSize**(DIMENSIONS+1))
  double* lQi = DataHeap::getInstance().getData(p.getSpaceTimePredictor()).data();
  double* lFi = DataHeap::getInstance().getData(p.getSpaceTimeVolumeFlux()).data();

  // volume DoF (basisSize**(DIMENSIONS))
  double* luh  = DataHeap::getInstance().getData(p.getSolution()).data();
  double* lduh = DataHeap::getInstance().getData(p.getUpdate()).data();
  double* lQhi = DataHeap::getInstance().getData(p.getPredictor()).data();
  double* lFhi = DataHeap::getInstance().getData(p.getVolumeFlux()).data();

  // face DoF (basisSize**(DIMENSIONS-1))
  double* lQhbnd = DataHeap::getInstance().getData(p.getExtrapolatedPredictor()).data();
  double* lFhbnd = DataHeap::getInstance().getData(p.getFluctuation()).data();

  solver->spaceTimePredictor(
      lQi, lFi, lQhi, lFhi,
      lQhbnd,
      lFhbnd,
      luh, p.getSize(),
      p.getPredictorTimeStepSize());
  solver->volumeIntegral(lduh, lFhi, p.getSize());
}

void exahype::mappings::SpaceTimePredictor::prepareAncestor(
    exahype::records::ADERDGCellDescription& p) {
  exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
      exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);

  assertion1(p.getType()==exahype::records::ADERDGCellDescription::Ancestor,p.toString());
  std::fill_n(DataHeap::getInstance().getData(p.getExtrapolatedPredictor()).data(),
      solver->getUnknownsPerCellBoundary(), 0.0);
  std::fill_n(DataHeap::getInstance().getData(p.getFluctuation()).data(),
      solver->getUnknownsPerCellBoundary(), 0.0);
}

void exahype::mappings::SpaceTimePredictor::prolongateADERDGFaceData(
    const exahype::records::ADERDGCellDescription& cellDescription,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  // todo not dynamic with respect to the solver registry
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(parentIndex),
             cellDescription.toString());

  exahype::records::ADERDGCellDescription& cellDescriptionParent =
      ADERDGCellDescriptionHeap::getInstance().getData(
          parentIndex)[cellDescription.getSolverNumber()];

#if defined(Debug)
  _parentOfDescendantFound += 1;
#endif

  assertion(cellDescriptionParent.getSolverNumber() ==
            cellDescription.getSolverNumber());
  assertion(cellDescriptionParent.getType() ==
                exahype::records::ADERDGCellDescription::Cell ||
            cellDescriptionParent.getType() ==
                exahype::records::ADERDGCellDescription::Descendant);

  const int levelFine = cellDescription.getLevel();
  const int levelCoarse = cellDescriptionParent.getLevel();
  assertion(levelCoarse < levelFine);
  const int levelDelta = levelFine - levelCoarse;

  for (int d = 0; d < DIMENSIONS; ++d) {
    // Check if cell is at "left" or "right" d face of parent
    if (subcellIndex[d] == 0) {
      const int faceIndex = 2 * d;

      exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
        exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);

      const int numberOfFaceDof = solver->getUnknownsPerFace();

      double* lQhbndFine = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
          (faceIndex * numberOfFaceDof);
      double* lFhbndFine = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
          (faceIndex * numberOfFaceDof);
      const double* lQhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getExtrapolatedPredictor()).data() +
          (faceIndex * numberOfFaceDof);
      const double* lFhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getFluctuation()).data() +
          (faceIndex * numberOfFaceDof);

      solver->faceUnknownsProlongation(lQhbndFine, lFhbndFine, lQhbndCoarse,
                                       lFhbndCoarse, levelCoarse, levelFine,
                                       getSubfaceIndex(subcellIndex, d));

    } else if (subcellIndex[d] == tarch::la::aPowI(levelDelta, 3) - 1) {
      const int faceIndex = 2 * d + 1;

      exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
          exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);

      const int numberOfFaceDof = solver->getUnknownsPerFace();

      double* lQhbndFine = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
           (faceIndex * numberOfFaceDof);
      double* lFhbndFine = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
           (faceIndex * numberOfFaceDof);
      const double* lQhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getExtrapolatedPredictor()).data() +
          (faceIndex * numberOfFaceDof);
      const double* lFhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getFluctuation()).data() +
          (faceIndex * numberOfFaceDof);

      solver->faceUnknownsProlongation(lQhbndFine, lFhbndFine, lQhbndCoarse,
                                       lFhbndCoarse, levelCoarse, levelFine,
                                       getSubfaceIndex(subcellIndex, d));
    }
  }
}

void exahype::mappings::SpaceTimePredictor::prolongateFiniteVolumesFaceData(
    const exahype::records::FiniteVolumesCellDescription& cellDescription,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  // todo not dynamic with respect to the solver registry
  assertion1(FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(parentIndex),
             cellDescription.toString());

  exahype::records::FiniteVolumesCellDescription& cellDescriptionParent =
      FiniteVolumesCellDescriptionHeap::getInstance().getData(
          parentIndex)[cellDescription.getSolverNumber()];

#if defined(Debug)
  _parentOfDescendantFound += 1;
#endif

  assertion(cellDescriptionParent.getSolverNumber() ==
            cellDescription.getSolverNumber());
//  assertion(cellDescriptionParent.getType() ==
//                exahype::records::FiniteVolumesCellDescription::Cell ||
//            cellDescriptionParent.getType() ==
//                exahype::records::FiniteVolumesCellDescription::Descendant);

  const int levelFine = cellDescription.getLevel();
  const int levelCoarse = cellDescriptionParent.getLevel();
  assertion(levelCoarse < levelFine);
  const int levelDelta = levelFine - levelCoarse;

  for (int d = 0; d < DIMENSIONS; ++d) {
    // Check if cell is at "left" or "right" d face of parent
    if (subcellIndex[d] == 0) {
//      const int faceIndex = 2 * d;
//
//      exahype::solvers::FiniteVolumesSolver* solver = static_cast<exahype::solvers::FiniteVolumesSolver*>(
//        exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);


      // todo @Dominic: implement !!!
//      const int numberOfFaceDof = solver->getUnknownsPerFace();
//
//      double* lQhbndFine = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
//          (faceIndex * numberOfFaceDof);
//      double* lFhbndFine = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
//          (faceIndex * numberOfFaceDof);
//      const double* lQhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getExtrapolatedPredictor()).data() +
//          (faceIndex * numberOfFaceDof);
//      const double* lFhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getFluctuation()).data() +
//          (faceIndex * numberOfFaceDof);
//
//      solver->faceUnknownsProlongation(lQhbndFine, lFhbndFine, lQhbndCoarse,
//                                       lFhbndCoarse, levelCoarse, levelFine,
//                                       getSubfaceIndex(subcellIndex, d));

    } else if (subcellIndex[d] == tarch::la::aPowI(levelDelta, 3) - 1) {
//      const int faceIndex = 2 * d + 1;
//
//      exahype::solvers::FiniteVolumesSolver* solver = static_cast<exahype::solvers::FiniteVolumesSolver*>(
//          exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);

      // todo implement
//      const int numberOfFaceDof = solver->getUnknownsPerFace();
//
//      double* lQhbndFine = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
//           (faceIndex * numberOfFaceDof);
//      double* lFhbndFine = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
//           (faceIndex * numberOfFaceDof);
//      const double* lQhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getExtrapolatedPredictor()).data() +
//          (faceIndex * numberOfFaceDof);
//      const double* lFhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getFluctuation()).data() +
//          (faceIndex * numberOfFaceDof);
//
//      solver->faceUnknownsProlongation(lQhbndFine, lFhbndFine, lQhbndCoarse,
//                                       lFhbndCoarse, levelCoarse, levelFine,
//                                       getSubfaceIndex(subcellIndex, d));
    }
  }
}

tarch::la::Vector<DIMENSIONS - 1, int>
exahype::mappings::SpaceTimePredictor::getSubfaceIndex(
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex, const int d) {
  tarch::la::Vector<DIMENSIONS - 1, int> subfaceIndex;

  int i = 0;
  for (int j = 0; j < DIMENSIONS; j++) {
    if (j != d) {
      subfaceIndex[i] = subcellIndex[j];
      i++;
    }
  }

  return subfaceIndex;
}

void exahype::mappings::SpaceTimePredictor::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("leaveCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
          fineGridCell.getCellDescriptionsIndex())) {
    const int numberOfADERDGCellDescriptions = static_cast<int>(
        ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).size());
    // please use a different UserDefined per mapping/event
    const peano::datatraversal::autotuning::MethodTrace methodTrace = peano::datatraversal::autotuning::UserDefined10;
    const int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfADERDGCellDescriptions, methodTrace);
    pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
        records::ADERDGCellDescription& pFine = fineGridCell.getADERDGCellDescription(i);

    // if we have at least one parent
    if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(pFine.getParentIndex())) {
      for (auto& pParent : ADERDGCellDescriptionHeap::getInstance().getData(pFine.getParentIndex())) {
        exahype::Cell::SubcellPosition subcellPosition;

        if (pFine.getSolverNumber() == pParent.getSolverNumber()) {
          switch (pFine.getType()) {
            case exahype::records::ADERDGCellDescription::Cell:
            case exahype::records::ADERDGCellDescription::Ancestor:
              assertion1(pFine.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,pFine.toString());
              subcellPosition = fineGridCell.computeSubcellPositionOfCellOrAncestor(pFine);
              restrictADERDGFaceData(pFine,subcellPosition.parentIndex,subcellPosition.subcellIndex);
              break;
            default:
              break;
          }
        }
      }
    }
    endpfor peano::datatraversal::autotuning::Oracle::getInstance()
        .parallelSectionHasTerminated(methodTrace);
  }

  logTraceOutWith1Argument("leaveCell(...)", fineGridCell);
}

void exahype::mappings::SpaceTimePredictor::restrictADERDGFaceData(
    const exahype::records::ADERDGCellDescription& cellDescription,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  // todo not dynamic with respect to the solver registry
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(parentIndex),
             cellDescription.toString());

  exahype::records::ADERDGCellDescription& cellDescriptionParent =
      ADERDGCellDescriptionHeap::getInstance().getData(parentIndex)[cellDescription.getSolverNumber()];

  assertion(cellDescriptionParent.getSolverNumber() == cellDescription.getSolverNumber());

  // Only do something if parent is an ancestor that holds data.
  if (cellDescriptionParent.getType() == exahype::records::ADERDGCellDescription::Ancestor) {
    #if defined(Debug)
    _parentOfCellOrAncestorFound += 1;
    #endif

    const int levelFine = cellDescription.getLevel();
    const int levelCoarse = cellDescriptionParent.getLevel();
    assertion(levelCoarse < levelFine);
    const int levelDelta = levelFine - levelCoarse;

    for (int d = 0; d < DIMENSIONS; d++) {
      // Check if cell is at "left" or "right" d face of parent
      if (subcellIndex[d] == 0) {
        const int faceIndex = 2 * d;

        exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
          exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);

        const int numberOfFaceDof = solver->getUnknownsPerFace();

        const double* lQhbndFine = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
            (faceIndex * numberOfFaceDof);
        const double* lFhbndFine =
            DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
            (faceIndex * numberOfFaceDof);
        double* lQhbndCoarse =DataHeap::getInstance().getData(cellDescriptionParent.getExtrapolatedPredictor()).data() +
            (faceIndex * numberOfFaceDof);
        double* lFhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getFluctuation()).data() +
            (faceIndex * numberOfFaceDof);

        tarch::multicore::Lock lock(_semaphoreForRestriction); // Is unlocked if lock gets out of scope.
        solver->faceUnknownsRestriction(lQhbndCoarse, lFhbndCoarse, lQhbndFine,
                                        lFhbndFine, levelCoarse, levelFine,
                                        getSubfaceIndex(subcellIndex, d));

      } else if (subcellIndex[d] == tarch::la::aPowI(levelDelta, 3) - 1) {
        const int faceIndex = 2 * d + 1;

        exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
          exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);

        const int numberOfFaceDof = solver->getUnknownsPerFace();

        const double* lQhbndFine = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
            (faceIndex * numberOfFaceDof);
        const double* lFhbndFine = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
            (faceIndex * numberOfFaceDof);
        double* lQhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getExtrapolatedPredictor()).data() +
            (faceIndex * numberOfFaceDof);
        double* lFhbndCoarse =
            DataHeap::getInstance().getData(cellDescriptionParent.getFluctuation()).data() +
            (faceIndex * numberOfFaceDof);

        tarch::multicore::Lock lock(_semaphoreForRestriction); // Is unlocked if lock gets out of scope.
        solver->faceUnknownsRestriction(lQhbndCoarse, lFhbndCoarse, lQhbndFine,
                                        lFhbndFine, levelCoarse, levelFine,
                                        getSubfaceIndex(subcellIndex, d));
      }
    }
  } else {
    #if defined(Debug)
    _parentOfCellOrAncestorNotFound += 1;
    #endif
  }
}

void exahype::mappings::SpaceTimePredictor::restrictFiniteVolumesSolution(
    const exahype::records::FiniteVolumesCellDescription& cellDescription,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  // todo not dynamic with respect to the solver registry
  assertion1(FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(parentIndex),
             cellDescription.toString());

  exahype::records::FiniteVolumesCellDescription& cellDescriptionParent =
      FiniteVolumesCellDescriptionHeap::getInstance().getData(parentIndex)[cellDescription.getSolverNumber()];

  assertion(cellDescriptionParent.getSolverNumber() == cellDescription.getSolverNumber());

  // todo @Dominic: implement!!!

  // Only do something if parent is an ancestor that holds data.
//  if (cellDescriptionParent.getType() == exahype::records::FiniteVolumesCellDescription::Ancestor) {
//    _parentOfCellOrAncestorFound += 1;
//
//    const int levelFine = cellDescription.getLevel();
//    const int levelCoarse = cellDescriptionParent.getLevel();
//    assertion(levelCoarse < levelFine);
//    const int levelDelta = levelFine - levelCoarse;
//
//    for (int d = 0; d < DIMENSIONS; d++) {
//      // Check if cell is at "left" or "right" d face of parent
//      if (subcellIndex[d] == 0) {
//        const int faceIndex = 2 * d;
//
//        exahype::solvers::FiniteVolumesSolver* solver = static_cast<exahype::solvers::FiniteVolumesSolver*>(
//          exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);
//
//        const int numberOfFaceDof = solver->getUnknownsPerFace();
//
//        const double* lQhbndFine = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
//            (faceIndex * numberOfFaceDof);
//        const double* lFhbndFine =
//            DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
//            (faceIndex * numberOfFaceDof);
//        double* lQhbndCoarse =DataHeap::getInstance().getData(cellDescriptionParent.getExtrapolatedPredictor()).data() +
//            (faceIndex * numberOfFaceDof);
//        double* lFhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getFluctuation()).data() +
//            (faceIndex * numberOfFaceDof);
//
//        solver->faceUnknownsRestriction(lQhbndCoarse, lFhbndCoarse, lQhbndFine,
//                                        lFhbndFine, levelCoarse, levelFine,
//                                        getSubfaceIndex(subcellIndex, d));
//
//      } else if (subcellIndex[d] == tarch::la::aPowI(levelDelta, 3) - 1) {
//        const int faceIndex = 2 * d + 1;
//
//        exahype::solvers::FiniteVolumesSolver* solver = static_cast<exahype::solvers::FiniteVolumesSolver*>(
//          exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);
//
//        const int numberOfFaceDof = solver->getUnknownsPerFace();
//
//        const double* lQhbndFine = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
//            (faceIndex * numberOfFaceDof);
//        const double* lFhbndFine = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
//            (faceIndex * numberOfFaceDof);
//        double* lQhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getExtrapolatedPredictor()).data() +
//            (faceIndex * numberOfFaceDof);
//        double* lFhbndCoarse =
//            DataHeap::getInstance().getData(cellDescriptionParent.getFluctuation()).data() +
//            (faceIndex * numberOfFaceDof);
//
//        solver->faceUnknownsRestriction(lQhbndCoarse, lFhbndCoarse, lQhbndFine,
//                                        lFhbndFine, levelCoarse, levelFine,
//                                        getSubfaceIndex(subcellIndex, d));
//      }
//    }
//  } else {
//    _parentOfCellOrAncestorNotFound += 1;
//  }
}

#ifdef Parallel
int exahype::mappings::SpaceTimePredictor::countListingsOfRemoteRankAtFace(
    const int faceIndex,
    exahype::Vertex* const verticesAroundCell,
    const peano::grid::VertexEnumerator& verticesEnumerator) {
  int result = 0;

  const int f = faceIndex % 2;   // "0" indicates a left face, "1" indicates a right face.
  const int d = (faceIndex-f)/2; // The normal direction: 0: x, 1: y, 1: z.

  tarch::la::Vector<DIMENSIONS,int> pos(1); // This is now the center, i.e., (1,1,...,1).
  pos(d) = 2*f;                             // This is a shift from the center by one unit in direction d.

  int faceNeighbourRank = -1; // This variable is introduced to make sure that the adjacent remote rank is unique.
  // TODO(Dominic): Uniqueness is probably guaranteed by the SFC based DD.
  dfor2(v) // Loop over vertices.
    if (verticesAroundCell[ verticesEnumerator(v) ].isAdjacentToRemoteRank()) {
      dfor2(a) // Loop over adjacent ranks. Does also include own rank.
        if (tarch::la::equals(v+a,pos) &&
            verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar]!=
            tarch::parallel::Node::getInstance().getRank()) {
          // Increment
          if (faceNeighbourRank==-1) {
            faceNeighbourRank = verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar];
          }
          if (verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar]==faceNeighbourRank) {
            result++;
          }
        }
      enddforx // a
    }
  enddforx // v

  assertion2(result==0||result==TWO_POWER_D_DIVIDED_BY_TWO,result,faceIndex);

  return result;
}

void exahype::mappings::SpaceTimePredictor::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // TODO(Dominic): Add to docu.
  //
  // Sender:
  // 1. Send out face data if cell type
  // is Cell/Ancestor/Descendant AND if the counter says yes
  // 2. Send out metadata: Send out "empty" metadata if the counters says no
  // or the src cellDescriptionsIndex is not valid.
  //
  // Receiver:
  // 1. Receive metadata
  // 2. If neighbour cell type is Cell/Ancestor/Descendant and the metadata is not "empty"
  // receive face data for cell types Cell/Ancestor/Descendant.
  // Drop messages for EmptyAncestor/EmptyDescendant and if the
  // destCellDescriptionsIndex is not valid.

// TODO(Dominic): This does prevent couting down the vertices.
//#if !defined(PeriodicBC)
//  if (vertex.isBoundary()) return;
//#endif

  // TODO(Dominic): Add to docu why we remove the vertex.isInside() constraint here.
  // We might also consider to remove it from the grid setup mapping functions.
  // fineGridCell.isInside does not imply that all adjacent vertices are
  // inside. If we count down the counter only on
  // fineGridVertices that are inside we might not send out all faces
  // of a cell that is close to the boundary.

  tarch::la::Vector<TWO_POWER_D, int>& adjacentADERDGCellDescriptionsIndices =
      vertex.getCellDescriptionsIndex();

  dfor2(dest)
    dfor2(src)
      if (vertex.hasToSendMetadata(_state,src,dest,toRank)) {
        // we are solely exchanging faces
        const int srcCellDescriptionIndex = adjacentADERDGCellDescriptionsIndices(srcScalar);

        if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(srcCellDescriptionIndex)) {
          assertion1(FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(srcCellDescriptionIndex),srcCellDescriptionIndex);

          decrementCounters(src,dest,srcCellDescriptionIndex);
          if (needToSendFaceData(src,dest,srcCellDescriptionIndex)) {
            sendADERDGFaceData(toRank,x,level,src,dest,srcCellDescriptionIndex,adjacentADERDGCellDescriptionsIndices(destScalar));
            //            sentFiniteVolumesFaceData(toRank,x,level,src,dest,srcCellDescriptionIndex,adjacentADERDGCellDescriptionsIndices(destScalar));
            auto encodedMetadata = exahype::Cell::encodeMetadata(srcCellDescriptionIndex);
            MetadataHeap::getInstance().sendData(
                encodedMetadata, toRank, x, level,
                peano::heap::MessageType::NeighbourCommunication);
          } else {
            logDebug("prepareSendToNeighbour(...)","[empty] sent to rank "<<toRank<<", x:"<<
                x.toString() << ", level=" <<level << ", vertex.adjacentRanks: "
                << vertex.getAdjacentRanks());
            auto encodedMetadata = exahype::Cell::createEncodedMetadataSequenceForInvalidCellDescriptionsIndex();
            MetadataHeap::getInstance().sendData(
                encodedMetadata, toRank, x, level,
                peano::heap::MessageType::NeighbourCommunication);
          }
        } else {
          logDebug("prepareSendToNeighbour(...)","[empty] sent to rank "<<toRank<<", x:"<<
              x.toString() << ", level=" <<level << ", vertex.adjacentRanks: "
              << vertex.getAdjacentRanks());
          auto encodedMetadata = exahype::Cell::createEncodedMetadataSequenceForInvalidCellDescriptionsIndex();
          MetadataHeap::getInstance().sendData(
              encodedMetadata, toRank, x, level,
              peano::heap::MessageType::NeighbourCommunication);
        }
      }
    enddforx
  enddforx
}

bool exahype::mappings::SpaceTimePredictor::needToSendFaceData(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    int cellDescriptionsIndex) {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!

  // ADER-DG
  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)) {
    if (p.getFaceDataExchangeCounter(faceIndex)!=0) {
      return false;
    }
  }

//  // TODO(Dominic): Introduce counters to FiniteVolumesCellDescription.
//  // FV
//  for (auto& p : FiniteVolumesCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)) {
//    if (p.getFaceDataExchangeCounter(faceIndex)!=0) {
//      return false;
//    }
//  }

  return true;
}

void exahype::mappings::SpaceTimePredictor::decrementCounters(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    int cellDescriptionsIndex) {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!

  // ADER-DG
  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)) {
    int newCounterValue = p.getFaceDataExchangeCounter(faceIndex)-1;
    assertion1(newCounterValue<TWO_POWER_D,newCounterValue);
    p.setFaceDataExchangeCounter(faceIndex,newCounterValue);
  }

//  // TODO(Dominic): Introduce counters to FiniteVolumesCellDescription.
//  // FV
//  for (auto& p : FiniteVolumesCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)) {
//    int newCounterValue = p.getFaceDataExchangeCounter(faceIndex)-1;
//    assertion1(newCounterValue<TWO_POWER_D,newCounterValue);
//    p.setFaceDataExchangeCounter(faceIndex,newCounterValue);
//  }
}

void exahype::mappings::SpaceTimePredictor::sendADERDGFaceData(
    int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    int level,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    int srcCellDescriptionIndex,
    int destCellDescriptionIndex) {
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!

  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(srcCellDescriptionIndex)) {
    exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
        exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);
    assertion1(p.getFaceDataExchangeCounter(faceIndex)==0,p.getFaceDataExchangeCounter(faceIndex));

    if (p.getType() == exahype::records::ADERDGCellDescription::Cell     ||
        p.getType() == exahype::records::ADERDGCellDescription::Ancestor ||
        p.getType() == exahype::records::ADERDGCellDescription::Descendant) {

      assertion(DataHeap::getInstance().isValidIndex(p.getExtrapolatedPredictor()));
      assertion(DataHeap::getInstance().isValidIndex(p.getFluctuation()));

      const int numberOfFaceDof = solver->getUnknownsPerFace();
      const double* lQhbnd = DataHeap::getInstance().getData(
          p.getExtrapolatedPredictor()).data() +
          (faceIndex * numberOfFaceDof);
      const double* lFhbnd = DataHeap::getInstance().getData(
          p.getFluctuation()).data() +
          (faceIndex * numberOfFaceDof);

      // TODO(Dominic) This can only happen if we consider
      // vertices that are outside in the MPI communication,
      // i.e. they belong to everyone else. Need to check if we
      // can drop isInside() in some of the mappings.
      // Q: If rank 0 is neighbour to every other rank that is on the boundary,
      // can we use rank 0 to handle the periodic boundary conditions?
      // If rank 0 is forwarding this will however take one extra iteration.
      // Can we resolve the periodic neighbour rank and alter the vertex position+-domain size(x,y,z) before/after we push
      // the metadata and facedata to the heap? The state knows all active ranks? Does he know
      // their respective subdomains?
      // We further need augmentation on the periodic neighbour rank. This can
      // be done by the same periodic neighbour rank resolving and exchanging metadata.
      // PeriodicBC might be easier to implement in a MPI scenario than
      // in a non-MPI scenario.
      if (destCellDescriptionIndex == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex) {
#ifdef PeriodicBC
        assertionMsg(false, "Vasco, we have to implement this");
        DataHeap::getInstance().sendData(
            sentMinMax, toRank, x, level, peano::heap::MessageType::NeighbourCommunication);
        DataHeap::getInstance().sendData(
            lQhbnd, numberOfFaceDof, toRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
        DataHeap::getInstance().sendData(
            lFhbnd, numberOfFaceDof, toRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
#else
        assertionMsg(false, "should never been entered");
#endif
      } else {
        logDebug(
            "sendADERDGFaceData(...)",
            "send three arrays to rank " <<
            toRank << " for vertex x=" << x << ", level=" << level <<
            ", dest type=" << multiscalelinkedcell::indexToString(destCellDescriptionIndex) <<
            ", src=" << src << ", dest=" << dest <<
            ", counter=" << p.getFaceDataExchangeCounter(faceIndex)
          );

        std::vector<double> sentMinMax(2);
        sentMinMax[0] = p.getSolutionMin(faceIndex);
        sentMinMax[1] = p.getSolutionMax(faceIndex);

        // Send order: minMax,lQhbnd,lFhbnd
        // Receive order: lFhbnd,lQhbnd,minMax
        DataHeap::getInstance().sendData(
            sentMinMax, toRank, x, level, peano::heap::MessageType::NeighbourCommunication);
        DataHeap::getInstance().sendData(
            lQhbnd, numberOfFaceDof, toRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
        DataHeap::getInstance().sendData(
            lFhbnd, numberOfFaceDof, toRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);

        // TODO(Dominic): If anarchic time stepping send the time step over too.
      }
    }

  }
}

bool exahype::mappings::SpaceTimePredictor::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  for (auto& p : exahype::solvers::RegisteredSolvers) {
    p->sendToRank(worker, _mpiTag);
  }

  return false;
}

void exahype::mappings::SpaceTimePredictor::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  for (auto& p : exahype::solvers::RegisteredSolvers) {
    p->receiveFromMasterRank(tarch::parallel::NodePool::getInstance().getMasterRank(),
                       _mpiTag);
  }
}



//
// Below all methods are nop.
//
// ====================================



void exahype::mappings::SpaceTimePredictor::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithMaster(
    const exahype::Cell& workerGridCell,
    exahype::Vertex* const workerGridVertices,
    const peano::grid::VertexEnumerator& workerEnumerator,
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker, const exahype::State& workerState,
    exahype::State& masterState) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif


exahype::mappings::SpaceTimePredictor::SpaceTimePredictor() {
  // do nothing
}

exahype::mappings::SpaceTimePredictor::~SpaceTimePredictor() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::SpaceTimePredictor::SpaceTimePredictor(
    const SpaceTimePredictor& masterThread) {
}
void exahype::mappings::SpaceTimePredictor::mergeWithWorkerThread(
    const SpaceTimePredictor& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::SpaceTimePredictor::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}


void exahype::mappings::SpaceTimePredictor::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::beginIteration(
    exahype::State& solverState) {
  #ifdef Parallel
  _state = &solverState;
  #endif

  #if defined(Debug)
  _parentOfCellOrAncestorNotFound = 0;
  _parentOfCellOrAncestorFound    = 0;
  _parentOfDescendantFound        = 0;
  #endif
}

void exahype::mappings::SpaceTimePredictor::endIteration(
    exahype::State& solverState) {
  logDebug("endIteration(...)", "_parentOfCellOrAncestorNotFound: "
      << _parentOfCellOrAncestorNotFound);
  logDebug("endIteration(...)",
      "_parentOfCellOrAncestorFound: " << _parentOfCellOrAncestorFound);
  logDebug("endIteration(...)",
      "_parentOfDescendantFound: " << _parentOfDescendantFound);
}

void exahype::mappings::SpaceTimePredictor::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
