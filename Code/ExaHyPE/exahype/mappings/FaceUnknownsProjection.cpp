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
 
#include "exahype/mappings/FaceUnknownsProjection.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"
#include "kernels/KernelCalls.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "tarch/multicore/Lock.h"

#include "exahype/VertexOperations.h"

#include <cstring>


peano::CommunicationSpecification
exahype::mappings::FaceUnknownsProjection::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}

int exahype::mappings::FaceUnknownsProjection::_parentOfCellOrAncestorNotFound =
    0;
int exahype::mappings::FaceUnknownsProjection::_parentOfCellOrAncestorFound = 0;
int exahype::mappings::FaceUnknownsProjection::_parentOfDescendantFound = 0;

peano::MappingSpecification
exahype::mappings::FaceUnknownsProjection::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::FaceUnknownsProjection::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::FaceUnknownsProjection::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::FaceUnknownsProjection::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::FaceUnknownsProjection::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::FaceUnknownsProjection::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::FaceUnknownsProjection::_log(
    "exahype::mappings::FaceUnknownsProjection");

tarch::multicore::BooleanSemaphore exahype::mappings::FaceUnknownsProjection::_semaphoreForRestriction;

void exahype::mappings::FaceUnknownsProjection::enterCell(
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
    const int numberOfADERDGCellDescriptions = fineGridCell.getNumberOfADERDGCellDescriptions();
    // please use a different UserDefined per mapping/event
    const peano::datatraversal::autotuning::MethodTrace methodTrace = peano::datatraversal::autotuning::UserDefined1;
    const int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfADERDGCellDescriptions, methodTrace);
    pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
      records::ADERDGCellDescription& pFine = fineGridCell.getADERDGCellDescription(i);
      exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
        exahype::solvers::RegisteredSolvers[pFine.getSolverNumber()]);
      switch (pFine.getType()) {
        case exahype::records::ADERDGCellDescription::Ancestor:
          std::memset(DataHeap::getInstance().getData(pFine.getExtrapolatedPredictor()).data(), 0,
                 sizeof(double) * solver->getUnknownsPerCellBoundary());
          std::memset(DataHeap::getInstance().getData(pFine.getFluctuation()).data(), 0,
                 sizeof(double) * solver->getUnknownsPerCellBoundary());
          break;
        default:
          break;
      }

      // if we have at least one parent
      if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(pFine.getParentIndex())) {
        for (auto& pParent : ADERDGCellDescriptionHeap::getInstance().getData(pFine.getParentIndex())) {
          exahype::Cell::SubcellPosition subcellPosition;

          if (pFine.getSolverNumber() == pParent.getSolverNumber()) {
            switch (pFine.getType()) {
              case exahype::records::ADERDGCellDescription::Descendant:
                assertion1(pFine.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,pFine.toString());
                subcellPosition = fineGridCell.computeSubcellPositionOfDescendant(pFine);

                prolongateADERDGFaceData(pFine, subcellPosition.parentIndex,
                                   subcellPosition.subcellIndex);

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
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::FaceUnknownsProjection::prolongateADERDGFaceData(
    const exahype::records::ADERDGCellDescription& cellDescription,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) const {
  // todo not dynamic with respect to the solver registry
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(parentIndex),
             cellDescription.toString());

  exahype::records::ADERDGCellDescription& cellDescriptionParent =
      ADERDGCellDescriptionHeap::getInstance().getData(
          parentIndex)[cellDescription.getSolverNumber()];

  _parentOfDescendantFound += 1;

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

void exahype::mappings::FaceUnknownsProjection::prolongateFiniteVolumesFaceData(
    const exahype::records::FiniteVolumesCellDescription& cellDescription,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) const {
  // todo not dynamic with respect to the solver registry
  assertion1(FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(parentIndex),
             cellDescription.toString());

  exahype::records::FiniteVolumesCellDescription& cellDescriptionParent =
      FiniteVolumesCellDescriptionHeap::getInstance().getData(
          parentIndex)[cellDescription.getSolverNumber()];

  _parentOfDescendantFound += 1;

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
exahype::mappings::FaceUnknownsProjection::getSubfaceIndex(
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex, const int d) const {
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

void exahype::mappings::FaceUnknownsProjection::leaveCell(
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

void exahype::mappings::FaceUnknownsProjection::restrictADERDGFaceData(
    const exahype::records::ADERDGCellDescription& cellDescription,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) const {
  // todo not dynamic with respect to the solver registry
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(parentIndex),
             cellDescription.toString());

  exahype::records::ADERDGCellDescription& cellDescriptionParent =
      ADERDGCellDescriptionHeap::getInstance().getData(parentIndex)[cellDescription.getSolverNumber()];

  assertion(cellDescriptionParent.getSolverNumber() == cellDescription.getSolverNumber());

  // Only do something if parent is an ancestor that holds data.
  if (cellDescriptionParent.getType() == exahype::records::ADERDGCellDescription::Ancestor) {
    _parentOfCellOrAncestorFound += 1;

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
    _parentOfCellOrAncestorNotFound += 1;
  }
}

void exahype::mappings::FaceUnknownsProjection::restrictFiniteVolumesSolution(
    const exahype::records::FiniteVolumesCellDescription& cellDescription,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) const {
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

void exahype::mappings::FaceUnknownsProjection::beginIteration(
    exahype::State& solverState) {
  _parentOfCellOrAncestorNotFound = 0;
  _parentOfCellOrAncestorFound    = 0;
  _parentOfDescendantFound        = 0;
}

void exahype::mappings::FaceUnknownsProjection::endIteration(
    exahype::State& solverState) {
  logDebug("endIteration(...)", "_parentOfCellOrAncestorNotFound: "
                                    << _parentOfCellOrAncestorNotFound);
  logDebug("endIteration(...)",
           "_parentOfCellOrAncestorFound: " << _parentOfCellOrAncestorFound);
  logDebug("endIteration(...)",
           "_parentOfDescendantFound: " << _parentOfDescendantFound);
}

exahype::mappings::FaceUnknownsProjection::FaceUnknownsProjection() {
  // do nothing
}

exahype::mappings::FaceUnknownsProjection::~FaceUnknownsProjection() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::FaceUnknownsProjection::FaceUnknownsProjection(
    const FaceUnknownsProjection& masterThread) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::mergeWithWorkerThread(
    const FaceUnknownsProjection& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::FaceUnknownsProjection::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::FaceUnknownsProjection::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::FaceUnknownsProjection::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  // do nothing
  return false;
}

void exahype::mappings::FaceUnknownsProjection::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::mergeWithMaster(
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

void exahype::mappings::FaceUnknownsProjection::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::FaceUnknownsProjection::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::FaceUnknownsProjection::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
