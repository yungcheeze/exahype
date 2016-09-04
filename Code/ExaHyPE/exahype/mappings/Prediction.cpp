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
 
#include "exahype/mappings/Prediction.h"
#include "exahype/solvers/ADERDGSolver.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/utils/Globals.h"

#include "tarch/multicore/Loop.h"
#include "tarch/multicore/Lock.h"

#include "peano/utils/Loop.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include <algorithm>

peano::CommunicationSpecification
exahype::mappings::Prediction::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}

peano::MappingSpecification
exahype::mappings::Prediction::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::Prediction::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

// The remainder specs all are nop
peano::MappingSpecification
exahype::mappings::Prediction::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::Prediction::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::Prediction::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::Prediction::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::Prediction::_log(
    "exahype::mappings::Prediction");

#ifdef Parallel
int exahype::mappings::Prediction::_mpiTag =
    tarch::parallel::Node::reserveFreeTag(
        "exahype::mappings::Prediction");
#endif


tarch::multicore::BooleanSemaphore exahype::mappings::Prediction::_semaphoreForRestriction;

#if defined(Debug)
int exahype::mappings::Prediction::
_parentOfCellOrAncestorNotFound = 0;
int exahype::mappings::Prediction::
_parentOfCellOrAncestorFound    = 0;
int exahype::mappings::Prediction::
_parentOfDescendantFound        = 0;
#endif

exahype::mappings::Prediction::Prediction() {
  initTemporaryVariables();
}

exahype::mappings::Prediction::~Prediction() {
  deleteTemporaryVariables();
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::Prediction::Prediction(
    const Prediction& masterThread) {
  initTemporaryVariables();
}
void exahype::mappings::Prediction::mergeWithWorkerThread(
    const Prediction& workerThread) {
}
#endif

void exahype::mappings::Prediction::initTemporaryVariables() {
  assertion1(_lQi ==0,_lQi );
  assertion1(_lFi ==0,_lFi );
  assertion1(_lQhi==0,_lQhi);
  assertion1(_lFhi==0,_lFhi);

  int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _lQi  = new double*[numberOfSolvers];
  _lFi  = new double*[numberOfSolvers];
  _lQhi = new double*[numberOfSolvers];
  _lFhi = new double*[numberOfSolvers];

  int solverNumber=0;
  for (auto& solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::ADER_DG) {
      _lQi [solverNumber]  =
          new double[static_cast<exahype::solvers::ADERDGSolver*>(solver)->getSpaceTimeUnknownsPerCell()];
      _lFi [solverNumber]  =
          new double[static_cast<exahype::solvers::ADERDGSolver*>(solver)->getSpaceTimeFluxUnknownsPerCell()];
      _lQhi[solverNumber] =
          new double[static_cast<exahype::solvers::ADERDGSolver*>(solver)->getUnknownsPerCell()];
      _lFhi[solverNumber] =
          new double[static_cast<exahype::solvers::ADERDGSolver*>(solver)->getFluxUnknownsPerCell()];
    } else {
      _lQi[solverNumber]  = 0;
      _lFi[solverNumber]  = 0;
      _lQhi[solverNumber] = 0;
      _lFhi[solverNumber] = 0;
    }
    ++solverNumber;
  }
}

void exahype::mappings::Prediction::deleteTemporaryVariables() {
  int solverNumber=0;
  for (auto& solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::ADER_DG) {
      delete[]_lQi[solverNumber];
      delete[]_lFi[solverNumber];
      delete[]_lQhi[solverNumber];
      delete[]_lFhi[solverNumber];
    }
    _lQi [solverNumber] = 0;
    _lFi [solverNumber] = 0;
    _lQhi[solverNumber] = 0;
    _lFhi[solverNumber] = 0;

    ++solverNumber;
  }

  delete[] _lQi;
  delete[] _lFi;
  delete[]_lQhi;
  delete[]_lFhi;
}

void exahype::mappings::Prediction::beginIteration(
    exahype::State& solverState) {
  #if defined(Debug)
  _parentOfCellOrAncestorNotFound = 0;
  _parentOfCellOrAncestorFound    = 0;
  _parentOfDescendantFound        = 0;
  #endif
}

void exahype::mappings::Prediction::endIteration(
    exahype::State& solverState) {
  logDebug("endIteration(...)", "_parentOfCellOrAncestorNotFound: "
      << _parentOfCellOrAncestorNotFound);
  logDebug("endIteration(...)",
      "_parentOfCellOrAncestorFound: " << _parentOfCellOrAncestorFound);
  logDebug("endIteration(...)",
      "_parentOfDescendantFound: " << _parentOfDescendantFound);
}

void exahype::mappings::Prediction::enterCell(
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
        exahype::solvers::ADERDGSolver::Heap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).size());

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
      // For static AMR we assume always that the refinement event is None.

      // Update helper variables
      for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
        pFine.setRiemannSolvePerformed(faceIndex,false);
        // TODO(Dominic): Normally, isFaceInside has not to be called everytime here
        // but only once when the cell is initialised. Problem: Call addNewCellDescr.. from  merge..DueToForkOrJoin(...),
        // where no vertices are given. [Solved] - We send out the cellDescriptions from
        // the ADERDG/FV cell descr. heaps.
        pFine.setIsInside(faceIndex,isFaceInside(faceIndex,fineGridVertices,fineGridVerticesEnumerator));
        // TODO(Dominic): Add to docu.
        // 0 - no connection: no send. Set to unreachable value.
        // 2^{d-2} - full face connection where cell is inside but half of face vertices are outside:
        // send at time of 2^{d-2}-th touch of face.
        // 4^{d-2} - full face connection where cell is inside and face vertices are all inside:
        // send at time of 2^{d-2}-th touch of face.
        // We require that #ifdef vertex.isBoundary() is set.
        // These are counter values that have to be reinitialised every
        // iteration. They can stay here.
        #ifdef Parallel
        int listingsOfRemoteRank = countListingsOfRemoteRankByInsideVerticesAtFace(faceIndex,fineGridVertices,fineGridVerticesEnumerator);
        if (listingsOfRemoteRank==0) {
          listingsOfRemoteRank = TWO_POWER_D;
        }
        pFine.setFaceDataExchangeCounter(faceIndex,listingsOfRemoteRank);
        #endif
      }

      exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
        exahype::solvers::RegisteredSolvers[pFine.getSolverNumber()]);
      solver->synchroniseTimeStepping(pFine); // Time step synchr. might be done multiple times per traversal; but this is no issue.
      exahype::Cell::SubcellPosition subcellPosition;
      switch (pFine.getType()) {
      case exahype::records::ADERDGCellDescription::Cell:
        assertion1(pFine.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,pFine.toString());
        fineGridCell.validateNoNansInADERDGSolver(i,fineGridVerticesEnumerator,"exahype::mappings::Prediction::enterCell[pre]");
        performPredictionAndVolumeIntegral(pFine);
        fineGridCell.validateNoNansInADERDGSolver(i,fineGridVerticesEnumerator,"exahype::mappings::Prediction::enterCell[post]");
        break;
      case exahype::records::ADERDGCellDescription::Ancestor:
        assertion1(pFine.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,pFine.toString());
        prepareAncestor(pFine);
        break;
      case exahype::records::ADERDGCellDescription::Descendant:
        assertion1(pFine.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,pFine.toString());
        assertion1(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(pFine.getParentIndex()),pFine.toString());
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

bool exahype::mappings::Prediction::isFaceInside(
    const int faceIndex,
    exahype::Vertex* const verticesAroundCell,
    const peano::grid::VertexEnumerator& verticesEnumerator) {
  const int f = faceIndex % 2;   // "0" indicates a left face, "1" indicates a right face.
  const int d = (faceIndex-f)/2; // The normal direction: 0: x, 1: y, 1: z.

  dfor2(v) // Loop over vertices.
    if (v(d) == f && verticesAroundCell[ verticesEnumerator(v) ].isInside()) {
      return true;
    }
  enddforx // v
  return false;
}

void exahype::mappings::Prediction::performPredictionAndVolumeIntegral(
    exahype::records::ADERDGCellDescription& p) {
  exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
              exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);

  assertion1(p.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getSolution()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getUpdate()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getExtrapolatedPredictor()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getFluctuation()),p.toString());

  // persistent fields
  // volume DoF (basisSize**(DIMENSIONS))
  double* luh  = DataHeap::getInstance().getData(p.getSolution()).data();
  double* lduh = DataHeap::getInstance().getData(p.getUpdate()).data();
  // face DoF (basisSize**(DIMENSIONS-1))
  double* lQhbnd = DataHeap::getInstance().getData(p.getExtrapolatedPredictor()).data();
  double* lFhbnd = DataHeap::getInstance().getData(p.getFluctuation()).data();

  solver->spaceTimePredictor(
      _lQi[p.getSolverNumber()],
      _lFi[p.getSolverNumber()],
      _lQhi[p.getSolverNumber()],
      _lFhi[p.getSolverNumber()],
      lQhbnd,
      lFhbnd,
      luh, p.getSize(),
      p.getPredictorTimeStepSize());

  solver->volumeIntegral(
      lduh,
      _lFhi[p.getSolverNumber()],
      p.getSize());

    for (int i=0; i<solver->getSpaceTimeUnknownsPerCell(); i++) {
      assertion3(std::isfinite(_lQi[p.getSolverNumber()][i]),p.toString(),"performPredictionAndVolumeIntegral(...)",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.

    for (int i=0; i<solver->getSpaceTimeFluxUnknownsPerCell(); i++) {
      assertion3(std::isfinite(_lFi[p.getSolverNumber()][i]), p.toString(),"performPredictionAndVolumeIntegral",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
      for (int i=0; i<solver->getUnknownsPerCell(); i++) {
      assertion3(std::isfinite(_lQhi[p.getSolverNumber()][i]),p.toString(),"performPredictionAndVolumeIntegral(...)",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.

    for (int i=0; i<solver->getFluxUnknownsPerCell(); i++) {
      assertion3(std::isfinite(_lFhi[p.getSolverNumber()][i]),p.toString(),"performPredictionAndVolumeIntegral(...)",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
}

void exahype::mappings::Prediction::prepareAncestor(
    exahype::records::ADERDGCellDescription& p) {
  exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
      exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);

  assertion1(p.getType()==exahype::records::ADERDGCellDescription::Ancestor,p.toString());
  std::fill_n(DataHeap::getInstance().getData(p.getExtrapolatedPredictor()).begin(),
      solver->getUnknownsPerCellBoundary(), 0.0);
  std::fill_n(DataHeap::getInstance().getData(p.getFluctuation()).begin(),
      solver->getUnknownsPerCellBoundary(), 0.0);
}

void exahype::mappings::Prediction::prolongateADERDGFaceData(
    const exahype::records::ADERDGCellDescription& cellDescription,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  // todo not dynamic with respect to the solver registry
  assertion1(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(parentIndex),
             cellDescription.toString());

  exahype::records::ADERDGCellDescription& cellDescriptionParent =
      exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
          parentIndex)[cellDescription.getSolverNumber()];

#if defined(Debug)
  _parentOfDescendantFound += 1;
#endif

  assertion(cellDescriptionParent.getSolverNumber() == cellDescription.getSolverNumber());
  assertion(cellDescriptionParent.getType() == exahype::records::ADERDGCellDescription::Cell ||
            cellDescriptionParent.getType() == exahype::records::ADERDGCellDescription::Descendant);

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

void exahype::mappings::Prediction::prolongateFiniteVolumesFaceData(
    const exahype::records::FiniteVolumesCellDescription& cellDescription,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  // todo not dynamic with respect to the solver registry
  assertion1(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(parentIndex),
             cellDescription.toString());

  exahype::records::FiniteVolumesCellDescription& cellDescriptionParent =
      exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(
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
exahype::mappings::Prediction::getSubfaceIndex(
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

void exahype::mappings::Prediction::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("leaveCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(
          fineGridCell.getCellDescriptionsIndex())) {
    const int numberOfADERDGCellDescriptions = static_cast<int>(
        exahype::solvers::ADERDGSolver::Heap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).size());
    // please use a different UserDefined per mapping/event
    const peano::datatraversal::autotuning::MethodTrace methodTrace = peano::datatraversal::autotuning::UserDefined10;
    const int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfADERDGCellDescriptions, methodTrace);
    pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
        records::ADERDGCellDescription& pFine = fineGridCell.getADERDGCellDescription(i);

    // if we have at least one parent
    if (exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(pFine.getParentIndex())) {
      for (auto& pParent : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(pFine.getParentIndex())) {
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

void exahype::mappings::Prediction::restrictADERDGFaceData(
    const exahype::records::ADERDGCellDescription& cellDescription,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  // todo not dynamic with respect to the solver registry
  assertion1(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(parentIndex),
             cellDescription.toString());

  exahype::records::ADERDGCellDescription& cellDescriptionParent =
      exahype::solvers::ADERDGSolver::Heap::getInstance().getData(parentIndex)[cellDescription.getSolverNumber()];

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

void exahype::mappings::Prediction::restrictFiniteVolumesSolution(
    const exahype::records::FiniteVolumesCellDescription& cellDescription,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  // todo not dynamic with respect to the solver registry
  assertion1(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(parentIndex),
             cellDescription.toString());

  exahype::records::FiniteVolumesCellDescription& cellDescriptionParent =
      exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(parentIndex)[cellDescription.getSolverNumber()];

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
int exahype::mappings::Prediction::countListingsOfRemoteRankByInsideVerticesAtFace(
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
        if (verticesAroundCell[ verticesEnumerator(v) ].isInside() &&
            tarch::la::equals(v+a,pos) &&
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

  // result must be either no connection, edge connection, or whole face connection.
  assertion2(result==0||result==TWO_POWER_D_DIVIDED_BY_TWO/2||result==TWO_POWER_D_DIVIDED_BY_TWO,result,faceIndex);

  return result;
}


void exahype::mappings::Prediction::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
#if !defined(PeriodicBC)
  if (vertex.isBoundary()) return;
#endif

  dfor2(dest)
    dfor2(src)
    if (vertex.hasToSendMetadata(src,dest,toRank)) {
      vertex.tryDecrementFaceDataExchangeCountersOfSource(src,dest);
      if (vertex.hasToSendDataToNeighbour(src,dest)) {
        sendSolverData(
            toRank,src,dest,
            vertex.getCellDescriptionsIndex()[srcScalar],
            vertex.getCellDescriptionsIndex()[destScalar],
            x,level);
      } else {
        sendEmptySolverData(toRank,src,dest,x,level);
      }
    }
    enddforx
  enddforx
}

void exahype::mappings::Prediction::sendEmptySolverData(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  for (exahype::solvers::Solver* solver : exahype::solvers::RegisteredSolvers) {
    solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
  }

  auto encodedMetadata = exahype::Vertex::createEncodedMetadataSequenceWithInvalidEntries();
  MetadataHeap::getInstance().sendData(
      encodedMetadata, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}

void exahype::mappings::Prediction::sendSolverData(
    const int                                    toRank,
    const tarch::la::Vector<DIMENSIONS,int>&     src,
    const tarch::la::Vector<DIMENSIONS,int>&     dest,
    const int                                    srcCellDescriptionIndex,
    const int                                    destCellDescriptionIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  assertion(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionIndex));
  assertion(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(srcCellDescriptionIndex));

  int solverNumber=0;
  for (exahype::solvers::Solver* solver : exahype::solvers::RegisteredSolvers) {
    bool solverIsRegistered = false;
    int element = 0;
    for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(srcCellDescriptionIndex)) {
      if (p.getSolverNumber()==solverNumber) {
        assertionMsg(!solverIsRegistered,"The solverNumber is a unique key. Every solver is allowed to "
            "be registered only once on each face description.");
        solverIsRegistered=true;
        solver->sendDataToNeighbour(toRank,srcCellDescriptionIndex,element,src,dest,x,level);
      }
      ++element;
    }

    if (!solverIsRegistered)
      solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);

    ++solverNumber;
  }

  auto encodedMetadata = exahype::Vertex::encodeMetadata(srcCellDescriptionIndex);
  MetadataHeap::getInstance().sendData(
      encodedMetadata, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}

bool exahype::mappings::Prediction::prepareSendToWorker(
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

void exahype::mappings::Prediction::receiveDataFromMaster(
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



void exahype::mappings::Prediction::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::Prediction::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Prediction::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Prediction::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithMaster(
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

void exahype::mappings::Prediction::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::Prediction::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Prediction::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}


void exahype::mappings::Prediction::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::Prediction::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
