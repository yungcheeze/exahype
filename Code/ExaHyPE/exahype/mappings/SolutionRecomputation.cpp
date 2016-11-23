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
 
#include "exahype/mappings/SolutionRecomputation.h"

#include "tarch/multicore/Loop.h"

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"
#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/VertexOperations.h"
#include "multiscalelinkedcell/HangingVertexBookkeeper.h"


peano::CommunicationSpecification
exahype::mappings::SolutionRecomputation::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}

peano::MappingSpecification
exahype::mappings::SolutionRecomputation::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces);
}

// Below specs are all nop

peano::MappingSpecification
exahype::mappings::SolutionRecomputation::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

peano::MappingSpecification
exahype::mappings::SolutionRecomputation::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification
exahype::mappings::SolutionRecomputation::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::SolutionRecomputation::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}


peano::MappingSpecification
exahype::mappings::SolutionRecomputation::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::SolutionRecomputation::_log(
    "exahype::mappings::SolutionRecomputation");

void exahype::mappings::SolutionRecomputation::prepareTemporaryVariables() {
  assertion(_tempStateSizedVectors       ==nullptr);
  assertion(_tempStateSizedSquareMatrices==nullptr);
  assertion(_tempFaceUnknowns            ==nullptr);
  assertion(_tempUnknowns                ==nullptr);

  int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _tempStateSizedVectors        = new double**[numberOfSolvers];
  _tempStateSizedSquareMatrices = new double**[numberOfSolvers];
  _tempFaceUnknowns             = new double**[numberOfSolvers];
  _tempUnknowns                 = new double**[numberOfSolvers];
//    _tempSpaceTimeFaceUnknownsArray = new double* [numberOfSolvers]; todo

  int solverNumber=0;
  for (auto solver : exahype::solvers::RegisteredSolvers) {
    int numberOfStateSizedVectors  = 0; // TOOD(Dominic): Check if we need number of parameters too
    int numberOfStateSizedMatrices = 0;
    int numberOfFaceUnknowns       = 0;
    int lengthOfFaceUnknowns       = 0;
    int numberOfUnknowns           = 0;
    int lengthOfUnknowns           = 0;
    switch (solver->getType()) {
      case exahype::solvers::Solver::Type::ADERDG:
        numberOfStateSizedVectors  = 5; // See riemannSolverLinear
        numberOfStateSizedMatrices = 3; // See riemannSolverLinear
        numberOfFaceUnknowns       = 3; // See exahype::solvers::ADERDGSolver::applyBoundaryConditions
        lengthOfFaceUnknowns       =
            static_cast<exahype::solvers::ADERDGSolver*>(solver)->getUnknownsPerFace();
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        // Needs the same temporary variables as the normal ADER-DG scheme
        // plus the ones for the Finite Volume scheme.
        numberOfStateSizedVectors  = 5;
        numberOfStateSizedMatrices = 3;
        numberOfFaceUnknowns       = 3;
        lengthOfFaceUnknowns       = std::max(
                    static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->_solver->getUnknownsPerFace(),
                    static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->_limiter->getUnknownsPerFace() );
        numberOfUnknowns     = 0; // TODO(Dominic): We do not consider high-order FV methods yet; numberOfUnknowns is thus set to zero.
        lengthOfUnknowns     = 0;
        break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        numberOfUnknowns     = 0; // TODO(Dominic): We do not consider high-order FV methods yet; numberOfUnknowns is thus set to zero.
        lengthOfUnknowns     = 0;
        numberOfFaceUnknowns = 2; // See exahype::solvers::FiniteVolumesSolver::mergeWithBoundaryData
        lengthOfFaceUnknowns =
            static_cast<exahype::solvers::FiniteVolumesSolver*>(solver)->getUnknownsPerFace();
        numberOfStateSizedVectors = 5;
        break;
    }

    _tempStateSizedVectors[solverNumber] = nullptr;
    if (numberOfStateSizedVectors>0) {
      _tempStateSizedVectors[solverNumber] = new double*[numberOfStateSizedVectors];
      _tempStateSizedVectors[solverNumber][0] = new double[numberOfStateSizedVectors*solver->getNumberOfVariables()];
      for (int i=1; i<numberOfStateSizedVectors; ++i) { // see riemanSolverLinear
        _tempStateSizedVectors[solverNumber][i] = _tempStateSizedVectors[solverNumber][i-1] + solver->getNumberOfVariables();
      }
    }
    //
    _tempStateSizedSquareMatrices[solverNumber] = nullptr;
    if (numberOfStateSizedMatrices>0) {
      _tempStateSizedSquareMatrices[solverNumber] = new double*[numberOfStateSizedMatrices];
      _tempStateSizedSquareMatrices[solverNumber][0] =
          new double[numberOfStateSizedMatrices* solver->getNumberOfVariables() * solver->getNumberOfVariables()];
      for (int i=1; i<numberOfStateSizedMatrices; ++i) { // see riemanSolverLinear
        _tempStateSizedSquareMatrices[solverNumber][i] =
            _tempStateSizedSquareMatrices[solverNumber][i-1] +
            solver->getNumberOfVariables() * solver->getNumberOfVariables();
      }
    }
    //
    _tempFaceUnknowns[solverNumber] = nullptr;
    if (numberOfFaceUnknowns>0) {
      _tempFaceUnknowns[solverNumber] = new double*[numberOfFaceUnknowns];
      _tempFaceUnknowns[solverNumber][0] = new double[numberOfFaceUnknowns*lengthOfFaceUnknowns];
      for (int i=1; i<numberOfFaceUnknowns; ++i) { // see ADERDGSolver::applyBoundaryConditions(...)
        _tempFaceUnknowns[solverNumber][i] = _tempFaceUnknowns[solverNumber][i-1] + lengthOfFaceUnknowns;
      }
    }
    //
    _tempUnknowns[solverNumber] = nullptr;
    if (numberOfUnknowns>0) {
      _tempUnknowns[solverNumber] = new double*[numberOfUnknowns];
      _tempUnknowns[solverNumber][0] = new double[numberOfUnknowns*lengthOfUnknowns];
      for (int i=1; i<numberOfUnknowns; ++i) { // see ADERDGSolver::applyBoundaryConditions(...)
        _tempUnknowns[solverNumber][i] = _tempUnknowns[solverNumber][i-1] + lengthOfUnknowns;
      }
    }

    ++solverNumber;
  }
}

void exahype::mappings::SolutionRecomputation::deleteTemporaryVariables() {
  if (_tempStateSizedVectors!=nullptr) {
    assertion(_tempStateSizedSquareMatrices!=nullptr);

    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      int numberOfStateSizedVectors  = 0;
      int numberOfStateSizedMatrices = 0;
      int numberOfFaceUnknowns       = 0;
      int numberOfUnknowns           = 0;
      switch (solver->getType()) {
        case exahype::solvers::Solver::Type::ADERDG:
        case exahype::solvers::Solver::Type::LimitingADERDG:
          numberOfStateSizedVectors  = 5; // See riemannSolverLinear
          numberOfStateSizedMatrices = 3; // See riemannSolverLinear
          numberOfFaceUnknowns       = 3; // See exahype::solvers::ADERDGSolver::applyBoundaryConditions
          break;
        case exahype::solvers::Solver::Type::FiniteVolumes:
          // TODO(Dominic): We do not consider high-order FV methods yet;
          // numberOfUnknowns is thus set to zero.
          numberOfUnknowns = 0;
          numberOfFaceUnknowns = 2; // See exahype::solvers::FiniteVolumesSolver::mergeWithBoundaryData
          numberOfStateSizedVectors = 5;
          break;
      }

      if (numberOfStateSizedVectors>0) {
        delete[] _tempStateSizedVectors[solverNumber][0];
        delete[] _tempStateSizedVectors[solverNumber];
        _tempStateSizedVectors[solverNumber] = nullptr;
      }
      //
      if (numberOfStateSizedMatrices>0) {
        delete[] _tempStateSizedSquareMatrices[solverNumber][0];
        delete[] _tempStateSizedSquareMatrices[solverNumber];
        _tempStateSizedSquareMatrices[solverNumber] = nullptr;
      }
      //
      if (numberOfFaceUnknowns>0) {
        delete[] _tempFaceUnknowns[solverNumber][0];
        delete[] _tempFaceUnknowns[solverNumber];
        _tempFaceUnknowns[solverNumber] = nullptr;
      }
      //
      if (numberOfUnknowns>0) {
        delete[] _tempUnknowns[solverNumber][0];
        delete[] _tempUnknowns[solverNumber];
        _tempUnknowns[solverNumber] = nullptr;
      }
      //
      // _tempSpaceTimeFaceUnknownsArray[solverNumber] = nullptr; // todo

      ++solverNumber;
    }

    delete[] _tempStateSizedVectors;
    delete[] _tempStateSizedSquareMatrices;
    delete[] _tempFaceUnknowns;
    delete[] _tempUnknowns;
//    delete[] _tempSpaceTimeFaceUnknownsArray; todo
    _tempStateSizedVectors        = nullptr;
    _tempStateSizedSquareMatrices = nullptr;
    _tempFaceUnknowns             = nullptr;
    _tempUnknowns                 = nullptr;
//    _tempSpaceTimeFaceUnknownsArray  = nullptr; todo
  }
}

exahype::mappings::SolutionRecomputation::SolutionRecomputation()
 #ifdef Debug
 :
 _interiorFaceMerges(0),
 _boundaryFaceMerges(0)
 #endif
{
  // do nothing
}

exahype::mappings::SolutionRecomputation::~SolutionRecomputation() {
  deleteTemporaryVariables();
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::SolutionRecomputation::SolutionRecomputation(
    const SolutionRecomputation& masterThread)
: _localState(masterThread._localState),
  _tempStateSizedVectors(nullptr),
  _tempUnknowns(nullptr) {
  prepareTemporaryVariables();
}

void exahype::mappings::SolutionRecomputation::mergeWithWorkerThread(
    const SolutionRecomputation& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::SolutionRecomputation::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  _localState = solverState;

  prepareTemporaryVariables();

  #ifdef Debug // TODO(Dominic): And not parallel and not shared memory
  _interiorFaceMerges = 0;
  _boundaryFaceMerges = 0;
  #endif

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::SolutionRecomputation::endIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("endIteration(State)", solverState);

  deleteTemporaryVariables();

  #if defined(Debug) // TODO(Dominic): Use logDebug if it works with filters
  logInfo("endIteration(...)","interiorFaceSolves: " << _interiorFaceMerges);
  logInfo("endIteration(...)","boundaryFaceSolves: " << _boundaryFaceMerges);
  #endif

  logTraceOutWith1Argument("endIteration(State)", solverState);
}

void exahype::mappings::SolutionRecomputation::enterCell(
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
    const int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
    auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfSolvers, peano::datatraversal::autotuning::MethodTrace::UserDefined14);
    pfor(i, 0, numberOfSolvers, grainSize.getGrainSize())
      auto solver = exahype::solvers::RegisteredSolvers[i];

      const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),i);
      if (element!=exahype::solvers::Solver::NotFound) {
        if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
            && static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->limiterDomainHasChanged()) {
          auto* limitingADERSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);

          // TODO(Dominic): Solver time step data merging is happening before.
          //  TODO(Dominic): remove
          if (fineGridCell.getCellDescriptionsIndex()==78) {
            logInfo("enterCell(...)","[pre] solverPatch="<<exahype::solvers::ADERDGSolver::Heap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[0].toString());

            // TODO(Dominic): TimeStamp is smaller than zero already here. Why?

            if (exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).size()>0) {
              logInfo("enterCell(...)","[pre] limiterPatch="<<exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[0].toString());
            }

            logInfo("enterCell(...)","[pre] solver="<<static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->toString());
            logInfo("enterCell(...)","[pre] limiter="<<static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiter()->toString());
          }


          limitingADERSolver->updateMergedLimiterStatus(fineGridCell.getCellDescriptionsIndex(),element); // update before recomputation

          limitingADERSolver->recomputeSolution(
              fineGridCell.getCellDescriptionsIndex(),
              element,
              _tempStateSizedVectors[i],
              _tempUnknowns[i],
              fineGridVertices,
              fineGridVerticesEnumerator);

          // It is important that we update the limiter status only after the recomputation since we use
          // the previous and current limiter status in the recomputation.
          limitingADERSolver->updateLimiterStatus(fineGridCell.getCellDescriptionsIndex(),element);

          static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
                        determineMinAndMax(fineGridCell.getCellDescriptionsIndex(),element);

          //  TODO(Dominic): remove
          if (fineGridCell.getCellDescriptionsIndex()==78) {
            logInfo("enterCell(...)","[post] solverPatch="<<exahype::solvers::ADERDGSolver::Heap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[0].toString());

            // TODO(Dominic): TimeStamp is smaller than zero already here. Why?

            if (exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).size()>0) {
              logInfo("enterCell(...)","[post] limiterPatch="<<exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[0].toString());
            }

            logInfo("enterCell(...)","[post] solver="<<static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->toString());
            logInfo("enterCell(...)","[post] limiter="<<static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiter()->toString());
          }
        }

        solver->prepareNextNeighbourMerging(
            fineGridCell.getCellDescriptionsIndex(),element,
            fineGridVertices,fineGridVerticesEnumerator); // !!! Has to be done for all solvers (cf. touchVertexFirstTime etc.)
      }
    endpfor
    grainSize.parallelSectionHasTerminated();
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::SolutionRecomputation::touchVertexFirstTime(
  exahype::Vertex& fineGridVertex,
  const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
  const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
  exahype::Vertex* const coarseGridVertices,
  const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
  exahype::Cell& coarseGridCell,
  const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
    dfor2(pos1)
      dfor2(pos2)
        if (fineGridVertex.hasToMergeNeighbours(pos1,pos2)) { // Assumes that we have to valid indices
          auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().
              parallelise(solvers::RegisteredSolvers.size(), peano::datatraversal::autotuning::MethodTrace::UserDefined15);
          pfor(solverNumber, 0, static_cast<int>(solvers::RegisteredSolvers.size()),grainSize.getGrainSize())
            auto solver = exahype::solvers::RegisteredSolvers[solverNumber];

            if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
                && static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->_limiterDomainHasChanged) {
              const int cellDescriptionsIndex1 = fineGridVertex.getCellDescriptionsIndex()[pos1Scalar];
              const int cellDescriptionsIndex2 = fineGridVertex.getCellDescriptionsIndex()[pos2Scalar];
              const int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
              const int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
              if (element2>=0 && element1>=0) {
                auto& solverPatch1 = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
                    getSolver()->getCellDescription(cellDescriptionsIndex1,element1);
                auto& solverPatch2 = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
                    getSolver()->getCellDescription(cellDescriptionsIndex2,element2);

                static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
                    mergeNeighboursBasedOnLimiterStatus(
                    cellDescriptionsIndex1,element1,
                    cellDescriptionsIndex2,element2,
                    solverPatch1.getMergedLimiterStatus(0),
                    solverPatch2.getMergedLimiterStatus(0), // !!! We assume here that we have already unified the merged limiter status values
                    pos1,pos2,
                    _tempFaceUnknowns[solverNumber],
                    _tempStateSizedVectors[solverNumber],
                    _tempStateSizedSquareMatrices[solverNumber]);
              }
            }

            #ifdef Debug // TODO(Dominic)
            _interiorFaceMerges++;
            #endif
          endpfor
          grainSize.parallelSectionHasTerminated();

          fineGridVertex.setMergePerformed(pos1,pos2,true);
        }
        if (fineGridVertex.hasToMergeWithBoundaryData(pos1,pos2)) {
          auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().
              parallelise(solvers::RegisteredSolvers.size(), peano::datatraversal::autotuning::MethodTrace::UserDefined16);
          pfor(solverNumber, 0, static_cast<int>(solvers::RegisteredSolvers.size()),grainSize.getGrainSize())
            auto solver = exahype::solvers::RegisteredSolvers[solverNumber];
            const int cellDescriptionsIndex1 = fineGridVertex.getCellDescriptionsIndex()[pos1Scalar];
            const int cellDescriptionsIndex2 = fineGridVertex.getCellDescriptionsIndex()[pos2Scalar];
            int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
            int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
            assertion4((element1==exahype::solvers::Solver::NotFound &&
                        element2==exahype::solvers::Solver::NotFound)
                       || (element1 >= 0 && element2==exahype::solvers::Solver::NotFound)
                       || (element2 >= 0 && element1==exahype::solvers::Solver::NotFound),
                       cellDescriptionsIndex1,cellDescriptionsIndex2,element1,element2);

            if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
                && static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->_limiterDomainHasChanged) {
              if (element1 >= 0) {
                auto& solverPatch1 = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
                    getSolver()->getCellDescription(cellDescriptionsIndex1,element1);

                static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
                    mergeWithBoundaryDataBasedOnLimiterStatus( // !!! Be aware of indices "2" and "1" and the order of the arguments.
                                              cellDescriptionsIndex1,element1,
                                              solverPatch1.getMergedLimiterStatus(0), // !!! We assume here that we have already unified the merged limiter status values.
                                              pos1,pos2,                              // The cell-based limiter status is still holding the old value though.
                                              _tempFaceUnknowns[solverNumber],
                                              _tempStateSizedVectors[solverNumber],
                                              _tempStateSizedSquareMatrices[solverNumber]);

                #ifdef Debug
                _boundaryFaceMerges++;
                #endif
              }
              if (element2 >= 0){
                auto& solverPatch2 = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
                    getSolver()->getCellDescription(cellDescriptionsIndex2,element2);

                static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
                    mergeWithBoundaryDataBasedOnLimiterStatus( // !!! Be aware of indices "2" and "1" and the order of the arguments.
                                              cellDescriptionsIndex2,element2,
                                              solverPatch2.getMergedLimiterStatus(0), // !!! We assume here that we have already unified the merged limiter status values
                                              pos2,pos1,                              // The cell-based limiter status is still holding the old value though.
                                              _tempFaceUnknowns[solverNumber],
                                              _tempStateSizedVectors[solverNumber],
                                              _tempStateSizedSquareMatrices[solverNumber]);
                #ifdef Debug
                _boundaryFaceMerges++;
                #endif
              }
            }
          endpfor
          grainSize.parallelSectionHasTerminated();

          fineGridVertex.setMergePerformed(pos1,pos2,true);
        }
      enddforx
    enddforx
}


//
// Below all methods are nop.
//
//=====================================


void exahype::mappings::SolutionRecomputation::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::SolutionRecomputation::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::SolutionRecomputation::prepareSendToWorker(
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

void exahype::mappings::SolutionRecomputation::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::mergeWithMaster(
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

void exahype::mappings::SolutionRecomputation::receiveDataFromMaster(
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

void exahype::mappings::SolutionRecomputation::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::SolutionRecomputation::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
