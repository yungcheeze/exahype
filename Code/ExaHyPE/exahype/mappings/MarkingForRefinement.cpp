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
 
#include "exahype/mappings/MarkingForRefinement.h"

#include "peano/utils/Globals.h"

#include "peano/utils/Loop.h"

#include "tarch/la/VectorScalarOperations.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"
#include "kernels/KernelCalls.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

peano::CommunicationSpecification
exahype::mappings::MarkingForRefinement::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::MarkingForRefinement::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::MarkingForRefinement::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::MarkingForRefinement::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial);
}
peano::MappingSpecification
exahype::mappings::MarkingForRefinement::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::MarkingForRefinement::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::MarkingForRefinement::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::MarkingForRefinement::_log(
    "exahype::mappings::MarkingForRefinement");


void exahype::mappings::MarkingForRefinement::beginIteration(
    exahype::State& solverState) {
  _state = solverState;

  #ifdef Parallel
  MetadataHeap::getInstance().startToSendSynchronousData();
  DataHeap::getInstance().startToSendSynchronousData();
  #endif
}

void exahype::mappings::MarkingForRefinement::endIteration(
    exahype::State& solverState) {

  #ifdef Parallel
  MetadataHeap::getInstance().finishedToSendSynchronousData();
  DataHeap::getInstance().finishedToSendSynchronousData();
  #endif
}

void exahype::mappings::MarkingForRefinement::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);


  if ( fineGridCell.isInitialised() ) {
    logDebug( "enterCell(...)", "enter with " <<
        fineGridCell.getNumberOfADERDGCellDescriptions() << " ADER-DG solver(s) and " << fineGridCell.getNumberOfFiniteVolumeCellDescriptions() << " FV solver(s)")

    bool refineFineGridCell = false;
    for (int i=0; i<fineGridCell.getNumberOfADERDGCellDescriptions(); i++) {
      auto* uncastedSolver = exahype::solvers::RegisteredSolvers[ fineGridCell.getADERDGCellDescription(i).getSolverNumber() ];
      assertion( uncastedSolver->getType()==exahype::solvers::Solver::Type::ADER_DG );
      exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(uncastedSolver);

      // TODO(Dominic): Initialise new fine grid cell and coarse grid cell relationships
      // for cells initialised in forking events.
      // if (fineGridCell.parentCell==Remote && coarseGridCell.parentCellIndex==Remote)
      #ifdef Parallel
      if (fineGridCell.getADERDGCellDescription(i).getParentIndex()==multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex &&
          exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(coarseGridCell.getCellDescriptionsIndex())) {
        for (auto& pCoarse : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(coarseGridCell.getCellDescriptionsIndex())) {
          if (pCoarse.getSolverNumber()==fineGridCell.getADERDGCellDescription(i).getSolverNumber()) {
            fineGridCell.getADERDGCellDescription(i).setParentIndex(coarseGridCell.getCellDescriptionsIndex());
          }
        }
      }
      #endif

      if (tarch::la::allSmallerEquals(fineGridVerticesEnumerator.getCellSize(),solver->getMaximumMeshSize()) &&
          fineGridCell.getADERDGCellDescription(i).getType()==exahype::records::ADERDGCellDescription::Cell) {
        double*  solution = 0;
        exahype::solvers::Solver::RefinementControl refinementControl;

        switch (fineGridCell.getADERDGCellDescription(i).getRefinementEvent()) {
          case exahype::records::ADERDGCellDescription::RefiningRequested:
            refineFineGridCell = true;
            break;
          case exahype::records::ADERDGCellDescription::None:
             solution = DataHeap::getInstance().getData(fineGridCell.getADERDGCellDescription(i).getSolution()).data();
             refinementControl = solver->refinementCriterion(
                solution, fineGridVerticesEnumerator.getCellCenter(),
                fineGridVerticesEnumerator.getCellSize(),
                fineGridCell.getADERDGCellDescription(i).getCorrectorTimeStamp(), fineGridCell.getADERDGCellDescription(i).getLevel());

            switch (refinementControl) {
              case exahype::solvers::Solver::RefinementControl::Refine:
                fineGridCell.getADERDGCellDescription(i).setRefinementEvent(exahype::records::ADERDGCellDescription::RefiningRequested);
                refineFineGridCell = true;
              break;
              case exahype::solvers::Solver::RefinementControl::Erase:
                fineGridCell.getADERDGCellDescription(i).setRefinementEvent(exahype::records::ADERDGCellDescription::ErasingRequested);
              break;
              default:
                break;
            }
            break;
          default:
            break;
          }
      }
    }

//    logDebug( "enterCell(...)", "continue  with " << fineGridCell.getNumberOfFiniteVolumeCellDescriptions() << " FV solver(s)")
//    for (int i=0; i<fineGridCell.getNumberOfFiniteVolumeCellDescriptions(); i++) {
//      auto* uncastedSolver = exahype::solvers::RegisteredSolvers[ fineGridCell.getFiniteVolumesCellDescription(i).getSolverNumber() ];
//      assertion( uncastedSolver->getType()==exahype::solvers::Solver::Type::FiniteVolumes );
//      exahype::solvers::FiniteVolumesSolver* solver = static_cast<exahype::solvers::FiniteVolumesSolver*>(uncastedSolver);
//
//      if (
//        tarch::la::allSmallerEquals(fineGridVerticesEnumerator.getCellSize(),solver->getMaximumMeshSize())
//        &&
//        fineGridCell.getFiniteVolumesCellDescription(i).getType()==exahype::records::FiniteVolumesCellDescription::Cell
///*
//        &&
//        fineGridCell.FiniteVolumesSolver(i).getRefinementEvent()==exahype::records::ADERDGCellDescription::None
//*/
//      ) {
//        double*  solution = DataHeap::getInstance().getData(fineGridCell.getFiniteVolumesCellDescription(i).getSolution()).data();
//
//        exahype::solvers::Solver::RefinementControl refinementControl = solver->refinementCriterion(
//            solution, fineGridVerticesEnumerator.getCellCenter(),
//            fineGridVerticesEnumerator.getCellSize(),
//            fineGridCell.getFiniteVolumesCellDescription(i).getTimeStamp(),
//            fineGridCell.getFiniteVolumesCellDescription(i).getTimeStepSize());
//
//        switch (refinementControl) {
//          case exahype::solvers::Solver::RefinementControl::Refine:
////            fineGridCell.getADERDGCellDescription(i).setRefinementEvent(exahype::records::ADERDGCellDescription::RefiningRequested);
//
//            dfor2(v)
//              if (fineGridVertices[ fineGridVerticesEnumerator(v) ].getRefinementControl()==exahype::Vertex::Records::RefinementControl::Unrefined) {
//                fineGridVertices[ fineGridVerticesEnumerator(v) ].refine();
//              }
//            enddforx
//          break;
//          case exahype::solvers::Solver::RefinementControl::Erase:
//            assertionMsg( false, "not implemented yet" );
//            //fineGridCell.getADERDGCellDescription(i).setRefinementEvent(exahype::records::ADERDGCellDescription::ErasingRequested);
//          break;
//          default:
//            break;
//        }
//      }
//    }

//    if (refineFineGridCell && _state.refineInitialGridInTouchVertexLastTime()) {
    if (refineFineGridCell) {
      dfor2(v)
        if (
          fineGridVertices[ fineGridVerticesEnumerator(v) ].getRefinementControl()==exahype::Vertex::Records::RefinementControl::Unrefined
	  &&
	  _state.refineInitialGridInTouchVertexLastTime()
	) {
          fineGridVertices[ fineGridVerticesEnumerator(v) ].refine();
        }
      enddforx
    }

    logDebug( "enterCell(...)", "left with " << fineGridCell.getNumberOfADERDGCellDescriptions() << " ADER-DG solver(s) and " << fineGridCell.getNumberOfFiniteVolumeCellDescriptions() << " FV solver(s)")
  }

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

#ifdef Parallel
void exahype::mappings::MarkingForRefinement::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  return;

  if (localCell.isInside() && localCell.isInitialised()) {
    exahype::solvers::ADERDGSolver::sendCellDescriptions(toRank,localCell.getCellDescriptionsIndex(),
        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);
//    exahype::solvers::FiniteVolumesSolver::sendCellDescriptions(toRank,localCell.getCellDescriptionsIndex(),
//        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);

    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      int element = solver->tryGetElement(localCell.getCellDescriptionsIndex(),solverNumber);

      if(element!=exahype::solvers::Solver::NotFound) {
        solver->sendDataToWorkerOrMasterDueToForkOrJoin(
            toRank,localCell.getCellDescriptionsIndex(),element,cellCentre,level);
      } else {
        solver->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(toRank,cellCentre,level);
      }

      ++solverNumber;
    }
  } else if (localCell.isInside() && !localCell.isInitialised()){
    exahype::solvers::ADERDGSolver::sendEmptyCellDescriptions(toRank,
        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);
    exahype::solvers::FiniteVolumesSolver::sendEmptyCellDescriptions(toRank,
        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);

    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      solver->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(toRank,cellCentre,level);
      ++solverNumber;
    }
  }
}

//// TODO(Dominic): Add to docu: Make sure that ancestors and descendants
//// at fork boundaries always hold face data.
//// TODO(Dominic): Move change of type out of here to MarkingForAugmentation::enterCell(...)
//void exahype::mappings::MarkingForRefinement::sendADERDGDataToMasterOrWorker(
//    int cellDescriptionsIndex,
//    int toRank,
//    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
//    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
//  for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex)) {
//    double* solution    = 0;
//    int unknownsPerCell = 0;
//
//    switch(p.getType()) {
//      case exahype::records::ADERDGCellDescription::Descendant:
//      case exahype::records::ADERDGCellDescription::EmptyDescendant:
//        p.setType(exahype::records::ADERDGCellDescription::Descendant);
//        exahype::Cell::ensureNecessaryMemoryIsAllocated(p);
//        break;
//      case exahype::records::ADERDGCellDescription::Ancestor:
//      case exahype::records::ADERDGCellDescription::EmptyAncestor:
//        p.setType(exahype::records::ADERDGCellDescription::Ancestor);
//        exahype::Cell::ensureNecessaryMemoryIsAllocated(p);
//        break;
//      case exahype::records::ADERDGCellDescription::Cell:
//         solution        = DataHeap::getInstance().getData(p.getSolution()).data();
//         unknownsPerCell = static_cast<exahype::solvers::ADERDGSolver*>(
//                exahype::solvers::RegisteredSolvers[p.getSolverNumber()])->getUnknownsPerCell();
//
//         logDebug("sendADERDGDataToMasterOrWorker(...)","[solution] of solver " << p.getSolverNumber() << " sent to rank "<<toRank<<
//             ", cell: "<< cellCentre << ", level: " << level);
//
//        DataHeap::getInstance().sendData(
//            solution, unknownsPerCell, toRank, cellCentre, level,
//            peano::heap::MessageType::ForkOrJoinCommunication);
//        break;
//      default:
//        break;
//      }
//  }
//}

void exahype::mappings::MarkingForRefinement::receiveDataFromMaster(
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

// TODO(Dominic): How to deal with cell descriptions index that
// is copied from the remote rank but is a valid index on the local
// remote rank? Currently use geometryInfoDoesMatch! Not best idea.
void exahype::mappings::MarkingForRefinement::mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  return; // TODO(Dominic): Remove

  if (localCell.isInside()) {
    if (!geometryInfoDoesMatch(localCell.getCellDescriptionsIndex(),cellCentre,cellSize,level)) {
      localCell.setCellDescriptionsIndex(
          multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
    }
    if ( !localCell.isInitialised() ) {
      localCell.setupMetaData();
    }

    exahype::solvers::ADERDGSolver::mergeCellDescriptionsWithRemoteData(
        fromRank,localCell.getCellDescriptionsIndex(),
        peano::heap::MessageType::ForkOrJoinCommunication,
        cellCentre,level);
//    exahype::solvers::FiniteVolumesSolver::mergeCellDescriptionsWithRemoteData(
//        fromRank,localCell.getCellDescriptionsIndex(),
//        peano::heap::MessageType::ForkOrJoinCommunication,
//        cellCentre,level);

    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      int element = solver->tryGetElement(localCell.getCellDescriptionsIndex(),solverNumber);

      if (element!=exahype::solvers::Solver::NotFound) {
        solver->mergeWithWorkerOrMasterDataDueToForkOrJoin(
            fromRank,localCell.getCellDescriptionsIndex(),element,cellCentre,level);
      } else {
        solver->dropWorkerOrMasterDataDueToForkOrJoin(fromRank,cellCentre,level);
      }
      ++solverNumber;
    }
  }
}

bool exahype::mappings::MarkingForRefinement::geometryInfoDoesMatch(
    const int cellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>& cellSize,
    const int level) {
  if (!exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex)) {
    assertion1(!exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex),
        cellDescriptionsIndex);
    return false;
  }
  if (exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex).empty() &&
      exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex).empty()) {
    return false;
  }
  // TODO(Dominic): Optimisation for multi-solver runs: Only check the first element of each.
  for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex)) {
    if (!tarch::la::equals(cellCentre,p.getOffset()+0.5*p.getSize()) ||
        p.getLevel()!=level) {
      return false;
    }
  }
  for (auto& p : exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex)) {
    if (!tarch::la::equals(cellCentre,p.getOffset()+0.5*p.getSize()) ||
        p.getLevel()!=level) {
      return false;
    }
  }

  return true;
}

void exahype::mappings::MarkingForRefinement::receiveADERDGDataFromMasterOrWorker(
    const int cellDescriptionsIndex,
    const int fromRank,
    const peano::heap::MessageType& messageType,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize,
    const int level,
    const int receivedMetadataIndex) {
  assertion1(messageType==peano::heap::MessageType::MasterWorkerCommunication ||
            messageType==peano::heap::MessageType::ForkOrJoinCommunication,peano::heap::toString(messageType));
  assertion(!MetadataHeap::getInstance().getData(receivedMetadataIndex).empty());
  exahype::MetadataHeap::HeapEntries neighbourCellTypes =
      MetadataHeap::getInstance().getData(receivedMetadataIndex);

  for (int solverNumber=exahype::solvers::RegisteredSolvers.size()-1; solverNumber>0; --solverNumber) {
    const int typeAsInt = neighbourCellTypes[solverNumber].getU();
    exahype::records::ADERDGCellDescription::Type neighbourType =
        static_cast<exahype::records::ADERDGCellDescription::Type>(typeAsInt);
    assertion1(neighbourType==exahype::records::ADERDGCellDescription::Cell
        ||neighbourType==exahype::records::ADERDGCellDescription::Ancestor
        ||neighbourType==exahype::records::ADERDGCellDescription::Descendant,neighbourType);
    // 2. Receive solution values if necessary
    if (neighbourType==exahype::records::ADERDGCellDescription::Cell) {
      for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex)) {
        if (p.getSolverNumber() == solverNumber) {
          assertion5(tarch::la::equals(cellCentre,p.getOffset()+0.5*p.getSize()),
              cellCentre,cellSize,p.getOffset()+0.5*p.getSize(),level,p.getLevel());
          assertion2(p.getLevel()==level,p.getLevel(),level);

          if (p.getType()==exahype::records::ADERDGCellDescription::Cell) {
            logDebug("mergeWithRemoteDataDueToForkOrJoin(...)","[solution] receive from rank "<<fromRank<<
                ", cell: "<< cellCentre << ", level: " << level);

            DataHeap::getInstance().receiveData(
                p.getSolution(), fromRank, cellCentre, level,
                peano::heap::MessageType::ForkOrJoinCommunication);
          }
        }
      }
    }
  }
}

//
// Below all methods are nop.
//
// ====================================



void exahype::mappings::MarkingForRefinement::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

bool exahype::mappings::MarkingForRefinement::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  // do nothing
  return true;
}

void exahype::mappings::MarkingForRefinement::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::mergeWithMaster(
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

void exahype::mappings::MarkingForRefinement::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::MarkingForRefinement::MarkingForRefinement() {
  // do nothing
}

exahype::mappings::MarkingForRefinement::~MarkingForRefinement() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::MarkingForRefinement::MarkingForRefinement(
    const MarkingForRefinement& masterThread):
    _state(masterThread._state) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::mergeWithWorkerThread(
    const MarkingForRefinement& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::MarkingForRefinement::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::MarkingForRefinement::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
