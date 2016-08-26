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

  MetadataHeap::getInstance().startToSendSynchronousData();
  DataHeap::getInstance().startToSendSynchronousData();
}

void exahype::mappings::MarkingForRefinement::endIteration(
    exahype::State& solverState) {
  MetadataHeap::getInstance().finishedToSendSynchronousData();
  DataHeap::getInstance().finishedToSendSynchronousData();
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
          ADERDGCellDescriptionHeap::getInstance().isValidIndex(coarseGridCell.getCellDescriptionsIndex())) {
        for (auto& pCoarse : ADERDGCellDescriptionHeap::getInstance().getData(coarseGridCell.getCellDescriptionsIndex())) {
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
        if (fineGridVertices[ fineGridVerticesEnumerator(v) ].getRefinementControl()==exahype::Vertex::Records::RefinementControl::Unrefined) {
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
//  return;

  if (localCell.isInside()) {
    if (localCell.isInitialised()) {
      // 1. Send out the solution values of this cell.
      //    Be careful with the receive order! The order must be reversed.
      // ADER-DG

      for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(localCell.getCellDescriptionsIndex())) {
        double* solution    = 0;
        int unknownsPerCell = 0;

        // TODO(Dominic): Add to docu: Make sure that ancestors and descendants
        // at fork boundaries always hold face data.
        switch(p.getType()) {
          case exahype::records::ADERDGCellDescription::Descendant:
          case exahype::records::ADERDGCellDescription::EmptyDescendant:
            p.setType(exahype::records::ADERDGCellDescription::Descendant);
            localCell.ensureNecessaryMemoryIsAllocated(p.getSolverNumber());
            break;
          case exahype::records::ADERDGCellDescription::Ancestor:
          case exahype::records::ADERDGCellDescription::EmptyAncestor:
            p.setType(exahype::records::ADERDGCellDescription::Ancestor);
            localCell.ensureNecessaryMemoryIsAllocated(p.getSolverNumber());
            break;
          case exahype::records::ADERDGCellDescription::Cell:
             solution        = DataHeap::getInstance().getData(p.getSolution()).data();
             unknownsPerCell = static_cast<exahype::solvers::ADERDGSolver*>(
                    exahype::solvers::RegisteredSolvers[p.getSolverNumber()])->getUnknownsPerCell();

             logDebug("prepareCopyToRemoteNode(...)","[solution] of solver " << p.getSolverNumber() << " sent to rank "<<toRank<<
                 ", cell: "<< cellCentre << ", level: " << level);

            DataHeap::getInstance().sendData(
                solution, unknownsPerCell, toRank, cellCentre, level,
                peano::heap::MessageType::ForkOrJoinCommunication);
            break;
          default:
            break;
          }
      }

      // Finite Volumes
      // TODO(Dominic):
//      for (std::vector<exahype::records::FiniteVolumesCellDescription>::const_iterator p =
//          FiniteVolumesCellDescriptionHeap::getInstance().getData(localCell.getCellDescriptionsIndex()).begin();
//          p!=FiniteVolumesCellDescriptionHeap::getInstance().getData(localCell.getCellDescriptionsIndex()).end(); ++p) {
//        // TODO(Dominic): Add to docu: Make sure that ancestors and descendants
//        // at fork boundaries always hold face data.
//
//        if (p->getType()==exahype::records::FiniteVolumesCellDescription::Cell) {
//          const double* solution    = DataHeap::getInstance().getData(p->getSolution()).data();
//          const int unknownsPerCell = static_cast<exahype::solvers::FiniteVolumesSolver*>(
//              exahype::solvers::RegisteredSolvers[p->getSolverNumber()])->getUnknownsPerCell();
//
//          DataHeap::getInstance().sendData(
//              solution, unknownsPerCell, toRank, cellCentre, level,
//              peano::heap::MessageType::ForkOrJoinCommunication);
//        }
//      }

      // 2. Finally, send out the metadata as last message.
      std::vector<peano::heap::records::IntegerHeapData> metadata =
          exahype::Cell::encodeMetadata(localCell.getCellDescriptionsIndex());

      MetadataHeap::getInstance().sendData(
          metadata, toRank, cellCentre, level,
          peano::heap::MessageType::ForkOrJoinCommunication);
    } else {
      if (localCell.getCellDescriptionsIndex()==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex) {
        assertion2(false,"Should not happen!",localCell.getCellDescriptionsIndex());
      } // Dead code elimination will remove this code section if Asserts/Debug flag not set.

      logDebug("prepareCopyToRemoteNode(...)","[empty] sent to rank "<<toRank<<", cell: "<<localCell.toString());

      std::vector<peano::heap::records::IntegerHeapData> metadata =
          exahype::Cell::createEncodedMetadataSequenceForInvalidCellDescriptionsIndex();

      MetadataHeap::getInstance().sendData(
          metadata, toRank, cellCentre, level,
          peano::heap::MessageType::ForkOrJoinCommunication);
    }
  }
}
// TODO(Dominic): How to deal with cell descriptions index that
// is copied from the remote rank but is a valid index on the local
// remote rank?

void exahype::mappings::MarkingForRefinement::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
//  return;

  if (localCell.isInside()) {
    // 1. First receive the metadata
    int receivedMetadataIndex = MetadataHeap::getInstance().createData(0,0);
    assertion(MetadataHeap::getInstance().getData(receivedMetadataIndex).empty());

    MetadataHeap::getInstance().receiveData(
        receivedMetadataIndex,
        fromRank, cellCentre, level,
        peano::heap::MessageType::ForkOrJoinCommunication);

    // TODO(Dominic): Check if the cellDescriptionsIndex does really belong to this cell
    // and is not a copy from the cell on the remote rank.
    // prepareCopyToRemoteNode(...) works not with copies according to the
    // documentation. Changes to the cellDescriptionsIndex would thus
    // result there in a modification of the master data which is not wanted.
    // This is why I have to do use the routine geometryInfoDoesNotMatch.
    if (!geometryInfoDoesMatch(localCell.getCellDescriptionsIndex(), cellCentre, cellSize, level)) {
      localCell.setCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
    }

    if ( !localCell.isInitialised() ) {
      logDebug("mergeWithRemoteDataDueToForkOrJoin(...)","[data] recv from rank "<<fromRank<<
          ", cell: "<< cellCentre << ", level: " << level);

      // 2. Receive finite volumes metadata. Order not interchangable with
      // receiving ADER-DG metadata.
//      receiveFiniteVolumesMetadataInMergeWithRemoteDataDueToForkOrJoin(
//          localCell,fromRank,cellCentre,cellSize,level,receivedMetadataIndex);

      // 3. Receive ADER-DG metadata. Order not interchangable with
      // receiving finite volumes metadata.
      addNewADERDGCellDescriptionsDueToForkOrJoin(
           localCell,fromRank,cellCentre,cellSize,level,receivedMetadataIndex);
    } else { // TODO(Dominic): Is this a problem?
      assertionMsg(false,"This should not happen at the moment since we do not consider joins.")
    }

    receiveADERDGCellDescriptionsDataDueToForkOrJoin(
        localCell,fromRank,cellCentre,cellSize,level,receivedMetadataIndex);

    // Clean up
    MetadataHeap::getInstance().deleteData(receivedMetadataIndex);
  }
}

bool exahype::mappings::MarkingForRefinement::geometryInfoDoesMatch(
    const int cellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>& cellSize,
    const int level) const {
  if (!ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndex)) {
    assertion1(!FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndex),
        cellDescriptionsIndex);
    return false;
  }
  if (ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex).empty() &&
      FiniteVolumesCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex).empty()) {
    return false;
  }
  // TODO(Dominic): Optimisation for multi-solver runs: Only check the first element of each.
  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)) {
    if (!tarch::la::equals(cellCentre,p.getOffset()+0.5*p.getSize()) ||
        p.getLevel()!=level) {
      return false;
    }
  }
  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)) {
    if (!tarch::la::equals(cellCentre,p.getOffset()+0.5*p.getSize()) ||
        p.getLevel()!=level) {
      return false;
    }
  }

  return true;
}

// TODO(Dominic): Untested
void exahype::mappings::MarkingForRefinement::addNewADERDGCellDescriptionsDueToForkOrJoin(
    exahype::Cell& localCell,
    const int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize,
    const int level,
    const int receivedMetadataIndex) const {
  std::vector<peano::heap::records::IntegerHeapData>::const_iterator metadataIterator=
      MetadataHeap::getInstance().getData(receivedMetadataIndex).begin();

  // 2. Read in the number of active solvers on the master rank.
  const int nADERDG = metadataIterator->getU(); ++metadataIterator;

  logDebug("mergeWithRemoteDataDueToForkOrJoin(...)","[metadata] receive from rank "<<fromRank<<
      ", cell: "<< cellCentre << ", level: " << level);


  // 3. Receive the finite volumes solution values since
  // they have been sent out first.
  // TODO(Dominic): Add to docu.
  if (nADERDG > 0) {
    const int endOfADERDGSection = 1+2*nADERDG;
    while (metadataIterator!=MetadataHeap::getInstance().getData(receivedMetadataIndex).begin()+endOfADERDGSection) {
      const int solverNumber = metadataIterator->getU(); ++metadataIterator;
      const int typeAsInt    = metadataIterator->getU(); ++metadataIterator;

      exahype::records::ADERDGCellDescription::Type type =
                    static_cast<exahype::records::ADERDGCellDescription::Type>(typeAsInt);
      assertion1(type==exahype::records::ADERDGCellDescription::Cell
              ||type==exahype::records::ADERDGCellDescription::Ancestor
              ||type==exahype::records::ADERDGCellDescription::Descendant,type);
      // We pass the lower left corner of the cell as offset.
      tarch::la::Vector<DIMENSIONS,double> offset = cellCentre - 0.5 * cellSize;

      localCell.addNewCellDescription(
          solverNumber,
          type,
          exahype::records::ADERDGCellDescription::None,
          level,
          multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex,
          cellSize,
          offset);
      localCell.ensureNecessaryMemoryIsAllocated(solverNumber);
    }
  }
}

void exahype::mappings::MarkingForRefinement::receiveADERDGCellDescriptionsDataDueToForkOrJoin(
    exahype::Cell& localCell,
    const int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize,
    const int level,
    const int receivedMetadataIndex) const {
  std::vector<peano::heap::records::IntegerHeapData>::const_iterator metadataIterator=
      MetadataHeap::getInstance().getData(receivedMetadataIndex).begin();

  // 2. Read in the number of active solvers on the master rank.
  const int nADERDG = metadataIterator->getU(); ++metadataIterator;

  // 3. Receive the finite volumes solution values since
  // they have been sent out first.
  // TODO(Dominic): Add to docu.
  if (nADERDG > 0) {
    const int endOfADERDGSection = 1+2*nADERDG;
    while (metadataIterator!=MetadataHeap::getInstance().getData(receivedMetadataIndex).begin()+endOfADERDGSection) {
      const int solverNumber = metadataIterator->getU(); ++metadataIterator;
      const int typeAsInt    = metadataIterator->getU(); ++metadataIterator;
      exahype::records::ADERDGCellDescription::Type type =
                    static_cast<exahype::records::ADERDGCellDescription::Type>(typeAsInt);
      assertion1(type==exahype::records::ADERDGCellDescription::Cell
              ||type==exahype::records::ADERDGCellDescription::Ancestor
              ||type==exahype::records::ADERDGCellDescription::Descendant,type);
      // 2. Receive solution values if necessary
      if (type==exahype::records::ADERDGCellDescription::Cell) {
        for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(localCell.getCellDescriptionsIndex())) {
          if (p.getSolverNumber() == solverNumber) {
            assertion5(tarch::la::equals(cellCentre,p.getOffset()+0.5*p.getSize()),
                cellCentre,cellSize,p.getOffset()+0.5*p.getSize(),level,p.getLevel());
            assertion2(p.getLevel()==level,p.getLevel(),level);

            // TODO(Dominic): Receive time stamps and time step sizes.
            if (p.getType()==exahype::records::ADERDGCellDescription::Cell) {
              logDebug("mergeWithRemoteDataDueToForkOrJoin(...)","[solution] recv from rank "<<fromRank<<
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
}

// TODO(Dominic): Untested
void exahype::mappings::MarkingForRefinement::receiveFiniteVolumesMetadataInMergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell,
    const int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize,
    const int level,
    const int receivedMetadataIndex) const {
  std::vector<peano::heap::records::IntegerHeapData>::const_iterator metadataIterator=
      MetadataHeap::getInstance().getData(receivedMetadataIndex).begin();

  // 2. Read in the number of active solvers on the master rank.
  const int nADERDG = metadataIterator->getU();
  metadataIterator += 1+2*nADERDG;
  const int nFV = metadataIterator->getU();

  // 3. Receive the finite volumes solution values since
  // they have been sent out first.
  // TODO(Dominic): Add to docu.
  if (nFV > 0) {
    const int endOfFiniteVolumesSection = 2+2*nADERDG+2*nFV;

    while (metadataIterator!=MetadataHeap::getInstance().getData(receivedMetadataIndex).begin()+endOfFiniteVolumesSection) {
      const int solverNumber = metadataIterator->getU(); ++metadataIterator;
      const int typeAsInt    = metadataIterator->getU(); ++metadataIterator;

      exahype::records::FiniteVolumesCellDescription::Type type =
                    static_cast<exahype::records::FiniteVolumesCellDescription::Type>(typeAsInt);
      assertion1(type==exahype::records::FiniteVolumesCellDescription::Cell
//              ||type==exahype::records::FiniteVolumesCellDescription::Ancestor
//              ||type==exahype::records::FiniteVolumesCellDescription::Descendant
                ,type);
      // We pass the lower left corner of the cell as offset.
      tarch::la::Vector<DIMENSIONS,double> offset = cellCentre - 0.5 * cellSize;

      // 1. Initiliase the cells.
      localCell.addNewCellDescription(
          solverNumber,
          type,
//              exahype::records::ADERDGCellDescription::None,
          level,
          multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex, // TODO(Dominic): Consider that in parent lookup.
          cellSize,
          offset);

      // 2. Receive solution values if necessary
      for (auto& p : FiniteVolumesCellDescriptionHeap::getInstance().getData(localCell.getCellDescriptionsIndex())) {
        if (p.getSolverNumber() == solverNumber) {
          localCell.ensureNecessaryMemoryIsAllocated(solverNumber);

          // TODO(Dominic): Receive time stamps and time step sizes.

          if (type==exahype::records::FiniteVolumesCellDescription::Cell) {
            DataHeap::getInstance().receiveData(p.getSolution(), fromRank, cellCentre, level,
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
    const MarkingForRefinement& masterThread) {
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
