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
 
#include "exahype/mappings/RiemannSolver.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/utils/Globals.h"

#include "tarch/la/VectorScalarOperations.h"

#include "tarch/multicore/Loop.h"

#include "exahype/solvers/ADERDGSolver.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"


peano::CommunicationSpecification
exahype::mappings::RiemannSolver::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}


peano::MappingSpecification
exahype::mappings::RiemannSolver::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification
exahype::mappings::RiemannSolver::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification
exahype::mappings::RiemannSolver::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification
exahype::mappings::RiemannSolver::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification
exahype::mappings::RiemannSolver::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}


peano::MappingSpecification
exahype::mappings::RiemannSolver::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}



tarch::logging::Log exahype::mappings::RiemannSolver::_log( "exahype::mappings::RiemannSolver");


void exahype::mappings::RiemannSolver::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("touchVertexFirstTime(...)", fineGridVertex,
                           fineGridX, fineGridH,
                           coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  tarch::la::Vector<TWO_POWER_D, int>& adjacentCellDescriptionsIndices =
      fineGridVertex.getCellDescriptionsIndex();
  logDebug(
      "touchVertexFirstTime(...)",
      "cell descriptions around vertex. "
          << "coarse grid level: " << coarseGridVerticesEnumerator.getLevel()
          << ", fine grid position:" << fineGridPositionOfVertex
          << ", adjacent cell descriptions indices:"
          << adjacentCellDescriptionsIndices);
  logDebug("touchVertexFirstTime(...)", "cell descriptions around vertex. "
                                            << "fine grid x " << fineGridX);

  /* Right cell-left cell   pair indices: 0,1; 2,3;   4,5; 6;7
   * Front cell-back cell   pair indices: 0,2; 1,3;   4,6; 5;7
   * Top   cell-bottom cell pair indices: 0,4; 1,5;   2,6; 3;7
   *
   * Note that from the viewpoint of a cell, the face
   * has always the "opposite" index, i.e., we solve a Riemann
   * problem on the left face of the right cell (which
   * is the right face of the left cell).
   */
  constexpr int cellIndicesLeft[4] = {0, 2, 4, 6};
  constexpr int cellIndicesRight[4] = {1, 3, 5, 7};
  constexpr int cellIndicesFront[4] = {0, 1, 4, 5};
  constexpr int cellIndicesBack[4] = {2, 3, 6, 7};
#if DIMENSIONS == 3
  constexpr int cellIndicesBottom[4] = {0, 1, 2, 3};
  constexpr int cellIndicesTop[4] = {4, 5, 6, 7};
#endif
  for (int i = 0; i < TWO_POWER_D_DIVIDED_BY_TWO; i++) {
    logDebug("touchVertexFirstTime()","fineGridVertex.isBoundary(): " << fineGridVertex.isBoundary());

    solveRiemannProblemAtInterface(
      adjacentCellDescriptionsIndices[cellIndicesLeft[i]],
      adjacentCellDescriptionsIndices[cellIndicesRight[i]],
      EXAHYPE_FACE_RIGHT, EXAHYPE_FACE_LEFT,
      0);

    Cell::mergeSolutionMinMaxOnFace(
      adjacentCellDescriptionsIndices[cellIndicesLeft[i]],
      adjacentCellDescriptionsIndices[cellIndicesRight[i]],
      EXAHYPE_FACE_RIGHT, EXAHYPE_FACE_LEFT
    );

    solveRiemannProblemAtInterface(
      adjacentCellDescriptionsIndices[cellIndicesFront[i]],
      adjacentCellDescriptionsIndices[cellIndicesBack[i]],
      EXAHYPE_FACE_BACK, EXAHYPE_FACE_FRONT,
      1);

    Cell::mergeSolutionMinMaxOnFace(
      adjacentCellDescriptionsIndices[cellIndicesFront[i]],
      adjacentCellDescriptionsIndices[cellIndicesBack[i]],
      EXAHYPE_FACE_BACK, EXAHYPE_FACE_FRONT
    );

    #if DIMENSIONS == 3
    solveRiemannProblemAtInterface(
      adjacentCellDescriptionsIndices[cellIndicesBottom[i]],
      adjacentCellDescriptionsIndices[cellIndicesTop[i]],
      EXAHYPE_FACE_TOP, EXAHYPE_FACE_BOTTOM,
      2);

    Cell::mergeSolutionMinMaxOnFace(
      adjacentCellDescriptionsIndices[cellIndicesBottom[i]],
      adjacentCellDescriptionsIndices[cellIndicesTop[i]],
      EXAHYPE_FACE_TOP, EXAHYPE_FACE_BOTTOM
    );
    #endif
  }

  logTraceOutWith1Argument("touchVertexFirstTime(...)", fineGridVertex);
}

void exahype::mappings::RiemannSolver::solveRiemannProblemAtInterface(
    const int cellDescriptionsIndexOfLeftCell,
    const int cellDescriptionsIndexOfRightCell,
    const int faceIndexForLeftCell,
    const int faceIndexForRightCell,
    const int normalNonZero) {
  // Only continue if this is an internal face, i.e.,
  // both cell description indices are valid
  if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndexOfLeftCell) &&
      ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndexOfRightCell)) {
    logDebug("touchVertexLastTime(...)::solveRiemannProblemAtInterface(...)",
             "Performing Riemann solve. "
                 << "faceIndexForLeftCell:" << faceIndexForLeftCell
                 << " faceIndexForRightCell:" << faceIndexForRightCell);
                 // << " indexOfLeftCell:"
                 // << adjacentADERDGCellDescriptionsIndices[indexOfLeftCell]
                 // << " indexOfRightCell:"
                 // << adjacentADERDGCellDescriptionsIndices[indexOfRightCell]);
    std::vector<records::ADERDGCellDescription>& cellDescriptionsOfLeftCell =
        ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndexOfLeftCell);
    int numberOfADERDGCellDescriptionsLeft =
        ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndexOfLeftCell).size();
    const peano::datatraversal::autotuning::MethodTrace methodTrace = peano::datatraversal::autotuning::UserDefined4;
    const int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfADERDGCellDescriptionsLeft, methodTrace);

    pfor(i, 0, numberOfADERDGCellDescriptionsLeft,grainSize)
      for (auto& cellDescriptionOfRightCell :
          ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndexOfRightCell)) {
        if (cellDescriptionsOfLeftCell [i].getType() == exahype::records::ADERDGCellDescription::Cell ||
            cellDescriptionOfRightCell.getType() == exahype::records::ADERDGCellDescription::Cell) {
        assertion1(cellDescriptionsOfLeftCell[i].getType() == exahype::records::ADERDGCellDescription::Cell ||
                   cellDescriptionsOfLeftCell[i].getType() == exahype::records::ADERDGCellDescription::Ancestor ||
                   cellDescriptionsOfLeftCell[i].getType() == exahype::records::ADERDGCellDescription::Descendant,
                   cellDescriptionsOfLeftCell[i].toString());
        assertion1(cellDescriptionOfRightCell.getType() == exahype::records::ADERDGCellDescription::Cell ||
                   cellDescriptionOfRightCell.getType() == exahype::records::ADERDGCellDescription::Ancestor ||
                   cellDescriptionOfRightCell.getType() == exahype::records::ADERDGCellDescription::Descendant,
                   cellDescriptionOfRightCell.toString());
        assertion1(cellDescriptionsOfLeftCell[i].getRefinementEvent()==exahype::records::ADERDGCellDescription::None,
                   cellDescriptionsOfLeftCell[i].toString());
        assertion1(cellDescriptionOfRightCell.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,
                   cellDescriptionOfRightCell.toString());
        assertionEquals4(cellDescriptionsOfLeftCell[i].getRiemannSolvePerformed(faceIndexForLeftCell),
                         cellDescriptionOfRightCell.getRiemannSolvePerformed(faceIndexForRightCell),
                         faceIndexForLeftCell, faceIndexForRightCell,
                         cellDescriptionsOfLeftCell[i].toString(),
                         cellDescriptionOfRightCell.toString());
        exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*> (
            exahype::solvers::RegisteredSolvers[cellDescriptionsOfLeftCell[i].getSolverNumber()] );


        if (!cellDescriptionsOfLeftCell[i].getRiemannSolvePerformed(faceIndexForLeftCell)) {
          cellDescriptionsOfLeftCell [i].setRiemannSolvePerformed(faceIndexForLeftCell,  true);
          cellDescriptionOfRightCell.setRiemannSolvePerformed(faceIndexForRightCell, true);
          #ifdef Debug
          _interiorFaceSolves++;
          #endif

          const int numberOfFaceDof = solver->getUnknownsPerFace();

          double* QL = DataHeap::getInstance() .getData(cellDescriptionsOfLeftCell[i].getExtrapolatedPredictor()).data() +
              (faceIndexForLeftCell * numberOfFaceDof);
          double* QR = DataHeap::getInstance().getData(cellDescriptionOfRightCell.getExtrapolatedPredictor()).data() +
              (faceIndexForRightCell * numberOfFaceDof);
          double* FL = DataHeap::getInstance().getData(cellDescriptionsOfLeftCell[i].getFluctuation()).data() +
              (faceIndexForLeftCell * numberOfFaceDof);
          double* FR = DataHeap::getInstance().getData(cellDescriptionOfRightCell.getFluctuation()).data() +
              (faceIndexForRightCell * numberOfFaceDof);

          for(int ii=0; ii<numberOfFaceDof; ++ii) {
            assertion(std::isfinite(QL[ii]));
            assertion(std::isfinite(QR[ii]));
            assertion(std::isfinite(FL[ii]));
            assertion(std::isfinite(FR[ii]));
          }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

          // Synchronise time stepping.
          solver->synchroniseTimeStepping(cellDescriptionsOfLeftCell[i]);
          solver->synchroniseTimeStepping(cellDescriptionOfRightCell);

          // todo Time step must be interpolated in local time stepping case
          // both time step sizes are the same, so the min has no effect here.
          solver->riemannSolver(
              FL, FR, QL, QR,
              std::min(cellDescriptionsOfLeftCell[i].getCorrectorTimeStepSize(),
                       cellDescriptionOfRightCell.getCorrectorTimeStepSize()),
              normalNonZero);

          for(int i=0; i<numberOfFaceDof; ++i) {
            assertion(std::isfinite(FL[i]));
            assertion(std::isfinite(FR[i]));
          }  // Dead code elimination will get rid of this loop if Asserts flag is not set.
        }
      }
    }
    endpfor peano::datatraversal::autotuning::Oracle::getInstance()
        .parallelSectionHasTerminated(methodTrace);
  } else if (
      (ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndexOfLeftCell) &&
          (cellDescriptionsIndexOfRightCell == multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex ||
          cellDescriptionsIndexOfRightCell == multiscalelinkedcell::HangingVertexBookkeeper::RemoteAndDomainBoundaryAdjacencyIndex))
          ||
          (ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndexOfRightCell) &&
              (cellDescriptionsIndexOfLeftCell == multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex ||
              cellDescriptionsIndexOfLeftCell == multiscalelinkedcell::HangingVertexBookkeeper::RemoteAndDomainBoundaryAdjacencyIndex))

  ) {
    #if defined(PeriodicBC)
      assertionMsg(false,"PeriodicBC: Please implement!");
      return;
    #endif
      // TODO(Dominic): Previous code for reference:
  //  } else if ((ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndexOfLeftCell) &&
  //      cellDescriptionsIndexOfRightCell == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex) ||
  //          (ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndexOfRightCell) &&
  //              cellDescriptionsIndexOfLeftCell == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex)) {

    int cellDescriptionsIndex = cellDescriptionsIndexOfLeftCell;
    int faceIndex             = faceIndexForLeftCell;

    if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndexOfRightCell)) {
      cellDescriptionsIndex = cellDescriptionsIndexOfRightCell;
      faceIndex             = faceIndexForRightCell;

      assertion1(!ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndexOfLeftCell),
          cellDescriptionsIndexOfLeftCell);
    } else {
      assertion1(!ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                cellDescriptionsIndexOfRightCell),cellDescriptionsIndexOfRightCell);
    }  // Dead code elimination will get rid of this loop if Asserts flag is not set.
    assertion(ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndex));

    std::vector<records::ADERDGCellDescription>& cellDescriptions =
        ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex);

    const int numberOfADERDGCellDescriptions =
        static_cast<int>(cellDescriptions.size());
    const peano::datatraversal::autotuning::MethodTrace methodTrace =
        peano::datatraversal::autotuning::UserDefined0;
    const int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
            numberOfADERDGCellDescriptions, methodTrace);

    pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
    if (!cellDescriptions[i].getRiemannSolvePerformed(faceIndex) &&
        cellDescriptions[i].getType() == exahype::records::ADERDGCellDescription::Cell) {
      #ifdef Debug
      _boundaryFaceSolves++;
      #endif
      cellDescriptions[i].setRiemannSolvePerformed(faceIndex, true);
      applyBoundaryConditions(cellDescriptions[i], faceIndex, normalNonZero);
    }
    endpfor
  } // else do nothing
}

// Verified correct calling of this method for 9x9 grid on [0,1]x[0,1].
void exahype::mappings::RiemannSolver::applyBoundaryConditions(
    records::ADERDGCellDescription& cellDescription, const int faceIndexForCell,
    const int normalNonZero) {
  assertion1(cellDescription.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,cellDescription.toString());
  exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
      exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);

  const int numberOfFaceDof = solver->getUnknownsPerFace();

  double* stateIn = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
      (faceIndexForCell * numberOfFaceDof);
  double* fluxIn = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
      (faceIndexForCell * numberOfFaceDof);

  for(int ii=0; ii<numberOfFaceDof; ++ii) {
    assertion5(std::isfinite(stateIn[ii]), cellDescription.toString(),
        faceIndexForCell, normalNonZero, ii, stateIn[ii]);
    assertion5(std::isfinite(fluxIn[ii]), cellDescription.toString(),
        faceIndexForCell, normalNonZero, ii, fluxIn[ii]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  double* stateOut = new double[numberOfFaceDof];
  double* fluxOut  = new double[numberOfFaceDof];

  tarch::la::Vector<DIMENSIONS,double> cellCentre = cellDescription.getOffset() + 0.5*cellDescription.getSize();

  logDebug("applyBoundaryConditions(...)", "face index: " << faceIndexForCell <<
      ", cell lower-left corner " << cellDescription.getOffset() <<
      ", cell upper-right corner " << (cellDescription.getOffset()+cellDescription.getSize()));

  // Synchronise time stepping.
  solver->synchroniseTimeStepping(cellDescription);

  solver->boundaryConditions(fluxOut,stateOut,
                             fluxIn,stateIn,
                             cellCentre,cellDescription.getSize(),
                             cellDescription.getCorrectorTimeStamp(),
                             cellDescription.getCorrectorTimeStepSize(),
                             faceIndexForCell,normalNonZero);

  for(int ii=0; ii<numberOfFaceDof; ++ii) {
    assertion5(std::isfinite(stateOut[ii]), cellDescription.toString(), faceIndexForCell, normalNonZero, ii, stateOut[ii]);
    assertion5(std::isfinite(fluxOut[ii]), cellDescription.toString(), faceIndexForCell, normalNonZero, ii, fluxOut[ii]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  // @todo(Dominic): Add to docu why we need this.
  if (faceIndexForCell % 2 == 0) {
    solver->riemannSolver(fluxOut, fluxIn, stateOut, stateIn,
        cellDescription.getCorrectorTimeStepSize(),
        normalNonZero);
  } else {
    solver->riemannSolver(fluxIn, fluxOut, stateIn, stateOut,
        cellDescription.getCorrectorTimeStepSize(),
        normalNonZero);
  }

  for(int ii=0; ii<numberOfFaceDof; ++ii) {
    assertion5(std::isfinite(fluxIn[ii]), cellDescription.toString(),
               faceIndexForCell, normalNonZero, ii, fluxIn[ii]);
    assertion5(std::isfinite(fluxOut[ii]), cellDescription.toString(),
               faceIndexForCell, normalNonZero, ii, fluxOut[ii]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  delete[] stateOut;
  delete[] fluxOut;
}

/**
 * TODO(Dominic): Add docu. Returns zero otherwise.
 */

void exahype::mappings::RiemannSolver::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::RiemannSolver::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  #ifdef Parallel
  _state = &solverState;
  #endif

  #ifdef Debug
  _interiorFaceSolves = 0;
  _boundaryFaceSolves = 0;
  #endif

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::RiemannSolver::endIteration(
    exahype::State& solverState) {
  logDebug("endIteration(...)","interiorFaceSolves: " << _interiorFaceSolves);
  logDebug("endIteration(...)","boundaryFaceSolves: " << _boundaryFaceSolves);
}

#ifdef Parallel
// TODO(Dominic): Add to docu: We receive metadata on all vertices.
// We receive only for Cell/Ancestor/Descendants cell descriptions face data.
// EmptyAncestor/EmptyDescendants/InvalidAdjacencyIndices drop face data
// that was sent to them by Cells/Ancestors/Descendants.
// TODO(Dominic): Add to docu: The following invariant must hold:
// A cell holding Cell/Ancestor/Descendant cell descriptions
// next to a remote cell with cellDescriptionsIndex==InvalidAdjacencyIndex.
void exahype::mappings::RiemannSolver::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {

// TODO(Dominic): This is a bug or needs to be documented.
// see discussion in SpaceTimePredictor
//#if !defined(PeriodicBC)
//  if (vertex.isBoundary()) return;
//#endif

  // TODO(Dominic): Add to docu why we remove the vertex.isInside() constraint here.
  // We might also consider to remove it from the grid setup mapping functions.
  // fineGridCell.isInside does not imply that all adjacent vertices are
  // inside. If we count down the counter only on
  // fineGridVertices that are inside we might not send out all faces
  // of a cell that is close to the boundary.

  dfor2(myDest)
    dfor2(mySrc)
      tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
      tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;

      int destScalar = TWO_POWER_D - myDestScalar - 1;
      int srcScalar  = TWO_POWER_D - mySrcScalar  - 1;

      if (vertex.isInside() && // TODO(Dominic): Discuss with Tobias what to do for PeriodicBC
          vertex.hasToReceiveMetadata(_state,src,dest,fromRank)) {
        // 1. Receive metadata.
        int receivedMetadataIndex = MetadataHeap::getInstance().createData(0,0);
        assertion(MetadataHeap::getInstance().getData(receivedMetadataIndex).empty());
        MetadataHeap::getInstance().receiveData(
            receivedMetadataIndex,
            fromRank, fineGridX, level,
            peano::heap::MessageType::NeighbourCommunication);
        MetadataHeap::HeapEntries receivedMetadata = MetadataHeap::getInstance().getData(receivedMetadataIndex);

        const int destCellDescriptionIndex = vertex.getCellDescriptionsIndex()[destScalar];

        if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(destCellDescriptionIndex)) {
          logDebug(
              "mergeWithNeighbour(...)", "[data] " <<
              fromRank << " for vertex x=" << fineGridX << ", level=" << level <<
              ", src type=" << multiscalelinkedcell::indexToString(destCellDescriptionIndex) <<
              ", src=" << src << ", dest=" << dest);

          // TODO(Dominic): Add to docu: Order is important.
//            receiveFiniteVolumesFaceData(
//                fromRank,fineGridX,level,src,dest,
//                vertex.getCellDescriptionsIndex()[srcScalar],
//                vertex.getCellDescriptionsIndex()[destScalar],
//                receivedMetadataIndex);

          receiveADERDGFaceData(
              fromRank,fineGridX,level,src,dest,
              vertex.getCellDescriptionsIndex()[srcScalar],
              vertex.getCellDescriptionsIndex()[destScalar],
              receivedMetadata);


        } else {
          // TODO (Dominic): Drop the incoming messages if necessary.
          // We only send non trivial metadata out if we
          // send face data. See SpaceTimePredictor for details.
          // We can thus interpret from the received metadata
          // if the current vertex needs to receive and drop
          // facedata.
          if (!exahype::Cell::isEncodedMetadataSequenceForInvalidCellDescriptionsIndex(receivedMetadata)) {
            assertion1(destCellDescriptionIndex==multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex,
                destCellDescriptionIndex);

            dropADERDGFaceData(
                fromRank,fineGridX,level,src,dest,
                vertex.getCellDescriptionsIndex()[srcScalar],
                vertex.getCellDescriptionsIndex()[destScalar],
                receivedMetadata);
          }
        }
        // Clean up
        MetadataHeap::getInstance().deleteData(receivedMetadataIndex);
      }
    enddforx
  enddforx
}

void exahype::mappings::RiemannSolver::receiveADERDGFaceData(
    int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    int level,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    int srcCellDescriptionIndex,
    int destCellDescriptionIndex,
    exahype::MetadataHeap::HeapEntries& receivedMetadata) {
  ADERDGCellDescriptionHeap::getInstance().isValidIndex(destCellDescriptionIndex);
  FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(destCellDescriptionIndex);

  MetadataHeap::HeapEntries::const_iterator metadataIterator=receivedMetadata.begin();
  const int nADERDG = metadataIterator->getU();

  if (nADERDG > 0) {
    logDebug("receiveADERDGFaceData(...)","nADERDG: " << nADERDG);
    const int endOfADERDGMetadata = 2*nADERDG;
    metadataIterator += endOfADERDGMetadata;

    const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
    assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
    const int faceIndex = 2 * normalOfExchangedFace +
        (src(normalOfExchangedFace) > dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

    // TODO(Dominic): Add to docu why we need to invert the order.
    while (metadataIterator!=receivedMetadata.begin()) {
      const int neighbourTypeAsInt    = metadataIterator->getU(); --metadataIterator;
      const int neighbourSolverNumber = metadataIterator->getU(); --metadataIterator;
      exahype::records::ADERDGCellDescription::Type neighbourType =
          static_cast<exahype::records::ADERDGCellDescription::Type>(neighbourTypeAsInt);
      logDebug("receiveADERDGFaceData(...)","neighbourSolverNumber: " << neighbourSolverNumber);
      logDebug("receiveADERDGFaceData(...)","neighbourTypeAsInt: "    << neighbourTypeAsInt);

      if (neighbourType==exahype::records::ADERDGCellDescription::Cell ||
          neighbourType==exahype::records::ADERDGCellDescription::Ancestor ||
          neighbourType==exahype::records::ADERDGCellDescription::Descendant) {
        for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(destCellDescriptionIndex)) {
          if (p.getSolverNumber()==neighbourSolverNumber) {
            if (p.getFaceDataExchangeCounter(faceIndex)==0) {
              p.setFaceDataExchangeCounter(faceIndex,TWO_POWER_D); // TODO(Dominic): Add to docu what we do here with the counter.
              assertion1(!p.getRiemannSolvePerformed(faceIndex),p.toString());

              exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
                  exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);

              if (p.getType()==exahype::records::ADERDGCellDescription::Cell ||
                  p.getType()==exahype::records::ADERDGCellDescription::Ancestor ||
                  p.getType()==exahype::records::ADERDGCellDescription::Descendant){
                assertion(DataHeap::getInstance().isValidIndex(p.getExtrapolatedPredictor()));
                assertion(DataHeap::getInstance().isValidIndex(p.getFluctuation()));

                if (srcCellDescriptionIndex == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex) {
                  #if defined(PeriodicBC)
                  assertionMsg(false, "Vasco, we have to implement this");
                  #else
                  assertionMsg(false, "should never been entered");
                  #endif
                } else {
                  logDebug(
                      "receiveADERDGFaceData(...)", "receive three arrays from rank " <<
                       fromRank << " for vertex x=" << x << ", level=" << level <<
                       ", src type=" << multiscalelinkedcell::indexToString(srcCellDescriptionIndex) <<
                      ", src=" << src << ", dest=" << dest <<
                      ", counter=" << p.getFaceDataExchangeCounter(faceIndex)
                    );
                  const int numberOfFaceDof = solver->getUnknownsPerFace();
                  int receivedlQhbndIndex = DataHeap::getInstance().createData(0, numberOfFaceDof);
                  int receivedlFhbndIndex = DataHeap::getInstance().createData(0, numberOfFaceDof);
                  int receivedMinMax      = DataHeap::getInstance().createData(0, 2*solver->getNumberOfVariables());

                  assertion(DataHeap::getInstance().getData(receivedlQhbndIndex).empty());
                  assertion(DataHeap::getInstance().getData(receivedlFhbndIndex).empty());
                  assertion(DataHeap::getInstance().getData(receivedMinMax).empty());

                  // Send order: minMax,lQhbnd,lFhbnd
                  // Receive order: lFhbnd,lQhbnd,minMax
                  DataHeap::getInstance().receiveData(receivedlFhbndIndex, fromRank, x, level,
                      peano::heap::MessageType::NeighbourCommunication);
                  DataHeap::getInstance().receiveData(receivedlQhbndIndex, fromRank, x, level,
                      peano::heap::MessageType::NeighbourCommunication);
                  DataHeap::getInstance().receiveData(receivedMinMax,  fromRank, x, level,
                      peano::heap::MessageType::NeighbourCommunication);

                  logDebug(
                      "receiveADERDGFaceData(...)", "[pre] solve Riemann problem with received data." <<
                      " cellDescription=" << p.toString() <<
                      ",faceIndexForCell=" << faceIndex <<
                      ",normalOfExchangedFac=" << normalOfExchangedFace <<
                      ",x=" << x.toString() << ", level=" << level <<
                      ", counter=" << p.getFaceDataExchangeCounter(faceIndex)
                  );

                  solveRiemannProblemAtInterface(
                      p,
                      faceIndex,
                      normalOfExchangedFace,
                      receivedlQhbndIndex,
                      receivedlFhbndIndex);

                  Cell::mergeSolutionMinMaxOnFace(
                      p,
                      faceIndex,
                      DataHeap::getInstance().getData(receivedMinMax).data(),
                      DataHeap::getInstance().getData(receivedMinMax).data() + solver->getNumberOfVariables() );

                  // TODO(Dominic): If anarchic time stepping, receive the time step too.

                  DataHeap::getInstance().deleteData(receivedlQhbndIndex);
                  DataHeap::getInstance().deleteData(receivedlFhbndIndex);
                  DataHeap::getInstance().deleteData(receivedMinMax);
                }
              } else { // TODO(Dominic): just drop the data. Add to docu.
                if (srcCellDescriptionIndex == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex) {
                  #if defined(PeriodicBC)
                  assertionMsg(false, "Vasco, we have to implement this");
                  #else
                  assertionMsg(false, "should never been entered");
                  #endif
                } else {
                  assertion(p.getType()!=exahype::records::ADERDGCellDescription::Erased);
                  // TODO(Dominic): Revise this for dyn. AMR.

                  logDebug(
                      "receiveADERDGFaceData(...)", "drop three arrays from rank " <<
                      fromRank << " for vertex x=" << x << ", level=" << level <<
                      ", src type=" << multiscalelinkedcell::indexToString(srcCellDescriptionIndex) <<
                      ", src=" << src << ", dest=" << dest <<
                      ", counter=" << p.getFaceDataExchangeCounter(faceIndex)
                  );

                  DataHeap::getInstance().receiveData(fromRank, x, level,
                      peano::heap::MessageType::NeighbourCommunication);
                  DataHeap::getInstance().receiveData(fromRank, x, level,
                      peano::heap::MessageType::NeighbourCommunication);
                  DataHeap::getInstance().receiveData(fromRank, x, level,
                      peano::heap::MessageType::NeighbourCommunication);
                }
              }
            }
          }
        }
      }
    }
  }
}

void exahype::mappings::RiemannSolver::dropADERDGFaceData(
    int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    int level,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    int srcCellDescriptionIndex,
    int destCellDescriptionIndex,
    exahype::MetadataHeap::HeapEntries& receivedMetadata) {
  ADERDGCellDescriptionHeap::getInstance().isValidIndex(destCellDescriptionIndex);
  FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(destCellDescriptionIndex);

  exahype::MetadataHeap::HeapEntries::const_iterator metadataIterator=
      receivedMetadata.begin();
  const int nADERDG = metadataIterator->getU();

  if (nADERDG > 0) {
    logDebug("dropADERDGFaceData(...)","nADERDG: " << nADERDG);
    const int endOfADERDGMetadata = 2*nADERDG;
    metadataIterator += endOfADERDGMetadata;

    // TODO(Dominic): Add to docu why we need to invert the order.
    while (metadataIterator!=receivedMetadata.begin()) {
      const int neighbourTypeAsInt    = metadataIterator->getU();
      metadataIterator-=2;
      exahype::records::ADERDGCellDescription::Type neighbourType =
          static_cast<exahype::records::ADERDGCellDescription::Type>(neighbourTypeAsInt);
      logDebug("dropADERDGFaceData(...)","neighbourTypeAsInt: "    << neighbourTypeAsInt);

      if (neighbourType==exahype::records::ADERDGCellDescription::Cell ||
          neighbourType==exahype::records::ADERDGCellDescription::Ancestor ||
          neighbourType==exahype::records::ADERDGCellDescription::Descendant) {
        logDebug(
            "dropADERDGFaceData(...)", "drop three arrays from rank " <<
            fromRank << " for vertex x=" << x << ", level=" << level <<
            ", src type=" << multiscalelinkedcell::indexToString(srcCellDescriptionIndex) <<
            ", src=" << src << ", dest=" << dest
        );

        DataHeap::getInstance().receiveData(fromRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
        DataHeap::getInstance().receiveData(fromRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
        DataHeap::getInstance().receiveData(fromRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
      }
    }
  }
}

void exahype::mappings::RiemannSolver::solveRiemannProblemAtInterface(
    records::ADERDGCellDescription& cellDescription, const int faceIndexForCell,
    const int normalNonZero, const int indexOfQValues,
    const int indexOfFValues) {
  exahype::solvers::ADERDGSolver* solver =  static_cast<exahype::solvers::ADERDGSolver*>(
      exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);

  cellDescription.setRiemannSolvePerformed(faceIndexForCell, true);

  const int numberOfFaceDof = solver->getUnknownsPerFace();

  logDebug("solveRiemannProblemAtInterface(...)",
           "cell-description=" << cellDescription.toString());

  double* QL = 0;
  double* QR = 0;
  double* FL = 0;
  double* FR = 0;

  assertionEquals(DataHeap::getInstance().getData(indexOfQValues).size(),
                  static_cast<unsigned int>(numberOfFaceDof));
  assertionEquals(DataHeap::getInstance().getData(indexOfFValues).size(),
                  static_cast<unsigned int>(numberOfFaceDof));

  // @todo Doku im Header warum wir das hier brauchen,
  if (faceIndexForCell % 2 == 0) {
    QL = DataHeap::getInstance().getData(indexOfQValues).data();
    QR = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
         (faceIndexForCell * numberOfFaceDof);
    FL = DataHeap::getInstance().getData(indexOfFValues).data();
    FR = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
         (faceIndexForCell * numberOfFaceDof);
  } else {
    QR = DataHeap::getInstance().getData(indexOfQValues).data();
    QL = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
         (faceIndexForCell * numberOfFaceDof);
    FR = DataHeap::getInstance().getData(indexOfFValues).data();
    FL = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
         (faceIndexForCell * numberOfFaceDof);
  }

  // Synchronise time stepping.
  solver->synchroniseTimeStepping(cellDescription);

  solver->riemannSolver(FL, FR, QL, QR,
                        cellDescription.getCorrectorTimeStepSize(),
                        normalNonZero);

  for (int ii = 0; ii<numberOfFaceDof; ii++) {
    assertion8(std::isfinite(QR[ii]), cellDescription.toString(),
               faceIndexForCell, normalNonZero, indexOfQValues, indexOfFValues,
               ii, QR[ii], QL[ii]);
    assertion8(std::isfinite(QL[ii]), cellDescription.toString(),
               faceIndexForCell, normalNonZero, indexOfQValues, indexOfFValues,
               ii, QR[ii], QL[ii]);
    assertion8(std::isfinite(FR[ii]), cellDescription.toString(),
               faceIndexForCell, normalNonZero, indexOfQValues, indexOfFValues,
               ii, QR[ii], QL[ii]);
    assertion8(std::isfinite(FL[ii]), cellDescription.toString(),
               faceIndexForCell, normalNonZero, indexOfQValues, indexOfFValues,
               ii, QR[ii], QL[ii]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.
}



//
// Below all methods are nop.
//
//===================================



void exahype::mappings::RiemannSolver::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::RiemannSolver::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::RiemannSolver::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::RiemannSolver::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::RiemannSolver::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::RiemannSolver::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  // do nothing but please consult header documentation.
  return false;
}

void exahype::mappings::RiemannSolver::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::RiemannSolver::mergeWithMaster(
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

void exahype::mappings::RiemannSolver::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing but please consult header documentation
}

void exahype::mappings::RiemannSolver::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::RiemannSolver::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::RiemannSolver::RiemannSolver() {
  #ifdef Parallel
  _state = 0;
  #endif
}

exahype::mappings::RiemannSolver::~RiemannSolver() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::RiemannSolver::RiemannSolver(
    const RiemannSolver& masterThread) {}

void exahype::mappings::RiemannSolver::mergeWithWorkerThread(
    const RiemannSolver& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::RiemannSolver::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::RiemannSolver::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::RiemannSolver::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::RiemannSolver::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::RiemannSolver::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::RiemannSolver::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::RiemannSolver::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::RiemannSolver::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::RiemannSolver::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::RiemannSolver::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::RiemannSolver::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}
