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
  if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndexOfLeftCell)
      &&
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

          for(int ii=0; ii<numberOfFaceDof; ++ii) {
            assertion(std::isfinite(FL[ii]));
            assertion(std::isfinite(FR[ii]));
          }  // Dead code elimination will get rid of this loop if Asserts flag is not set.
        }
      }
    }
    endpfor peano::datatraversal::autotuning::Oracle::getInstance()
        .parallelSectionHasTerminated(methodTrace);
  } else if ((ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndexOfLeftCell) &&
      cellDescriptionsIndexOfRightCell == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex) ||
          (ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndexOfRightCell) &&
              cellDescriptionsIndexOfLeftCell == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex)) {
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
  }
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

  tarch::la::Vector<DIMENSIONS,double> cellCentre = cellDescription.getSize();
  cellCentre*=0.5;
  cellCentre += cellDescription.getOffset();

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

void exahype::mappings::RiemannSolver::enterCell(
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
        peano::datatraversal::autotuning::UserDefined3;
    const int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
            numberOfADERDGCellDescriptions, methodTrace);
    pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
      records::ADERDGCellDescription& p =
          ADERDGCellDescriptionHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[i];

      std::bitset<DIMENSIONS_TIMES_TWO> riemannSolvePerformed; // All bits are initialised to 'off'.
      p.setRiemannSolvePerformed(riemannSolvePerformed);
    endpfor peano::datatraversal::autotuning::Oracle::getInstance()
    .parallelSectionHasTerminated(methodTrace);
  }

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::RiemannSolver::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

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
void exahype::mappings::RiemannSolver::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
#if !defined(PeriodicBC)
  if (vertex.isBoundary()) return;
#endif

  tarch::la::Vector<TWO_POWER_D, int>& adjacentADERDGCellDescriptionsIndices =
      vertex.getCellDescriptionsIndex();

  dfor2(myDest)
    dfor2(mySrc)
      tarch::la::Vector<DIMENSIONS, int>
          dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
      tarch::la::Vector<DIMENSIONS, int> src =
          tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;

      int destScalar = TWO_POWER_D - myDestScalar - 1;
      int srcScalar = TWO_POWER_D - mySrcScalar - 1;

      if (vertex.getAdjacentRanks()(destScalar)   == tarch::parallel::Node::getInstance().getRank() &&
          vertex.getAdjacentRanks()(srcScalar)    == fromRank &&
          tarch::la::countEqualEntries(dest, src) == 1) {  // we are solely exchanging faces
        const int destCellDescriptionIndex =
            adjacentADERDGCellDescriptionsIndices(destScalar);

        if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(destCellDescriptionIndex)) {
          std::vector<records::ADERDGCellDescription>& cellDescriptions =
              ADERDGCellDescriptionHeap::getInstance().getData(
                  destCellDescriptionIndex);

          for (int currentSolver = 0;
               currentSolver < static_cast<int>(cellDescriptions.size());
               currentSolver++) {
            if (cellDescriptions[currentSolver].getType() == exahype::records::ADERDGCellDescription::Cell) {
              exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
                  exahype::solvers::RegisteredSolvers
                  [cellDescriptions[currentSolver].getSolverNumber()]);

              const int numberOfFaceDof = solver->getUnknownsPerFace();
              const int normalOfExchangedFace =
                  tarch::la::equalsReturnIndex(src, dest);
              assertion(normalOfExchangedFace >= 0 &&
                        normalOfExchangedFace < DIMENSIONS);

              assertion(DataHeap::getInstance().isValidIndex(
                  cellDescriptions[currentSolver].getExtrapolatedPredictor()));
              assertion(DataHeap::getInstance().isValidIndex(
                  cellDescriptions[currentSolver].getFluctuation()));

              if (adjacentADERDGCellDescriptionsIndices(destScalar) == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex) {
                #if defined(PeriodicBC)
                assertionMsg(false, "Vasco, we have to implement this");
                #else
                assertionMsg(false, "should never been entered");
                #endif
              } else {
                logDebug("mergeWithNeighbour(...)", "receive two arrays from rank "
                        << fromRank << " for vertex " << vertex.toString()
                        << ", src type=" << multiscalelinkedcell::indexToString(adjacentADERDGCellDescriptionsIndices(srcScalar))
                        << ", src=" << src << ", dest=" << dest);

                int receivedlQhbndIndex = DataHeap::getInstance().createData(0, numberOfFaceDof);
                int receivedlFhbndIndex = DataHeap::getInstance().createData(0, numberOfFaceDof);
                int receivedMinMax      = DataHeap::getInstance().createData(0, 2);

                assertion(DataHeap::getInstance().getData(receivedlQhbndIndex).empty());
                assertion(DataHeap::getInstance().getData(receivedlFhbndIndex).empty());

                DataHeap::getInstance().receiveData(receivedlFhbndIndex, fromRank, fineGridX, level,
                    peano::heap::MessageType::NeighbourCommunication);
                DataHeap::getInstance().receiveData(receivedlQhbndIndex, fromRank, fineGridX, level,
                    peano::heap::MessageType::NeighbourCommunication);
                DataHeap::getInstance().receiveData(receivedMinMax,  fromRank, fineGridX, level,
                    peano::heap::MessageType::NeighbourCommunication);

                int faceIndexForCell = -1;
                if ((normalOfExchangedFace == 0) &
                    (src(normalOfExchangedFace) < dest(normalOfExchangedFace))) {
                  faceIndexForCell = 0;
                } else if ((normalOfExchangedFace == 0) &
                           (src(normalOfExchangedFace) > dest(normalOfExchangedFace))) {
                  faceIndexForCell = 1;
                } else if ((normalOfExchangedFace == 1) &
                           (src(normalOfExchangedFace) < dest(normalOfExchangedFace))) {
                  faceIndexForCell = 2;
                } else if ((normalOfExchangedFace == 1) &
                           (src(normalOfExchangedFace) > dest(normalOfExchangedFace))) {
                  faceIndexForCell = 3;
                } else if ((normalOfExchangedFace == 2) &
                           (src(normalOfExchangedFace) < dest(normalOfExchangedFace))) {
                  faceIndexForCell = 4;
                } else if ((normalOfExchangedFace == 2) &
                           (src(normalOfExchangedFace) > dest(normalOfExchangedFace))) {
                  faceIndexForCell = 5;
                } else {
                  assertionMsg(false, "should not be entered");
                }

                if (!cellDescriptions[currentSolver].getRiemannSolvePerformed(faceIndexForCell)) {
                  logDebug("mergeWithNeighbour(...)", "solve Riemann problem with received data."
                               << " cellDescription=" << cellDescriptions[currentSolver].toString()
                               << ",faceIndexForCell=" << faceIndexForCell
                               << ",normalOfExchangedFac=" << normalOfExchangedFace
                               << ",vertex=" << vertex.toString());

                  solveRiemannProblemAtInterface(
                      cellDescriptions[currentSolver],
                      faceIndexForCell,
                      normalOfExchangedFace,
                      receivedlQhbndIndex,
                      receivedlFhbndIndex);

                  Cell::mergeSolutionMinMaxOnFace(
                      cellDescriptions[currentSolver],
                      faceIndexForCell,
                      DataHeap::getInstance().getData(receivedMinMax)[0],
                      DataHeap::getInstance().getData(receivedMinMax)[1]);
                }

                DataHeap::getInstance().deleteData(receivedlQhbndIndex);
                DataHeap::getInstance().deleteData(receivedlFhbndIndex);
                DataHeap::getInstance().deleteData(receivedMinMax);
              }
            } else {
              assertionMsg(false, "Dominic, please implement");
            }
          }
        }
      }
    enddforx
  enddforx
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

  return true;
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
  // do nothing
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
