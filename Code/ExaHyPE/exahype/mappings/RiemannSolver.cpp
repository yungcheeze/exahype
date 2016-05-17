#include "exahype/mappings/RiemannSolver.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/utils/Globals.h"

#include "tarch/multicore/Lock.h"
#include "tarch/multicore/Loop.h"

#include "exahype/solvers/Solver.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::RiemannSolver::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::RiemannSolver::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::RiemannSolver::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::RiemannSolver::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::RiemannSolver::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::RiemannSolver::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::RiemannSolver::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::RiemannSolver::_log(
    "exahype::mappings::RiemannSolver");

exahype::mappings::RiemannSolver::RiemannSolver() {
  // do nothing
}

exahype::mappings::RiemannSolver::~RiemannSolver() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::RiemannSolver::RiemannSolver(
    const RiemannSolver& masterThread)
    : _localState(masterThread._localState) {}

void exahype::mappings::RiemannSolver::mergeWithWorkerThread(
    const RiemannSolver& workerThread) {
  // do nothing
}
#endif

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

#ifdef Parallel
void exahype::mappings::RiemannSolver::mergeWithNeighbour(
  exahype::Vertex&                              vertex,
  const exahype::Vertex&                        neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS, double>&  fineGridX,
  const tarch::la::Vector<DIMENSIONS, double>&  fineGridH,
  int                                           level) {
  #if !defined(PeriodicBC)
  if (vertex.isBoundary()) return;
  #endif

  tarch::la::Vector<TWO_POWER_D, int>& adjacentADERDGCellDescriptionsIndices =
      vertex.getADERDGCellDescriptionsIndex();

  dfor2(myDest)
    dfor2(mySrc)
/*
 * doku
 */
      tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
      tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;

      int destScalar = TWO_POWER_D-myDestScalar-1;
      int srcScalar  = TWO_POWER_D-mySrcScalar-1;

      if (
        vertex.getAdjacentRanks()(destScalar) ==
        tarch::parallel::Node::getInstance().getRank() &&
        vertex.getAdjacentRanks()(srcScalar) == fromRank &&
        tarch::la::countEqualEntries(dest, src) == 1  // we are solely exchanging faces
      ) {
        const int destCellDescriptionIndex =
          adjacentADERDGCellDescriptionsIndices(destScalar);

        if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(destCellDescriptionIndex)) {
          std::vector<records::ADERDGCellDescription>& cellDescriptions =
              ADERDGCellDescriptionHeap::getInstance().getData(
                  destCellDescriptionIndex);

          for (
            int currentSolver = 0;
            currentSolver < static_cast<int>(cellDescriptions.size());
            currentSolver++
          ) {
            if (cellDescriptions[currentSolver].getType() == exahype::records::ADERDGCellDescription::Cell) {
              exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[cellDescriptions[currentSolver].getSolverNumber()];

              const int numberOfFaceDof = solver->getUnknownsPerFace();
              const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
              assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);

              assertion(DataHeap::getInstance().isValidIndex(cellDescriptions[currentSolver].getExtrapolatedPredictor()));
              assertion(DataHeap::getInstance().isValidIndex(cellDescriptions[currentSolver].getFluctuation()));

              if (adjacentADERDGCellDescriptionsIndices(destScalar) == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex) {
                #if defined(PeriodicBC)
                assertionMsg(false, "Vasco, we have to implement this");
                #else
                assertionMsg(false, "should never been entered");
                #endif
              } else {
                logDebug(
                  "mergeWithNeighbour(...)",
                  "receive two arrays from rank " << fromRank << " for vertex " << vertex.toString()
                  << ", src type=" << multiscalelinkedcell::indexToString(adjacentADERDGCellDescriptionsIndices(srcScalar)) <<
                  ", src=" << src << ", dest=" << dest
                );

                int receivedlQhbndIndex = DataHeap::getInstance().createData(0,numberOfFaceDof);
                int receivedlFhbndIndex = DataHeap::getInstance().createData(0,numberOfFaceDof);

                assertion( DataHeap::getInstance().getData(receivedlQhbndIndex).empty() );
                assertion( DataHeap::getInstance().getData(receivedlFhbndIndex).empty() );

                // @todo Reihenfolge dokumentieren! Auch umgedreht hier
                DataHeap::getInstance().receiveData(
                    receivedlFhbndIndex,
                    fromRank, fineGridX, level,
                    peano::heap::MessageType::NeighbourCommunication);
                DataHeap::getInstance().receiveData(
                    receivedlQhbndIndex,
                    fromRank, fineGridX, level,
                    peano::heap::MessageType::NeighbourCommunication);

                int faceIndexForCell = -1;
                     if ((normalOfExchangedFace==0) & (src(normalOfExchangedFace)<dest(normalOfExchangedFace))) {
                  faceIndexForCell = 0;
                }
                else if ((normalOfExchangedFace==0) & (src(normalOfExchangedFace)>dest(normalOfExchangedFace))) {
                  faceIndexForCell = 1;
                }
                else if ((normalOfExchangedFace==1) & (src(normalOfExchangedFace)<dest(normalOfExchangedFace))) {
                  faceIndexForCell = 2;
                }
                else if ((normalOfExchangedFace==1) & (src(normalOfExchangedFace)>dest(normalOfExchangedFace))) {
                  faceIndexForCell = 3;
                }
                else if ((normalOfExchangedFace==2) & (src(normalOfExchangedFace)<dest(normalOfExchangedFace))) {
                  faceIndexForCell = 4;
                }
                else if ((normalOfExchangedFace==2) & (src(normalOfExchangedFace)>dest(normalOfExchangedFace))) {
                  faceIndexForCell = 5;
                }
                else {
                  assertionMsg( false, "should not be entered" );
                }


                if (!cellDescriptions[currentSolver].getRiemannSolvePerformed(faceIndexForCell)) {
                    logDebug(
                      "mergeWithNeighbour(...)",
                      "solve Riemann problem with received data." <<
                      " cellDescription=" << cellDescriptions[currentSolver].toString() <<
                      ",faceIndexForCell=" << faceIndexForCell <<
                      ",normalOfExchangedFac=" << normalOfExchangedFace <<
                      ",vertex=" << vertex.toString()
                    );

                  solveRiemannProblemAtInterface(
                    cellDescriptions[currentSolver],
                    faceIndexForCell,
                    normalOfExchangedFace,
                    receivedlQhbndIndex,
                    receivedlFhbndIndex);
                }

                DataHeap::getInstance().deleteData(receivedlQhbndIndex);
                DataHeap::getInstance().deleteData(receivedlFhbndIndex);
              }
            } else {
              assertionMsg(false, "Dominic, please implement");
            }
          }
        }
      }
  enddforx enddforx
}

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

void exahype::mappings::RiemannSolver::touchVertexFirstTime(
    exahype::Vertex&                              fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>&  fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>&  fineGridH,
    exahype::Vertex* const                        coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&                                coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>&     fineGridPositionOfVertex) {
  logTraceInWith6Arguments("touchVertexFirstTime(...)", fineGridVertex,
                           fineGridX, fineGridH,
                           coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  tarch::la::Vector<TWO_POWER_D, int>& adjacentCellDescriptionsIndices =
      fineGridVertex.getADERDGCellDescriptionsIndex();
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
    solveRiemannProblemAtInterface(
      adjacentCellDescriptionsIndices[cellIndicesLeft[i]],
      adjacentCellDescriptionsIndices[cellIndicesRight[i]],
      EXAHYPE_FACE_RIGHT, EXAHYPE_FACE_LEFT, 0);

    solveRiemannProblemAtInterface(
      adjacentCellDescriptionsIndices[cellIndicesFront[i]],
      adjacentCellDescriptionsIndices[cellIndicesBack[i]],
      EXAHYPE_FACE_BACK, EXAHYPE_FACE_FRONT, 1);

    #if DIMENSIONS == 3
    solveRiemannProblemAtInterface(
      adjacentCellDescriptionsIndices[cellIndicesBottom[i]],
      adjacentCellDescriptionsIndices[cellIndicesTop[i]],
      EXAHYPE_FACE_TOP, EXAHYPE_FACE_BOTTOM, 2);
    #endif
  }
  logTraceOutWith1Argument("touchVertexFirstTime(...)", fineGridVertex);
}

void exahype::mappings::RiemannSolver::solveRiemannProblemAtInterface(
    const int cellDescriptionIndexOfLeftCell,
    const int cellDescriptionIndexOfRightCell,
    const int faceIndexForLeftCell,
    const int faceIndexForRightCell,
    const int normalNonZero) {
  // Only continue if this is an internal face, i.e.,
  // both cell description indices are valid
  if (
    ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionIndexOfLeftCell) &&
    ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionIndexOfRightCell)
  ) {
    logDebug("touchVertexLastTime(...)::solveRiemannProblemAtInterface(...)",
             "Performing Riemann solve. "
                 << "faceIndexForLeftCell:" << faceIndexForLeftCell
                 << " faceIndexForRightCell:" << faceIndexForRightCell
                 << " indexOfLeftCell:"
                 << adjacentADERDGCellDescriptionsIndices[indexOfLeftCell]
                 << " indexOfRightCell:"
                 << adjacentADERDGCellDescriptionsIndices[indexOfRightCell]);

    std::vector<records::ADERDGCellDescription>& cellDescriptionsOfLeftCell =
        ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionIndexOfLeftCell);
    std::vector<records::ADERDGCellDescription>& cellDescriptionsOfRightCell =
        ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionIndexOfRightCell);

    // @todo 08/02/16:Dominic Etienne Charrier
    // Assumes that the both cells hold the same number of cell descriptions
    assertion2(
        cellDescriptionsOfLeftCell.size() == cellDescriptionsOfRightCell.size(),
        cellDescriptionsOfLeftCell.size(),
        cellDescriptionsOfRightCell.size());

    const int numberOfADERDGCellDescriptions = static_cast<int>(
        ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionIndexOfLeftCell).size());
    const peano::datatraversal::autotuning::MethodTrace methodTrace =
        peano::datatraversal::autotuning::UserDefined4;
    const int grainSize =
        peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
            numberOfADERDGCellDescriptions, methodTrace);

    pfor(i, 0, numberOfADERDGCellDescriptions,grainSize)
      if (cellDescriptionsOfLeftCell[i].getType() ==
                              exahype::records::ADERDGCellDescription::Cell ||
                          cellDescriptionsOfRightCell[i].getType() ==
                              exahype::records::ADERDGCellDescription::Cell) {
        exahype::solvers::Solver* solver =
            exahype::solvers::RegisteredSolvers[cellDescriptionsOfLeftCell[i]
                                                    .getSolverNumber()];

        assertionEquals4(
            cellDescriptionsOfLeftCell[i].getRiemannSolvePerformed(faceIndexForLeftCell),
            cellDescriptionsOfRightCell[i].getRiemannSolvePerformed(faceIndexForRightCell),
            faceIndexForLeftCell, faceIndexForRightCell,
            cellDescriptionsOfLeftCell[i].toString(), cellDescriptionsOfRightCell[i].toString());

        if (!cellDescriptionsOfLeftCell[i].getRiemannSolvePerformed(faceIndexForLeftCell)) {
          cellDescriptionsOfLeftCell[i].setRiemannSolvePerformed( faceIndexForLeftCell,  true);
          cellDescriptionsOfRightCell[i].setRiemannSolvePerformed(faceIndexForRightCell, true);

          const int numberOfFaceDof = solver->getUnknownsPerFace();

          double* QL =
              DataHeap::getInstance()
                  .getData(cellDescriptionsOfLeftCell[i].getExtrapolatedPredictor())
                  .data() +
              (faceIndexForLeftCell * numberOfFaceDof);
          double* QR =
              DataHeap::getInstance()
                  .getData(cellDescriptionsOfRightCell[i].getExtrapolatedPredictor())
                  .data() +
              (faceIndexForRightCell * numberOfFaceDof);
          double* FL = DataHeap::getInstance()
                           .getData(cellDescriptionsOfLeftCell[i].getFluctuation())
                           .data() +
                       (faceIndexForLeftCell * numberOfFaceDof);
          double* FR = DataHeap::getInstance()
                           .getData(cellDescriptionsOfRightCell[i].getFluctuation())
                           .data() +
                       (faceIndexForRightCell * numberOfFaceDof);

          solver->synchroniseTimeStepping(cellDescriptionsOfLeftCell[i]);
          solver->synchroniseTimeStepping(cellDescriptionsOfRightCell[i]);

          logDebug("touchVertexLastTime(...)::debug::before::QL[0]*", QL[0]);
          logDebug("touchVertexLastTime(...)::debug::before::QR[0]*", QR[0]);
          logDebug("touchVertexLastTime(...)::debug::before::FL[0]", FL[0]);
          logDebug("touchVertexLastTime(...)::debug::before::FR[0]", FR[0]);

          solver->riemannSolver(
              FL, FR, QL, QR,
              std::min(cellDescriptionsOfLeftCell[i].getCorrectorTimeStepSize(),
                       cellDescriptionsOfRightCell[i].getCorrectorTimeStepSize()),
              normalNonZero);

          logDebug("touchVertexLastTime(...)::debug::after::QL[0]*", QL[0]);
          logDebug("touchVertexLastTime(...)::debug::after::QR[0]*", QR[0]);
          logDebug("touchVertexLastTime(...)::debug::after::FL[0]", FL[0]);
          logDebug("touchVertexLastTime(...)::debug::after::FR[0]", FR[0]);
        }
      }
    endpfor

    peano::datatraversal::autotuning::Oracle::getInstance().parallelSectionHasTerminated(methodTrace);
  } else if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionIndexOfLeftCell) !=
      ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionIndexOfRightCell)) { // XOR
    const int cellDescriptionsIndex =
        (ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionIndexOfLeftCell)) ?
            cellDescriptionIndexOfLeftCell :
            cellDescriptionIndexOfRightCell;
    const int faceIndex =
            (cellDescriptionsIndex==cellDescriptionIndexOfLeftCell) ?
                faceIndexForLeftCell :
                faceIndexForRightCell;

    assertion(ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndex));

    std::vector<records::ADERDGCellDescription>& cellDescriptions =
        ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex);

    const int numberOfADERDGCellDescriptions = static_cast<int>(cellDescriptions.size());
    const peano::datatraversal::autotuning::MethodTrace methodTrace =
        peano::datatraversal::autotuning::UserDefined0;
    const int grainSize =
        peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
            numberOfADERDGCellDescriptions, methodTrace);

    pfor(i, 0, numberOfADERDGCellDescriptions,grainSize)
      if (!cellDescriptions[i].getRiemannSolvePerformed(faceIndex)) {
        cellDescriptions[i].setRiemannSolvePerformed( faceIndex,  true);
        switch (cellDescriptions[i].getType()) {
        case exahype::records::ADERDGCellDescription::Cell:
          switch (cellDescriptions[i].getRefinementEvent()) {
          case exahype::records::ADERDGCellDescription::None:
            applyBoundaryConditions(cellDescriptions[i],faceIndex,normalNonZero);
            break;
          default:
            break;
          }
          break;
          default:
            break;
        }
      }
    endpfor
  }
}


void exahype::mappings::RiemannSolver::solveRiemannProblemAtInterface(
    records::ADERDGCellDescription& cellDescription,
    const int faceIndexForCell,
    const int normalNonZero,
    const int indexOfQValues,
    const int indexOfFValues
) {
  exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()];

  cellDescription.setRiemannSolvePerformed( faceIndexForCell, true);

  const int numberOfFaceDof = solver->getUnknownsPerFace();

  solver->synchroniseTimeStepping(cellDescription);
  logDebug( "solveRiemannProblemAtInterface(...)", "cell-description=" << cellDescription.toString() );

  double* QL = 0;
  double* QR = 0;
  double* FL = 0;
  double* FR = 0;

  assertionEquals(DataHeap::getInstance().getData(indexOfQValues).size(),static_cast<unsigned int>(numberOfFaceDof));
  assertionEquals(DataHeap::getInstance().getData(indexOfFValues).size(),static_cast<unsigned int>(numberOfFaceDof));

  // @todo Doku im Header warum wir das hier brauchen,
  if (faceIndexForCell%2==0) {
    QL = DataHeap::getInstance().getData(indexOfQValues).data();
    QR = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data()
                 + (faceIndexForCell * numberOfFaceDof);
    FL = DataHeap::getInstance().getData(indexOfFValues).data();
    FR = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data()
                 + (faceIndexForCell * numberOfFaceDof);
  }
  else {
    QR = DataHeap::getInstance().getData(indexOfQValues).data();
    QL = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data()
                  + (faceIndexForCell * numberOfFaceDof);
    FR = DataHeap::getInstance().getData(indexOfFValues).data();
    FL = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data()
                   + (faceIndexForCell * numberOfFaceDof);
  }

  solver->riemannSolver(
    FL, FR, QL, QR,
    cellDescription.getCorrectorTimeStepSize(),
    normalNonZero
    );


  for (int i=0; i<numberOfFaceDof; i++) {
      // @todo Groessen raus nehmen vermutlich
    assertion8( QR[i]==QR[i] && QR[i]<1e10, cellDescription.toString(), faceIndexForCell, normalNonZero, indexOfQValues, indexOfFValues, i, QR[i], QL[i] );
    assertion8( QL[i]==QL[i] && QL[i]<1e10, cellDescription.toString(), faceIndexForCell, normalNonZero, indexOfQValues, indexOfFValues, i, QR[i], QL[i] );
    assertion8( FR[i]==FR[i] && FR[i]<1e10, cellDescription.toString(), faceIndexForCell, normalNonZero, indexOfQValues, indexOfFValues, i, QR[i], QL[i] );
    assertion8( FL[i]==FL[i] && FL[i]<1e10, cellDescription.toString(), faceIndexForCell, normalNonZero, indexOfQValues, indexOfFValues, i, QR[i], QL[i] );
  }
}

void exahype::mappings::RiemannSolver::applyBoundaryConditions(
    records::ADERDGCellDescription& cellDescription,
    const int faceIndex,
    const int normalNonZero) {
  exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()];

  cellDescription.setRiemannSolvePerformed( faceIndex, true);

  const int numberOfFaceDof =
      solver->getUnknownsPerFace();

  double* lQhbnd =
      DataHeap::getInstance().getData(
          cellDescription.getExtrapolatedPredictor()).data() +
          (faceIndex * numberOfFaceDof);
  double* lFhbnd =
      DataHeap::getInstance().getData(
          cellDescription.getFluctuation()).data() +
          (faceIndex * numberOfFaceDof);

  solver->synchroniseTimeStepping(cellDescription);
  assertionEquals(lQhbnd[0],lQhbnd[0]); // assert no nan
  assertionEquals(lFhbnd[0],lFhbnd[0]); // assert no nan

  solver->riemannSolver(
      lFhbnd, lFhbnd, lQhbnd, lQhbnd,
      cellDescription.getCorrectorTimeStepSize(),
      normalNonZero);

  assertionEquals(lFhbnd[0],lFhbnd[0]); // assert no nan
}


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

void exahype::mappings::RiemannSolver::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
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

void exahype::mappings::RiemannSolver::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  _localState = solverState;

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::RiemannSolver::endIteration(
    exahype::State& solverState) {
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
