#include "EulerFlow/adapters/FaceDataExchange.h"

peano::CommunicationSpecification
exahype::adapters::FaceDataExchange::communicationSpecification() {
  return peano::CommunicationSpecification::getMinimalSpecification() &
         exahype::mappings::FaceDataExchange::communicationSpecification()

      ;
}

peano::MappingSpecification
exahype::adapters::FaceDataExchange::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::FaceDataExchange::touchVertexLastTimeSpecification()

      ;
}

peano::MappingSpecification
exahype::adapters::FaceDataExchange::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::FaceDataExchange::
             touchVertexFirstTimeSpecification()

      ;
}

peano::MappingSpecification
exahype::adapters::FaceDataExchange::enterCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::FaceDataExchange::enterCellSpecification()

      ;
}

peano::MappingSpecification
exahype::adapters::FaceDataExchange::leaveCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::FaceDataExchange::leaveCellSpecification()

      ;
}

peano::MappingSpecification
exahype::adapters::FaceDataExchange::ascendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::FaceDataExchange::ascendSpecification()

      ;
}

peano::MappingSpecification
exahype::adapters::FaceDataExchange::descendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::FaceDataExchange::descendSpecification()

      ;
}

exahype::adapters::FaceDataExchange::FaceDataExchange() {}

exahype::adapters::FaceDataExchange::~FaceDataExchange() {}

#if defined(SharedMemoryParallelisation)
exahype::adapters::FaceDataExchange::FaceDataExchange(
    const FaceDataExchange& masterThread)
    : _map2FaceDataExchange(masterThread._map2FaceDataExchange)

{}

void exahype::adapters::FaceDataExchange::mergeWithWorkerThread(
    const FaceDataExchange& workerThread) {
  _map2FaceDataExchange.mergeWithWorkerThread(
      workerThread._map2FaceDataExchange);
}
#endif

void exahype::adapters::FaceDataExchange::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2FaceDataExchange.createHangingVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::FaceDataExchange::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2FaceDataExchange.destroyHangingVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::FaceDataExchange::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2FaceDataExchange.createInnerVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::FaceDataExchange::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2FaceDataExchange.createBoundaryVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::FaceDataExchange::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2FaceDataExchange.destroyVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::FaceDataExchange::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2FaceDataExchange.createCell(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell);
}

void exahype::adapters::FaceDataExchange::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2FaceDataExchange.destroyCell(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell);
}

#ifdef Parallel
void exahype::adapters::FaceDataExchange::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  _map2FaceDataExchange.mergeWithNeighbour(vertex, neighbour, fromRank,
                                           fineGridX, fineGridH, level);
}

void exahype::adapters::FaceDataExchange::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2FaceDataExchange.prepareSendToNeighbour(vertex, toRank, x, h, level);
}

void exahype::adapters::FaceDataExchange::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2FaceDataExchange.prepareCopyToRemoteNode(localVertex, toRank, x, h,
                                                level);
}

void exahype::adapters::FaceDataExchange::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2FaceDataExchange.prepareCopyToRemoteNode(localCell, toRank, x, h, level);
}

void exahype::adapters::FaceDataExchange::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2FaceDataExchange.mergeWithRemoteDataDueToForkOrJoin(
      localVertex, masterOrWorkerVertex, fromRank, x, h, level);
}

void exahype::adapters::FaceDataExchange::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2FaceDataExchange.mergeWithRemoteDataDueToForkOrJoin(
      localCell, masterOrWorkerCell, fromRank, x, h, level);
}

bool exahype::adapters::FaceDataExchange::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  bool result = false;
  result |= _map2FaceDataExchange.prepareSendToWorker(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell, worker);

  return result;
}

void exahype::adapters::FaceDataExchange::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2FaceDataExchange.prepareSendToMaster(
      localCell, vertices, verticesEnumerator, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell);
}

void exahype::adapters::FaceDataExchange::mergeWithMaster(
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
  _map2FaceDataExchange.mergeWithMaster(
      workerGridCell, workerGridVertices, workerEnumerator, fineGridCell,
      fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell,
      worker, workerState, masterState);
}

void exahype::adapters::FaceDataExchange::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2FaceDataExchange.receiveDataFromMaster(
      receivedCell, receivedVertices, receivedVerticesEnumerator,
      receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator,
      receivedCoarseGridCell, workersCoarseGridVertices,
      workersCoarseGridVerticesEnumerator, workersCoarseGridCell,
      fineGridPositionOfCell);
}

void exahype::adapters::FaceDataExchange::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  _map2FaceDataExchange.mergeWithWorker(localCell, receivedMasterCell,
                                        cellCentre, cellSize, level);
}

void exahype::adapters::FaceDataExchange::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2FaceDataExchange.mergeWithWorker(localVertex, receivedMasterVertex, x, h,
                                        level);
}
#endif

void exahype::adapters::FaceDataExchange::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2FaceDataExchange.touchVertexFirstTime(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::FaceDataExchange::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2FaceDataExchange.touchVertexLastTime(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::FaceDataExchange::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2FaceDataExchange.enterCell(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell);
}

void exahype::adapters::FaceDataExchange::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2FaceDataExchange.leaveCell(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell);
}

void exahype::adapters::FaceDataExchange::beginIteration(
    exahype::State& solverState) {
  _map2FaceDataExchange.beginIteration(solverState);
}

void exahype::adapters::FaceDataExchange::endIteration(
    exahype::State& solverState) {
  _map2FaceDataExchange.endIteration(solverState);
}

void exahype::adapters::FaceDataExchange::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  _map2FaceDataExchange.descend(fineGridCells, fineGridVertices,
                                fineGridVerticesEnumerator, coarseGridVertices,
                                coarseGridVerticesEnumerator, coarseGridCell);
}

void exahype::adapters::FaceDataExchange::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  _map2FaceDataExchange.ascend(fineGridCells, fineGridVertices,
                               fineGridVerticesEnumerator, coarseGridVertices,
                               coarseGridVerticesEnumerator, coarseGridCell);
}
