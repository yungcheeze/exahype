#include "exahype/adapters/Plot.h"

peano::CommunicationSpecification
exahype::adapters::Plot::communicationSpecification() {
  return peano::CommunicationSpecification::getMinimalSpecification() &
         exahype::mappings::Plot::communicationSpecification()

      ;
}

peano::MappingSpecification
exahype::adapters::Plot::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::Plot::touchVertexLastTimeSpecification()

      ;
}

peano::MappingSpecification
exahype::adapters::Plot::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::Plot::touchVertexFirstTimeSpecification()

      ;
}

peano::MappingSpecification exahype::adapters::Plot::enterCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::Plot::enterCellSpecification()

      ;
}

peano::MappingSpecification exahype::adapters::Plot::leaveCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::Plot::leaveCellSpecification()

      ;
}

peano::MappingSpecification exahype::adapters::Plot::ascendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::Plot::ascendSpecification()

      ;
}

peano::MappingSpecification exahype::adapters::Plot::descendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::Plot::descendSpecification()

      ;
}

exahype::adapters::Plot::Plot() {}

exahype::adapters::Plot::~Plot() {}

#if defined(SharedMemoryParallelisation)
exahype::adapters::Plot::Plot(const Plot& masterThread)
    : _map2Plot(masterThread._map2Plot)

{}

void exahype::adapters::Plot::mergeWithWorkerThread(const Plot& workerThread) {
  _map2Plot.mergeWithWorkerThread(workerThread._map2Plot);
}
#endif

void exahype::adapters::Plot::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2Plot.createHangingVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::Plot::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2Plot.destroyHangingVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::Plot::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2Plot.createInnerVertex(fineGridVertex, fineGridX, fineGridH,
                              coarseGridVertices, coarseGridVerticesEnumerator,
                              coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::Plot::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2Plot.createBoundaryVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::Plot::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2Plot.destroyVertex(fineGridVertex, fineGridX, fineGridH,
                          coarseGridVertices, coarseGridVerticesEnumerator,
                          coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::Plot::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2Plot.createCell(fineGridCell, fineGridVertices,
                       fineGridVerticesEnumerator, coarseGridVertices,
                       coarseGridVerticesEnumerator, coarseGridCell,
                       fineGridPositionOfCell);
}

void exahype::adapters::Plot::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2Plot.destroyCell(fineGridCell, fineGridVertices,
                        fineGridVerticesEnumerator, coarseGridVertices,
                        coarseGridVerticesEnumerator, coarseGridCell,
                        fineGridPositionOfCell);
}

#ifdef Parallel
void exahype::adapters::Plot::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  _map2Plot.mergeWithNeighbour(vertex, neighbour, fromRank, fineGridX,
                               fineGridH, level);
}

void exahype::adapters::Plot::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2Plot.prepareSendToNeighbour(vertex, toRank, x, h, level);
}

void exahype::adapters::Plot::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2Plot.prepareCopyToRemoteNode(localVertex, toRank, x, h, level);
}

void exahype::adapters::Plot::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2Plot.prepareCopyToRemoteNode(localCell, toRank, x, h, level);
}

void exahype::adapters::Plot::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2Plot.mergeWithRemoteDataDueToForkOrJoin(
      localVertex, masterOrWorkerVertex, fromRank, x, h, level);
}

void exahype::adapters::Plot::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2Plot.mergeWithRemoteDataDueToForkOrJoin(localCell, masterOrWorkerCell,
                                               fromRank, x, h, level);
}

bool exahype::adapters::Plot::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  bool result = false;
  result |= _map2Plot.prepareSendToWorker(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell, worker);

  return result;
}

void exahype::adapters::Plot::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2Plot.prepareSendToMaster(
      localCell, vertices, verticesEnumerator, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell);
}

void exahype::adapters::Plot::mergeWithMaster(
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
  _map2Plot.mergeWithMaster(
      workerGridCell, workerGridVertices, workerEnumerator, fineGridCell,
      fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell,
      worker, workerState, masterState);
}

void exahype::adapters::Plot::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2Plot.receiveDataFromMaster(
      receivedCell, receivedVertices, receivedVerticesEnumerator,
      receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator,
      receivedCoarseGridCell, workersCoarseGridVertices,
      workersCoarseGridVerticesEnumerator, workersCoarseGridCell,
      fineGridPositionOfCell);
}

void exahype::adapters::Plot::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  _map2Plot.mergeWithWorker(localCell, receivedMasterCell, cellCentre, cellSize,
                            level);
}

void exahype::adapters::Plot::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2Plot.mergeWithWorker(localVertex, receivedMasterVertex, x, h, level);
}
#endif

void exahype::adapters::Plot::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2Plot.touchVertexFirstTime(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::Plot::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2Plot.touchVertexLastTime(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::Plot::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2Plot.enterCell(fineGridCell, fineGridVertices,
                      fineGridVerticesEnumerator, coarseGridVertices,
                      coarseGridVerticesEnumerator, coarseGridCell,
                      fineGridPositionOfCell);
}

void exahype::adapters::Plot::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2Plot.leaveCell(fineGridCell, fineGridVertices,
                      fineGridVerticesEnumerator, coarseGridVertices,
                      coarseGridVerticesEnumerator, coarseGridCell,
                      fineGridPositionOfCell);
}

void exahype::adapters::Plot::beginIteration(exahype::State& solverState) {
  _map2Plot.beginIteration(solverState);
}

void exahype::adapters::Plot::endIteration(exahype::State& solverState) {
  _map2Plot.endIteration(solverState);
}

void exahype::adapters::Plot::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  _map2Plot.descend(fineGridCells, fineGridVertices, fineGridVerticesEnumerator,
                    coarseGridVertices, coarseGridVerticesEnumerator,
                    coarseGridCell);
}

void exahype::adapters::Plot::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  _map2Plot.ascend(fineGridCells, fineGridVertices, fineGridVerticesEnumerator,
                   coarseGridVertices, coarseGridVerticesEnumerator,
                   coarseGridCell);
}
