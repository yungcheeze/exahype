#include "EulerFlow/adapters/PatchInitialisationAndExport.h"

peano::CommunicationSpecification
exahype::adapters::PatchInitialisationAndExport::communicationSpecification() {
  return peano::CommunicationSpecification::getMinimalSpecification() &
         exahype::mappings::PatchInitialisation::communicationSpecification() &
         exahype::adapters::
             PatchInitialisationAndExport2MultiscaleLinkedCell_1::
                 communicationSpecification() &
         exahype::mappings::VTKExport::communicationSpecification()

      ;
}

peano::MappingSpecification exahype::adapters::PatchInitialisationAndExport::
    touchVertexLastTimeSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::PatchInitialisation::
             touchVertexLastTimeSpecification() &
         exahype::adapters::
             PatchInitialisationAndExport2MultiscaleLinkedCell_1::
                 touchVertexLastTimeSpecification() &
         exahype::mappings::VTKExport::touchVertexLastTimeSpecification()

      ;
}

peano::MappingSpecification exahype::adapters::PatchInitialisationAndExport::
    touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::PatchInitialisation::
             touchVertexFirstTimeSpecification() &
         exahype::adapters::
             PatchInitialisationAndExport2MultiscaleLinkedCell_1::
                 touchVertexFirstTimeSpecification() &
         exahype::mappings::VTKExport::touchVertexFirstTimeSpecification()

      ;
}

peano::MappingSpecification
exahype::adapters::PatchInitialisationAndExport::enterCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::PatchInitialisation::enterCellSpecification() &
         exahype::adapters::
             PatchInitialisationAndExport2MultiscaleLinkedCell_1::
                 enterCellSpecification() &
         exahype::mappings::VTKExport::enterCellSpecification()

      ;
}

peano::MappingSpecification
exahype::adapters::PatchInitialisationAndExport::leaveCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::PatchInitialisation::leaveCellSpecification() &
         exahype::adapters::
             PatchInitialisationAndExport2MultiscaleLinkedCell_1::
                 leaveCellSpecification() &
         exahype::mappings::VTKExport::leaveCellSpecification()

      ;
}

peano::MappingSpecification
exahype::adapters::PatchInitialisationAndExport::ascendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::PatchInitialisation::ascendSpecification() &
         exahype::adapters::
             PatchInitialisationAndExport2MultiscaleLinkedCell_1::
                 ascendSpecification() &
         exahype::mappings::VTKExport::ascendSpecification()

      ;
}

peano::MappingSpecification
exahype::adapters::PatchInitialisationAndExport::descendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification() &
         exahype::mappings::PatchInitialisation::descendSpecification() &
         exahype::adapters::
             PatchInitialisationAndExport2MultiscaleLinkedCell_1::
                 descendSpecification() &
         exahype::mappings::VTKExport::descendSpecification()

      ;
}

exahype::adapters::PatchInitialisationAndExport::
    PatchInitialisationAndExport() {}

exahype::adapters::PatchInitialisationAndExport::
    ~PatchInitialisationAndExport() {}

#if defined(SharedMemoryParallelisation)
exahype::adapters::PatchInitialisationAndExport::PatchInitialisationAndExport(
    const PatchInitialisationAndExport& masterThread)
    : _map2PatchInitialisation(masterThread._map2PatchInitialisation),
      _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1(
          masterThread
              ._map2PatchInitialisationAndExport2MultiscaleLinkedCell_1),
      _map2VTKExport(masterThread._map2VTKExport)

{}

void exahype::adapters::PatchInitialisationAndExport::mergeWithWorkerThread(
    const PatchInitialisationAndExport& workerThread) {
  _map2PatchInitialisation.mergeWithWorkerThread(
      workerThread._map2PatchInitialisation);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1
      .mergeWithWorkerThread(
          workerThread
              ._map2PatchInitialisationAndExport2MultiscaleLinkedCell_1);
  _map2VTKExport.mergeWithWorkerThread(workerThread._map2VTKExport);
}
#endif

void exahype::adapters::PatchInitialisationAndExport::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2PatchInitialisation.createHangingVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.createHangingVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2VTKExport.createHangingVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::PatchInitialisationAndExport::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2PatchInitialisation.destroyHangingVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.destroyHangingVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2VTKExport.destroyHangingVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::PatchInitialisationAndExport::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2PatchInitialisation.createInnerVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.createInnerVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2VTKExport.createInnerVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::PatchInitialisationAndExport::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2PatchInitialisation.createBoundaryVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.createBoundaryVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2VTKExport.createBoundaryVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::PatchInitialisationAndExport::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2PatchInitialisation.destroyVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.destroyVertex(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2VTKExport.destroyVertex(fineGridVertex, fineGridX, fineGridH,
                               coarseGridVertices, coarseGridVerticesEnumerator,
                               coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::PatchInitialisationAndExport::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2PatchInitialisation.createCell(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.createCell(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell);
  _map2VTKExport.createCell(fineGridCell, fineGridVertices,
                            fineGridVerticesEnumerator, coarseGridVertices,
                            coarseGridVerticesEnumerator, coarseGridCell,
                            fineGridPositionOfCell);
}

void exahype::adapters::PatchInitialisationAndExport::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2PatchInitialisation.destroyCell(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.destroyCell(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell);
  _map2VTKExport.destroyCell(fineGridCell, fineGridVertices,
                             fineGridVerticesEnumerator, coarseGridVertices,
                             coarseGridVerticesEnumerator, coarseGridCell,
                             fineGridPositionOfCell);
}

#ifdef Parallel
void exahype::adapters::PatchInitialisationAndExport::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  _map2PatchInitialisation.mergeWithNeighbour(vertex, neighbour, fromRank,
                                              fineGridX, fineGridH, level);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.mergeWithNeighbour(
      vertex, neighbour, fromRank, fineGridX, fineGridH, level);
  _map2VTKExport.mergeWithNeighbour(vertex, neighbour, fromRank, fineGridX,
                                    fineGridH, level);
}

void exahype::adapters::PatchInitialisationAndExport::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2PatchInitialisation.prepareSendToNeighbour(vertex, toRank, x, h, level);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1
      .prepareSendToNeighbour(vertex, toRank, x, h, level);
  _map2VTKExport.prepareSendToNeighbour(vertex, toRank, x, h, level);
}

void exahype::adapters::PatchInitialisationAndExport::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2PatchInitialisation.prepareCopyToRemoteNode(localVertex, toRank, x, h,
                                                   level);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1
      .prepareCopyToRemoteNode(localVertex, toRank, x, h, level);
  _map2VTKExport.prepareCopyToRemoteNode(localVertex, toRank, x, h, level);
}

void exahype::adapters::PatchInitialisationAndExport::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2PatchInitialisation.prepareCopyToRemoteNode(localCell, toRank, x, h,
                                                   level);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1
      .prepareCopyToRemoteNode(localCell, toRank, x, h, level);
  _map2VTKExport.prepareCopyToRemoteNode(localCell, toRank, x, h, level);
}

void exahype::adapters::PatchInitialisationAndExport::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2PatchInitialisation.mergeWithRemoteDataDueToForkOrJoin(
      localVertex, masterOrWorkerVertex, fromRank, x, h, level);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1
      .mergeWithRemoteDataDueToForkOrJoin(localVertex, masterOrWorkerVertex,
                                          fromRank, x, h, level);
  _map2VTKExport.mergeWithRemoteDataDueToForkOrJoin(
      localVertex, masterOrWorkerVertex, fromRank, x, h, level);
}

void exahype::adapters::PatchInitialisationAndExport::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2PatchInitialisation.mergeWithRemoteDataDueToForkOrJoin(
      localCell, masterOrWorkerCell, fromRank, x, h, level);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1
      .mergeWithRemoteDataDueToForkOrJoin(localCell, masterOrWorkerCell,
                                          fromRank, x, h, level);
  _map2VTKExport.mergeWithRemoteDataDueToForkOrJoin(
      localCell, masterOrWorkerCell, fromRank, x, h, level);
}

bool exahype::adapters::PatchInitialisationAndExport::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  bool result = false;
  result |= _map2PatchInitialisation.prepareSendToWorker(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell, worker);
  result |= _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1
                .prepareSendToWorker(
                    fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
                    coarseGridVertices, coarseGridVerticesEnumerator,
                    coarseGridCell, fineGridPositionOfCell, worker);
  result |= _map2VTKExport.prepareSendToWorker(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell, worker);

  return result;
}

void exahype::adapters::PatchInitialisationAndExport::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2PatchInitialisation.prepareSendToMaster(
      localCell, vertices, verticesEnumerator, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.prepareSendToMaster(
      localCell, vertices, verticesEnumerator, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell);
  _map2VTKExport.prepareSendToMaster(
      localCell, vertices, verticesEnumerator, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell);
}

void exahype::adapters::PatchInitialisationAndExport::mergeWithMaster(
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
  _map2PatchInitialisation.mergeWithMaster(
      workerGridCell, workerGridVertices, workerEnumerator, fineGridCell,
      fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell,
      worker, workerState, masterState);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.mergeWithMaster(
      workerGridCell, workerGridVertices, workerEnumerator, fineGridCell,
      fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell,
      worker, workerState, masterState);
  _map2VTKExport.mergeWithMaster(
      workerGridCell, workerGridVertices, workerEnumerator, fineGridCell,
      fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell,
      worker, workerState, masterState);
}

void exahype::adapters::PatchInitialisationAndExport::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2PatchInitialisation.receiveDataFromMaster(
      receivedCell, receivedVertices, receivedVerticesEnumerator,
      receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator,
      receivedCoarseGridCell, workersCoarseGridVertices,
      workersCoarseGridVerticesEnumerator, workersCoarseGridCell,
      fineGridPositionOfCell);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1
      .receiveDataFromMaster(
          receivedCell, receivedVertices, receivedVerticesEnumerator,
          receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator,
          receivedCoarseGridCell, workersCoarseGridVertices,
          workersCoarseGridVerticesEnumerator, workersCoarseGridCell,
          fineGridPositionOfCell);
  _map2VTKExport.receiveDataFromMaster(
      receivedCell, receivedVertices, receivedVerticesEnumerator,
      receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator,
      receivedCoarseGridCell, workersCoarseGridVertices,
      workersCoarseGridVerticesEnumerator, workersCoarseGridCell,
      fineGridPositionOfCell);
}

void exahype::adapters::PatchInitialisationAndExport::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  _map2PatchInitialisation.mergeWithWorker(localCell, receivedMasterCell,
                                           cellCentre, cellSize, level);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.mergeWithWorker(
      localCell, receivedMasterCell, cellCentre, cellSize, level);
  _map2VTKExport.mergeWithWorker(localCell, receivedMasterCell, cellCentre,
                                 cellSize, level);
}

void exahype::adapters::PatchInitialisationAndExport::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  _map2PatchInitialisation.mergeWithWorker(localVertex, receivedMasterVertex, x,
                                           h, level);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.mergeWithWorker(
      localVertex, receivedMasterVertex, x, h, level);
  _map2VTKExport.mergeWithWorker(localVertex, receivedMasterVertex, x, h,
                                 level);
}
#endif

void exahype::adapters::PatchInitialisationAndExport::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2PatchInitialisation.touchVertexFirstTime(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.touchVertexFirstTime(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2VTKExport.touchVertexFirstTime(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::PatchInitialisationAndExport::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  _map2PatchInitialisation.touchVertexLastTime(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.touchVertexLastTime(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
  _map2VTKExport.touchVertexLastTime(
      fineGridVertex, fineGridX, fineGridH, coarseGridVertices,
      coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex);
}

void exahype::adapters::PatchInitialisationAndExport::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2PatchInitialisation.enterCell(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.enterCell(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell);
  _map2VTKExport.enterCell(fineGridCell, fineGridVertices,
                           fineGridVerticesEnumerator, coarseGridVertices,
                           coarseGridVerticesEnumerator, coarseGridCell,
                           fineGridPositionOfCell);
}

void exahype::adapters::PatchInitialisationAndExport::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  _map2PatchInitialisation.leaveCell(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.leaveCell(
      fineGridCell, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell,
      fineGridPositionOfCell);
  _map2VTKExport.leaveCell(fineGridCell, fineGridVertices,
                           fineGridVerticesEnumerator, coarseGridVertices,
                           coarseGridVerticesEnumerator, coarseGridCell,
                           fineGridPositionOfCell);
}

void exahype::adapters::PatchInitialisationAndExport::beginIteration(
    exahype::State& solverState) {
  _map2PatchInitialisation.beginIteration(solverState);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.beginIteration(
      solverState);
  _map2VTKExport.beginIteration(solverState);
}

void exahype::adapters::PatchInitialisationAndExport::endIteration(
    exahype::State& solverState) {
  _map2PatchInitialisation.endIteration(solverState);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.endIteration(
      solverState);
  _map2VTKExport.endIteration(solverState);
}

void exahype::adapters::PatchInitialisationAndExport::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  _map2PatchInitialisation.descend(
      fineGridCells, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.descend(
      fineGridCells, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell);
  _map2VTKExport.descend(fineGridCells, fineGridVertices,
                         fineGridVerticesEnumerator, coarseGridVertices,
                         coarseGridVerticesEnumerator, coarseGridCell);
}

void exahype::adapters::PatchInitialisationAndExport::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  _map2PatchInitialisation.ascend(
      fineGridCells, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell);
  _map2PatchInitialisationAndExport2MultiscaleLinkedCell_1.ascend(
      fineGridCells, fineGridVertices, fineGridVerticesEnumerator,
      coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell);
  _map2VTKExport.ascend(fineGridCells, fineGridVertices,
                        fineGridVerticesEnumerator, coarseGridVertices,
                        coarseGridVerticesEnumerator, coarseGridCell);
}
