#include "exahype/mappings/AugmentedAMRTreePlot2d.h"

#include <sstream>

#include "peano/grid/CellFlags.h"
#include "peano/utils/Loop.h"

#include "exahype/solvers/Solver.h"

#ifdef Parallel
#include "tarch/parallel/Node.h"
#endif

tarch::logging::Log exahype::mappings::AugmentedAMRTreePlot2d::_log(
    "exahype::mappings::AugmentedAMRTreePlot2d");

int exahype::mappings::AugmentedAMRTreePlot2d::_snapshotCounter =
    0;
double exahype::mappings::AugmentedAMRTreePlot2d::SqueezeZAxis =
    4.0;
double exahype::mappings::AugmentedAMRTreePlot2d::
    TreeConnectionsValue = -100.0;

peano::CommunicationSpecification exahype::mappings::
    AugmentedAMRTreePlot2d::communicationSpecification() {
  return peano::CommunicationSpecification::getPessimisticSpecification();
}

peano::MappingSpecification exahype::mappings::
    AugmentedAMRTreePlot2d::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

peano::MappingSpecification exahype::mappings::
    AugmentedAMRTreePlot2d::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,
                                     peano::MappingSpecification::Serial);
}

peano::MappingSpecification exahype::mappings::
    AugmentedAMRTreePlot2d::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,
                                     peano::MappingSpecification::Serial);
}

peano::MappingSpecification exahype::mappings::
    AugmentedAMRTreePlot2d::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}

peano::MappingSpecification
exahype::mappings::AugmentedAMRTreePlot2d::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

peano::MappingSpecification exahype::mappings::
    AugmentedAMRTreePlot2d::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

std::map<tarch::la::Vector<DIMENSIONS + 1, double>, int,
         tarch::la::VectorCompare<DIMENSIONS + 1> >
    exahype::mappings::AugmentedAMRTreePlot2d::_vertex2IndexMap;
std::map<tarch::la::Vector<DIMENSIONS + 1, double>, int,
         tarch::la::VectorCompare<DIMENSIONS + 1> >
    exahype::mappings::AugmentedAMRTreePlot2d::
        _cellCenter2IndexMap;

exahype::mappings::AugmentedAMRTreePlot2d::
    AugmentedAMRTreePlot2d()
    : _vtkWriter(0), _vertexWriter(0), _cellWriter(0), _cellNumberWriter(0),
      _cellTypeWriter(0), _cellDescriptionIndexWriter(0), _cellRefinementEventWriter(0),
      _cellDataWriter(0),_cellCounter(0){}

exahype::mappings::AugmentedAMRTreePlot2d::
    ~AugmentedAMRTreePlot2d() {}

#if defined(SharedMemoryParallelisation)
exahype::mappings::AugmentedAMRTreePlot2d::
    AugmentedAMRTreePlot2d(
        const AugmentedAMRTreePlot2d& masterThread)
    : _vtkWriter(masterThread._vtkWriter),
      _vertexWriter(masterThread._vertexWriter),
      _cellWriter(masterThread._cellWriter),
      _cellNumberWriter(masterThread._cellNumberWriter),
      _cellTypeWriter(masterThread._cellTypeWriter),
      _cellDescriptionIndexWriter(masterThread._cellDescriptionIndexWriter),
      _cellRefinementEventWriter(masterThread._cellRefinementEventWriter),
      _cellDataWriter(masterThread._cellDataWriter) {}

void exahype::mappings::AugmentedAMRTreePlot2d::
    mergeWithWorkerThread(
        const AugmentedAMRTreePlot2d& workerThread) {}
#endif

void exahype::mappings::AugmentedAMRTreePlot2d::plotVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x, int level) {
  tarch::la::Vector<DIMENSIONS + 1, double> y;
  for (int i = 0; i < DIMENSIONS; i++) y(i) = x(i);
  y(DIMENSIONS) = level;

  tarch::la::Vector<DIMENSIONS + 1, double> plotY = y;
  plotY(DIMENSIONS) = level / SqueezeZAxis;

#if DIMENSIONS==2
  if (_vertex2IndexMap.find(y) == _vertex2IndexMap.end()) {
    _vertex2IndexMap[y] = _vertexWriter->plotVertex(plotY);
  }
#else
  logError("plotVertex","This mapping can only be used for two-dimensional problems (DIMENSIONS==2).")
#endif
}

void exahype::mappings::AugmentedAMRTreePlot2d::
    createHangingVertex(
        exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  int minLevel = std::numeric_limits<int>::max();
  for (std::vector<exahype::solvers::Solver*>::const_iterator p =
      exahype::solvers::RegisteredSolvers.begin();
      p != exahype::solvers::RegisteredSolvers.end();
      ++p) {  // @todo replace by parloops?
    minLevel = std::min(minLevel,(*p)->getMinimumTreeDepth()+1);
  }

#if DIMENSIONS==2
  if (coarseGridVerticesEnumerator.getLevel() + 1 >= minLevel) {
    plotVertex(fineGridVertex, fineGridX,
        coarseGridVerticesEnumerator.getLevel() + 1 - minLevel + 1);
  }
#else
  logError("createHangingVertex","This mapping can only be used for two-dimensional problems (DIMENSIONS==2).")
#endif
}

void exahype::mappings::AugmentedAMRTreePlot2d::
    destroyHangingVertex(
        const exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::AugmentedAMRTreePlot2d::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::AugmentedAMRTreePlot2d::
    createBoundaryVertex(
        exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::AugmentedAMRTreePlot2d::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::AugmentedAMRTreePlot2d::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::AugmentedAMRTreePlot2d::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

#ifdef Parallel
void exahype::mappings::AugmentedAMRTreePlot2d::
    mergeWithNeighbour(exahype::Vertex& vertex,
                       const exahype::Vertex& neighbour, int fromRank,
                       const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
                       const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
                       int level) {}

void exahype::mappings::AugmentedAMRTreePlot2d::
    prepareSendToNeighbour(exahype::Vertex& vertex, int toRank,
                           const tarch::la::Vector<DIMENSIONS, double>& x,
                           const tarch::la::Vector<DIMENSIONS, double>& h,
                           int level) {}

void exahype::mappings::AugmentedAMRTreePlot2d::
    prepareCopyToRemoteNode(exahype::Vertex& localVertex, int toRank,
                            const tarch::la::Vector<DIMENSIONS, double>& x,
                            const tarch::la::Vector<DIMENSIONS, double>& h,
                            int level) {}

void exahype::mappings::AugmentedAMRTreePlot2d::
    prepareCopyToRemoteNode(
        exahype::Cell& localCell, int toRank,
        const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {}

void exahype::mappings::AugmentedAMRTreePlot2d::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}

void exahype::mappings::AugmentedAMRTreePlot2d::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}

bool exahype::mappings::AugmentedAMRTreePlot2d::
    prepareSendToWorker(
        exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
        int worker) {
  return false;
}

void exahype::mappings::AugmentedAMRTreePlot2d::
    prepareSendToMaster(
        exahype::Cell& localCell, exahype::Vertex* vertices,
        const peano::grid::VertexEnumerator& verticesEnumerator,
        const exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        const exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::AugmentedAMRTreePlot2d::mergeWithMaster(
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
    exahype::State& masterState) {}

void exahype::mappings::AugmentedAMRTreePlot2d::
    receiveDataFromMaster(
        exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
        const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
        exahype::Vertex* const receivedCoarseGridVertices,
        const peano::grid::VertexEnumerator&
            receivedCoarseGridVerticesEnumerator,
        exahype::Cell& receivedCoarseGridCell,
        exahype::Vertex* const workersCoarseGridVertices,
        const peano::grid::VertexEnumerator&
            workersCoarseGridVerticesEnumerator,
        exahype::Cell& workersCoarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::AugmentedAMRTreePlot2d::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {}

void exahype::mappings::AugmentedAMRTreePlot2d::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}
#endif

void exahype::mappings::AugmentedAMRTreePlot2d::
    touchVertexFirstTime(
        exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  int minLevel = std::numeric_limits<int>::max();
  for (std::vector<exahype::solvers::Solver*>::const_iterator p =
      exahype::solvers::RegisteredSolvers.begin();
      p != exahype::solvers::RegisteredSolvers.end();
      ++p) {  // @todo replace by parloops?
    minLevel = std::min(minLevel,(*p)->getMinimumTreeDepth()+1);
  }

#if DIMENSIONS==2
  if (coarseGridVerticesEnumerator.getLevel() + 1 >= minLevel) {
    plotVertex(fineGridVertex, fineGridX,
               coarseGridVerticesEnumerator.getLevel() + 1 - minLevel + 1);
  }
#else
  logError("touchVertexFirstTime","This mapping can only be used for two-dimensional problems (DIMENSIONS==2).")
#endif
}

void exahype::mappings::AugmentedAMRTreePlot2d::
    touchVertexLastTime(
        exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::AugmentedAMRTreePlot2d::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  int minLevel = std::numeric_limits<int>::max();
  for (std::vector<exahype::solvers::Solver*>::const_iterator p =
      exahype::solvers::RegisteredSolvers.begin();
      p != exahype::solvers::RegisteredSolvers.end();
      ++p) {  // @todo replace by parloops?
    minLevel = std::min(minLevel,(*p)->getMinimumTreeDepth()+1);
  }

  if (coarseGridVerticesEnumerator.getLevel() + 1 >= minLevel) {
    int vertexIndex[TWO_POWER_D];
    tarch::la::Vector<DIMENSIONS + 1, double> currentVertexPosition;
    currentVertexPosition(DIMENSIONS) = fineGridVerticesEnumerator.getLevel() - minLevel + 1;

    dfor2(i)
    for (int d = 0; d < DIMENSIONS; d++) {
      currentVertexPosition(d) =
          fineGridVerticesEnumerator.getVertexPosition(i)(d);
    }
    assertion2(
        _vertex2IndexMap.find(currentVertexPosition) != _vertex2IndexMap.end(),
        currentVertexPosition,
        fineGridVertices[fineGridVerticesEnumerator(i)].toString());
    vertexIndex[iScalar] = _vertex2IndexMap[currentVertexPosition];
    enddforx

    _cellWriter->plotQuadrangle(vertexIndex);

    const int cellIndex = _cellWriter->plotQuadrangle(vertexIndex);

    _cellNumberWriter->plotCell(cellIndex, _cellCounter);

    _cellDescriptionIndexWriter->plotCell(cellIndex, static_cast<int>(fineGridCell.getADERDGCellDescriptionsIndex()));

    if (ADERDGCellDescriptionHeap::getInstance().isValidIndex(
        fineGridCell.getADERDGCellDescriptionsIndex()) &&
        ADERDGCellDescriptionHeap::getInstance()
    .getData(fineGridCell.getADERDGCellDescriptionsIndex())
    .size() > 0) {
      int solverNumber = 0;
      bool solverFound  = false;

      for (std::vector<exahype::records::ADERDGCellDescription>::iterator pFine =
          ADERDGCellDescriptionHeap::getInstance()
      .getData(fineGridCell.getADERDGCellDescriptionsIndex())
      .begin();
          pFine != ADERDGCellDescriptionHeap::getInstance()
          .getData(fineGridCell.getADERDGCellDescriptionsIndex())
          .end();
          ++pFine) {
        if (pFine->getSolverNumber()==solverNumber) {
          _cellTypeWriter->plotCell(cellIndex, static_cast<int>(pFine->getType()));
          _cellRefinementEventWriter->plotCell(cellIndex, static_cast<int>(pFine->getRefinementEvent()));
          _cellDataWriter->plotCell(cellIndex,2*static_cast<int>(pFine->getSolution() > -1) +
              static_cast<int>(pFine->getExtrapolatedPredictor() > -1));
          solverFound = true;
        }
      }

      if (!solverFound) {
        _cellTypeWriter->plotCell(cellIndex, -1);
        _cellRefinementEventWriter->plotCell(cellIndex, -1);
        _cellDataWriter->plotCell(cellIndex,0);
      }

    } else {
      _cellTypeWriter->plotCell(cellIndex, static_cast<int>(fineGridCell.getADERDGCellDescriptionsIndex()));
      _cellRefinementEventWriter->plotCell(cellIndex, -1);
      _cellDataWriter->plotCell(cellIndex,0);
    }

    _cellCounter++;
  }
}

void exahype::mappings::AugmentedAMRTreePlot2d::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::AugmentedAMRTreePlot2d::beginIteration(
    exahype::State& solverState) {
  assertion(_vtkWriter == 0);

  _cellCounter = 0;

  _vtkWriter = new UsedWriter();

  _vertexWriter = _vtkWriter->createVertexWriter();
  _cellWriter = _vtkWriter->createCellWriter();

  _cellNumberWriter = _vtkWriter->createCellDataWriter("cell-number", 1);
  _cellTypeWriter  = _vtkWriter->createCellDataWriter("cell-type(NoPatch=-1,Erased=0,Ancestor=1,EmptyAncestor=2,Cell=3,Descendant=4,EmptyDescendant=5)", 1);
  _cellDescriptionIndexWriter = _vtkWriter->createCellDataWriter("NoPatch=-1,ValidPatch>=0", 1);
  _cellRefinementEventWriter = _vtkWriter->createCellDataWriter("None=0,ErasingRequested=1,Restricting=2,Erasing=3,AllocatingMemory=4,ErasingChildren=5,RefiningRequested=6,Refining=7,Prolongating=8,DeaugmentingRequested=9,AugmentingRequested=10,Augmenting=11)", 1);
  _cellDataWriter = _vtkWriter->createCellDataWriter("None=0,OnlyFaceData=1,VolumeAndFaceData=3)", 1);
}

void exahype::mappings::AugmentedAMRTreePlot2d::endIteration(
    exahype::State& solverState) {
  _vertexWriter->close();
  _cellWriter->close();

  _cellTypeWriter->close();
  _cellDescriptionIndexWriter->close();
  _cellRefinementEventWriter->close();
  _cellDataWriter->close();
  _cellNumberWriter->close();

  delete _vertexWriter;
  delete _cellWriter;

  delete _cellTypeWriter;
  delete _cellDescriptionIndexWriter;
  delete _cellRefinementEventWriter;
  delete _cellNumberWriter;
  delete _cellDataWriter;

  _vertexWriter = 0;
  _cellWriter = 0;

  _cellNumberWriter = 0;
  _cellTypeWriter = 0;
  _cellDescriptionIndexWriter = 0;
  _cellRefinementEventWriter = 0;
  _cellDataWriter = 0;

  std::ostringstream snapshotFileName;
  snapshotFileName << "tree"
#ifdef Parallel
                   << "-rank-" << tarch::parallel::Node::getInstance().getRank()
#endif
                   << "-" << _snapshotCounter << ".vtk";
  _vtkWriter->writeToFile(snapshotFileName.str());

  _snapshotCounter++;

  _vertex2IndexMap.clear();

  delete _vtkWriter;
  _vtkWriter = 0;
}

void exahype::mappings::AugmentedAMRTreePlot2d::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {}

void exahype::mappings::AugmentedAMRTreePlot2d::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {}
