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
 *
 * \author Dominic E. Charrier, Tobias Weinzierl, Jean-Matthieu Gallard, Fabian Güra
 **/
#include "exahype/solvers/ADERDGSolver.h"


#include <limits>

#include "exahype/Cell.h"
#include "exahype/Vertex.h"
#include "exahype/VertexOperations.h"

#include "tarch/la/VectorVectorOperations.h"
#include "tarch/multicore/Lock.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"

#include "peano/heap/CompressedFloatingPointNumbers.h"
#include "peano/datatraversal/TaskSet.h"

#include "exahype/solvers/LimitingADERDGSolver.h"


namespace {
  constexpr const char* tags[]{"solutionUpdate",
                             "volumeIntegral",
                             "surfaceIntegral",
                             "riemannSolver",
                             "spaceTimePredictor",
                             "stableTimeStepSize",
                             "solutionAdjustment",
                             "faceUnknownsProlongation",
                             "faceUnknownsRestriction",
                             "volumeUnknownsProlongation",
                             "volumeUnknownsRestriction",
                             "boundaryConditions",
                             "pointSource"
                             };
  constexpr const char* deepProfilingTags[]{
                             "spaceTimePredictor_PDEflux",
                             "spaceTimePredictor_PDEsource",
                             "spaceTimePredictor_overhead",
                             "spaceTimePredictor_PDEncp",
                             "riemannSolver_PDEmatrixb",
                             "riemannSolver_overhead"
  };
  typedef peano::heap::PlainCharHeap CompressedDataHeap;
}


#ifdef Asserts
/**
 * If you enable assertions, we have the option not to remove any entries from
 * any heap but to continue to store all unknown on the standard heap when we
 * compress. This allows us to validate that the data that is compressed in one
 * iteration and uncompressed in the next one does not differ too significantly
 * from the original data. There are however two drawbacks to this approach:
 *
 * - It is costly.
 * - It changes the code semantics - we actually work with other and more heap
 *   entries and thus cannot claim that a code with these assertions equals a
 *   code without any assertions.
 *
 * I thus decided to trigger the comparison of compressed vs. uncompressed data
 * through a special flag.
 */
//#define ValidateCompressedVsUncompressedData

double exahype::solvers::ADERDGSolver::PipedUncompressedBytes = 0;
double exahype::solvers::ADERDGSolver::PipedCompressedBytes = 0;
#endif



tarch::logging::Log exahype::solvers::ADERDGSolver::_log( "exahype::solvers::ADERDGSolver");


double exahype::solvers::ADERDGSolver::CompressionAccuracy = 0.0;

bool exahype::solvers::ADERDGSolver::SpawnCompressionAsBackgroundThread = false;

void exahype::solvers::ADERDGSolver::addNewCellDescription(
  const int cellDescriptionsIndex,
  const int                                      solverNumber,
  const exahype::records::ADERDGCellDescription::Type cellType,
  const exahype::records::ADERDGCellDescription::RefinementEvent refinementEvent,
  const int                                     level,
  const int                                     parentIndex,
  const tarch::la::Vector<DIMENSIONS, double>&  cellSize,
  const tarch::la::Vector<DIMENSIONS, double>&  cellOffset) {

  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion2(parentIndex == -1 || parentIndex != cellDescriptionsIndex, parentIndex, cellDescriptionsIndex);
  assertion2(parentIndex != cellDescriptionsIndex, parentIndex, cellDescriptionsIndex);

  assertion2(static_cast<unsigned int>(solverNumber) < solvers::RegisteredSolvers.size(),
             solverNumber, exahype::solvers::RegisteredSolvers.size());

  CellDescription newCellDescription;
  newCellDescription.setSolverNumber(solverNumber);

  // Default AMR settings
  newCellDescription.setType(cellType);
  newCellDescription.setParentIndex(parentIndex);
  newCellDescription.setLevel(level);
  newCellDescription.setRefinementEvent(refinementEvent);

  std::bitset<DIMENSIONS_TIMES_TWO>
      riemannSolvePerformed;  // default construction: no bit set
  newCellDescription.setRiemannSolvePerformed(riemannSolvePerformed);

  // Pass geometry information to the cellDescription description
  newCellDescription.setSize(cellSize);
  newCellDescription.setOffset(cellOffset);

  // Initialise MPI helper variables
  #ifdef Parallel
  newCellDescription.setHasToHoldDataForNeighbourCommunication(false);
  newCellDescription.setHasToHoldDataForMasterWorkerCommunication(false);
  for (int faceIndex = 0; faceIndex < DIMENSIONS_TIMES_TWO; faceIndex++) {
    newCellDescription.setFaceDataExchangeCounter(faceIndex,TWO_POWER_D);
  }
  #endif

  // Default field data indices
  newCellDescription.setSolution(-1);
  newCellDescription.setUpdate(-1);
  newCellDescription.setExtrapolatedPredictor(-1);
  newCellDescription.setFluctuation(-1);

  // Limiter meta data (oscillations identificator)
  newCellDescription.setLimiterStatus(0); // 0 is CellDescription::LimiterStatus::Ok
  newCellDescription.setPreviousLimiterStatus(0);
  newCellDescription.setFacewiseLimiterStatus(0);
  newCellDescription.setSolutionMin(-1);
  newCellDescription.setSolutionMax(-1);

  // Compression
  newCellDescription.setCompressionState(exahype::records::ADERDGCellDescription::CompressionState::Uncompressed);
  newCellDescription.setSolutionAverages(-1);
  newCellDescription.setUpdateAverages(-1);
  newCellDescription.setExtrapolatedPredictorAverages(-1);
  newCellDescription.setFluctuationAverages(-1);

  newCellDescription.setSolutionCompressed(-1);
  newCellDescription.setUpdateCompressed(-1);
  newCellDescription.setExtrapolatedPredictorCompressed(-1);
  newCellDescription.setFluctuationCompressed(-1);

  ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex).push_back(newCellDescription);
}

/**
 * Returns the ADERDGCellDescription heap vector
 * at address \p cellDescriptionsIndex.
 */
exahype::solvers::ADERDGSolver::Heap::HeapEntries& exahype::solvers::ADERDGSolver::getCellDescriptions(
    const int cellDescriptionsIndex) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);

  return Heap::getInstance().getData(cellDescriptionsIndex);
}

/**
 * Returns the ADERDGCellDescription with index \p element
 * in the heap vector at address \p cellDescriptionsIndex.
 */
exahype::solvers::ADERDGSolver::CellDescription& exahype::solvers::ADERDGSolver::getCellDescription(
    const int cellDescriptionsIndex,
    const int element) {
  assertion2(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex,element);
  assertion2(element>=0,cellDescriptionsIndex,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),cellDescriptionsIndex,element);

  return Heap::getInstance().getData(cellDescriptionsIndex)[element];
}

/**
 * Returns if a ADERDGCellDescription type holds face data.
 */
bool exahype::solvers::ADERDGSolver::holdsFaceData(const exahype::solvers::ADERDGSolver::CellDescription::Type& cellDescriptionType) {
  return cellDescriptionType==CellDescription::Cell       ||
      cellDescriptionType==CellDescription::Ancestor   ||
      cellDescriptionType==CellDescription::Descendant;
}

void exahype::solvers::ADERDGSolver::ensureNoUnnecessaryMemoryIsAllocated(exahype::records::ADERDGCellDescription& cellDescription) {

  if (DataHeap::getInstance().isValidIndex(cellDescription.getSolution())) {
    switch (cellDescription.getType()) {
      case exahype::records::ADERDGCellDescription::Erased:
      case exahype::records::ADERDGCellDescription::EmptyAncestor:
      case exahype::records::ADERDGCellDescription::Ancestor:
      case exahype::records::ADERDGCellDescription::EmptyDescendant:
      case exahype::records::ADERDGCellDescription::Descendant:
        {
        waitUntilAllBackgroundTasksHaveTerminated();

        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()));
        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolution()));
        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getUpdate()));

        tarch::multicore::Lock lock(_heapSemaphore);

        if (cellDescription.getUpdate()>=0) {
          DataHeap::getInstance().deleteData(cellDescription.getUpdate());
          assertion(cellDescription.getUpdateCompressed()==-1);
        }
        else {
          assertion(CompressionAccuracy>0.0);
          assertion(cellDescription.getUpdate()==-1);
          CompressedDataHeap::getInstance().deleteData(cellDescription.getUpdateCompressed());
        }

        if (cellDescription.getSolution()>=0) {
          DataHeap::getInstance().deleteData(cellDescription.getSolution());
          assertion(cellDescription.getSolutionCompressed()==-1);
        }
        else {
          assertion(CompressionAccuracy>0.0);
          assertion(cellDescription.getSolution()==-1);
          CompressedDataHeap::getInstance().deleteData(cellDescription.getSolutionCompressed());
        }

        DataHeap::getInstance().deleteData(cellDescription.getUpdateAverages());
        DataHeap::getInstance().deleteData(cellDescription.getSolutionAverages());
        DataHeap::getInstance().deleteData(cellDescription.getPreviousSolutionAverages());

        cellDescription.setPreviousSolution(-1);
        cellDescription.setSolution(-1);
        cellDescription.setUpdate(-1);

        cellDescription.setPreviousSolutionAverages(-1);
        cellDescription.setSolutionAverages(-1);
        cellDescription.setUpdateAverages(-1);

        cellDescription.setPreviousSolutionCompressed(-1);
        cellDescription.setSolutionCompressed(-1);
        cellDescription.setUpdateCompressed(-1);
        }
        break;
      default:
        break;
    }
  }

  if (DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor())) {
    switch (cellDescription.getType()) {
      case exahype::records::ADERDGCellDescription::Erased:
      case exahype::records::ADERDGCellDescription::EmptyAncestor:
      case exahype::records::ADERDGCellDescription::EmptyDescendant:
        {
        waitUntilAllBackgroundTasksHaveTerminated();

        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()));
        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionMin()));
        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionMax()));

        tarch::multicore::Lock lock(_heapSemaphore);

        if (cellDescription.getExtrapolatedPredictor()>=0) {
          DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictor());
          assertion(cellDescription.getExtrapolatedPredictorCompressed()==-1);
        }
        else {
          assertion(CompressionAccuracy>0.0);
          assertion(cellDescription.getExtrapolatedPredictor()==-1);
          CompressedDataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictorCompressed());
        }

        if (cellDescription.getFluctuation()>=0) {
          DataHeap::getInstance().deleteData(cellDescription.getFluctuation());
          assertion(cellDescription.getFluctuationCompressed()==-1);
        }
        else {
          assertion(CompressionAccuracy>0.0);
          assertion(cellDescription.getFluctuation()==-1);
          CompressedDataHeap::getInstance().deleteData(cellDescription.getFluctuationCompressed());
        }

        DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictorAverages());
        DataHeap::getInstance().deleteData(cellDescription.getFluctuationAverages());

        DataHeap::getInstance().deleteData(cellDescription.getSolutionMin());
        DataHeap::getInstance().deleteData(cellDescription.getSolutionMax());

        cellDescription.setExtrapolatedPredictor(-1);
        cellDescription.setFluctuation(-1);

        cellDescription.setExtrapolatedPredictorCompressed(-1);
        cellDescription.setFluctuationCompressed(-1);
        }
        break;
      default:
        break;
    }
  }
}

void exahype::solvers::ADERDGSolver::ensureNecessaryMemoryIsAllocated(exahype::records::ADERDGCellDescription& cellDescription) {
  switch (cellDescription.getType()) {
    case CellDescription::Cell:
      if (!DataHeap::getInstance().isValidIndex(cellDescription.getSolution())) {
        waitUntilAllBackgroundTasksHaveTerminated();

        tarch::multicore::Lock lock(_heapSemaphore);
        assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getUpdate()));
        // Allocate volume DoF for limiter
        const int dofPerCell        = getUnknownsPerCell();
        const int dataPointsPerCell = getDataPerCell(); // Only the solution and previousSolution store material parameters
        cellDescription.setPreviousSolution(DataHeap::getInstance().createData(dataPointsPerCell, dataPointsPerCell, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired));
        cellDescription.setSolution(DataHeap::getInstance().createData(dataPointsPerCell, dataPointsPerCell, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired));
        cellDescription.setUpdate(DataHeap::getInstance().createData(dofPerCell, dofPerCell, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired));

        assertionEquals(DataHeap::getInstance().getData(cellDescription.getPreviousSolution()).size(),static_cast<unsigned int>(dataPointsPerCell));
        assertionEquals(DataHeap::getInstance().getData(cellDescription.getUpdate()).capacity(),static_cast<unsigned int>(dofPerCell));
        assertionEquals(DataHeap::getInstance().getData(cellDescription.getUpdate()).size(),static_cast<unsigned int>(dofPerCell));
        assertionEquals(DataHeap::getInstance().getData(cellDescription.getSolution()).capacity(),static_cast<unsigned int>(dataPointsPerCell));

        cellDescription.setUpdateCompressed(-1);
        cellDescription.setSolutionCompressed(-1);
        cellDescription.setPreviousSolutionCompressed(-1);

        if (CompressionAccuracy>0.0) {
          CompressedDataHeap::getInstance().reserveHeapEntriesForRecycling(2);
        }

        cellDescription.setPreviousSolutionAverages( DataHeap::getInstance().createData( getNumberOfVariables()+getNumberOfParameters(), getNumberOfVariables()+getNumberOfParameters(), DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired ) );
        cellDescription.setUpdateAverages(           DataHeap::getInstance().createData( getNumberOfVariables(), getNumberOfVariables(), DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired ) );
        cellDescription.setSolutionAverages(         DataHeap::getInstance().createData( getNumberOfVariables()+getNumberOfParameters(), getNumberOfVariables()+getNumberOfParameters(), DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired ) );

        assertionEquals3(
            DataHeap::getInstance().getData(cellDescription.getPreviousSolutionAverages()).size(),static_cast<unsigned int>(getNumberOfVariables() + getNumberOfParameters()),
            DataHeap::getInstance().getData(cellDescription.getPreviousSolutionAverages()).size(),static_cast<unsigned int>(getNumberOfVariables() + getNumberOfParameters()),
            getNumberOfVariables()
        );
        assertionEquals3(
          DataHeap::getInstance().getData(cellDescription.getUpdateAverages()).size(),static_cast<unsigned int>(getNumberOfVariables()),
          DataHeap::getInstance().getData(cellDescription.getUpdateAverages()).size(),static_cast<unsigned int>(getNumberOfVariables()),
          getNumberOfVariables()
        );
        assertionEquals3(
          DataHeap::getInstance().getData(cellDescription.getSolutionAverages()).size(),static_cast<unsigned int>(getNumberOfVariables() + getNumberOfParameters()),
          DataHeap::getInstance().getData(cellDescription.getSolutionAverages()).size(),static_cast<unsigned int>(getNumberOfVariables() + getNumberOfParameters()),
          getNumberOfVariables()
        );

        cellDescription.setCompressionState(exahype::records::ADERDGCellDescription::Uncompressed);
      }
      break;
    default:
      break;
  }

  switch (cellDescription.getType()) {
    case CellDescription::Cell:
    case CellDescription::Ancestor:
    case CellDescription::Descendant:
      if (!DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor())) {
        waitUntilAllBackgroundTasksHaveTerminated();

        tarch::multicore::Lock lock(_heapSemaphore);
        assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()));

        // Allocate face DoF
        const int dataPerBnd = getBndTotalSize();
        const int dofPerBnd  = getBndFluxTotalSize();

        cellDescription.setExtrapolatedPredictor(DataHeap::getInstance().createData(dataPerBnd, dataPerBnd, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired));
        cellDescription.setFluctuation(          DataHeap::getInstance().createData(dofPerBnd, dofPerBnd, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired));

        assertionEquals3(
          DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).size(),static_cast<unsigned int>(dataPerBnd),
          cellDescription.getExtrapolatedPredictor(),
          cellDescription.toString(),
          toString()
        );
        assertionEquals3(
          DataHeap::getInstance().getData(cellDescription.getFluctuation()).size(),static_cast<unsigned int>(dofPerBnd),
          cellDescription.getExtrapolatedPredictor(),
          cellDescription.toString(),
          toString()
        );

        cellDescription.setExtrapolatedPredictorCompressed(-1);
        cellDescription.setFluctuationCompressed(-1);

        if (CompressionAccuracy>0.0) {
          CompressedDataHeap::getInstance().reserveHeapEntriesForRecycling(2);
        }

        //TODO JMG / Dominic adapt for padding with optimized kernels
        int faceAverageCardinality = getNumberOfVariables() * 2 * DIMENSIONS;
        cellDescription.setExtrapolatedPredictorAverages( DataHeap::getInstance().createData( faceAverageCardinality, faceAverageCardinality, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired ) );
        cellDescription.setFluctuationAverages(           DataHeap::getInstance().createData( faceAverageCardinality, faceAverageCardinality, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired ) );

        // Allocate volume DoF for limiter (we need for every of the 2*DIMENSIONS faces an array of min values
        // and array of max values of the neighbour at this face).
        const int numberOfVariables = getNumberOfVariables();
        cellDescription.setSolutionMin(DataHeap::getInstance().createData(
            numberOfVariables * 2 * DIMENSIONS, numberOfVariables * 2 * DIMENSIONS, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired));
        cellDescription.setSolutionMax(DataHeap::getInstance().createData(
            numberOfVariables * 2 * DIMENSIONS, numberOfVariables * 2 * DIMENSIONS, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired));

        for (int i=0; i<numberOfVariables * 2 * DIMENSIONS; i++) {
          DataHeap::getInstance().getData( cellDescription.getSolutionMax() )[i] =  std::numeric_limits<double>::max();
          DataHeap::getInstance().getData( cellDescription.getSolutionMin() )[i] = -std::numeric_limits<double>::max();
        }
      }
      break;
    default:
      break;
  }
}

exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus
exahype::solvers::ADERDGSolver::toLimiterStatusEnum(const int limiterStatusAsInt) {
  assertion1( limiterStatusAsInt >= 0, limiterStatusAsInt );
  assertion1( limiterStatusAsInt <=
      static_cast<int>(CellDescription::LimiterStatus::Troubled),
      limiterStatusAsInt);

  return static_cast<CellDescription::LimiterStatus>(limiterStatusAsInt);
}

/**
 * Iterate over the merged limiter statuses per face and
 * determine a unique value.
 */
int
exahype::solvers::ADERDGSolver::determineLimiterStatus(
    exahype::solvers::ADERDGSolver::CellDescription& cellDescription) {
  // Assumes increasing value of limiter status the closer we get to troubled cell
  int limiterStatusAsInt = 0;
  for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
    limiterStatusAsInt = std::max( limiterStatusAsInt, cellDescription.getFacewiseLimiterStatus(i) );
  }

  return limiterStatusAsInt;
}

void exahype::solvers::ADERDGSolver::writeLimiterStatusOnBoundary(
    exahype::solvers::ADERDGSolver::CellDescription& cellDescription) {
  for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
    cellDescription.setFacewiseLimiterStatus(faceIndex,0); // 0 is CellDescription::LimiterStatus::0k
    if (cellDescription.getIsInside(faceIndex)) {
      cellDescription.setFacewiseLimiterStatus(faceIndex,cellDescription.getLimiterStatus());
    }
  }
}

exahype::solvers::ADERDGSolver::ADERDGSolver(
    const std::string& identifier, int numberOfVariables,
    int numberOfParameters, int nodesPerCoordinateAxis,
    double maximumMeshSize, int maximumAdaptiveMeshDepth,
    exahype::solvers::Solver::TimeStepping timeStepping,
    std::unique_ptr<profilers::Profiler> profiler)
    : Solver(identifier, Solver::Type::ADERDG, numberOfVariables,
             numberOfParameters, nodesPerCoordinateAxis,
             maximumMeshSize, maximumAdaptiveMeshDepth,
             timeStepping, std::move(profiler)),
     _previousMinCorrectorTimeStamp( std::numeric_limits<double>::max() ),
     _previousMinCorrectorTimeStepSize( std::numeric_limits<double>::max() ),
     _minCorrectorTimeStamp( std::numeric_limits<double>::max() ),
     _minCorrectorTimeStepSize( std::numeric_limits<double>::max() ),
     _minPredictorTimeStamp( std::numeric_limits<double>::max() ),
     _minPredictorTimeStepSize( std::numeric_limits<double>::max() ),
     _minNextPredictorTimeStepSize( std::numeric_limits<double>::max() ),
     _dofPerFace( numberOfVariables * power(nodesPerCoordinateAxis, DIMENSIONS - 1) ),
     _dofPerCellBoundary( DIMENSIONS_TIMES_TWO * _dofPerFace ),
     _dofPerCell( numberOfVariables * power(nodesPerCoordinateAxis, DIMENSIONS + 0) ),
     _fluxDofPerCell( _dofPerCell * (DIMENSIONS + 1) ),  // +1 for sources
     _spaceTimeDofPerCell( numberOfVariables * power(nodesPerCoordinateAxis, DIMENSIONS + 1) ),
     _spaceTimeFluxDofPerCell( _spaceTimeDofPerCell * (DIMENSIONS + 1) ),  // +1 for sources
     _dataPointsPerCell( (numberOfVariables+numberOfParameters) * power(nodesPerCoordinateAxis, DIMENSIONS + 0)) {
  // register tags with profiler
  for (const char* tag : tags) {
    _profiler->registerTag(tag);
  }
  
  for (const char* tag : deepProfilingTags) {
    _profiler->registerTag(tag); //TODO JMG only if using deepProfiling
  }

  CompressedDataHeap::getInstance().setName("compressed-data");
}

int exahype::solvers::ADERDGSolver::getUnknownsPerFace() const {
  return _dofPerFace;
}

int exahype::solvers::ADERDGSolver::getUnknownsPerCellBoundary() const {
  return _dofPerCellBoundary;
}

int exahype::solvers::ADERDGSolver::getUnknownsPerCell() const {
  return _dofPerCell;
}

int exahype::solvers::ADERDGSolver::getFluxUnknownsPerCell() const {
  return _fluxDofPerCell;
}

int exahype::solvers::ADERDGSolver::getSpaceTimeUnknownsPerCell() const {
  return _spaceTimeDofPerCell;
}

int exahype::solvers::ADERDGSolver::getSpaceTimeFluxUnknownsPerCell() const {
  return _spaceTimeFluxDofPerCell;
}

int exahype::solvers::ADERDGSolver::getDataPerFace() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS - 1);
}

int exahype::solvers::ADERDGSolver::getDataPerCellBoundary() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS - 1) * DIMENSIONS_TIMES_TWO;
}

int exahype::solvers::ADERDGSolver::getDataPerCell() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS + 0);
}

int exahype::solvers::ADERDGSolver::getSpaceTimeDataPerCell() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS + 1);
}

void exahype::solvers::ADERDGSolver::synchroniseTimeStepping(
    CellDescription& p) const {
  switch (_timeStepping) {
    case TimeStepping::Global:
      p.setPreviousCorrectorTimeStamp(_previousMinCorrectorTimeStamp);
      p.setPreviousCorrectorTimeStepSize(_previousMinCorrectorTimeStepSize);

      p.setCorrectorTimeStamp(_minCorrectorTimeStamp);
      p.setCorrectorTimeStepSize(_minCorrectorTimeStepSize);

      p.setPredictorTimeStamp(_minPredictorTimeStamp);
      p.setPredictorTimeStepSize(_minPredictorTimeStepSize);
      break;
    case TimeStepping::GlobalFixed:
      p.setPreviousCorrectorTimeStamp(_previousMinCorrectorTimeStamp);
      p.setPreviousCorrectorTimeStepSize(_previousMinCorrectorTimeStepSize);

      p.setCorrectorTimeStamp(_minCorrectorTimeStamp);
      p.setCorrectorTimeStepSize(_minCorrectorTimeStepSize);

      p.setPredictorTimeStamp(_minPredictorTimeStamp);
      p.setPredictorTimeStepSize(_minPredictorTimeStepSize);
      break;
  }
}

void exahype::solvers::ADERDGSolver::synchroniseTimeStepping(
      const int cellDescriptionsIndex,
      const int element) {
  synchroniseTimeStepping(Heap::getInstance().getData(cellDescriptionsIndex)[element]);
}

void exahype::solvers::ADERDGSolver::startNewTimeStep() {
  switch (_timeStepping) {
    case TimeStepping::Global:
      // n-1
      _previousMinCorrectorTimeStamp            = _minCorrectorTimeStamp;
      _previousMinCorrectorTimeStepSize         = _minCorrectorTimeStepSize;
      // n
      _minCorrectorTimeStamp    = _minPredictorTimeStamp;
      _minCorrectorTimeStepSize = _minPredictorTimeStepSize;
      // n+1
      _minPredictorTimeStamp    = _minPredictorTimeStamp + _minPredictorTimeStepSize;
      _minPredictorTimeStepSize = _minNextPredictorTimeStepSize;

      _minNextPredictorTimeStepSize = std::numeric_limits<double>::max();
      break;
    case TimeStepping::GlobalFixed:
      // n-1
      _previousMinCorrectorTimeStamp            = _minCorrectorTimeStamp;
      _previousMinCorrectorTimeStepSize         = _minCorrectorTimeStepSize;
      // n
      _minCorrectorTimeStamp    = _minPredictorTimeStamp;
      _minCorrectorTimeStepSize = _minPredictorTimeStepSize;
      // n+1
      _minPredictorTimeStamp    = _minPredictorTimeStamp + _minPredictorTimeStepSize;
      _minPredictorTimeStepSize = _minNextPredictorTimeStepSize;
      break;
  }

  _minCellSize     = _nextMinCellSize;
  _maxCellSize     = _nextMaxCellSize;
  _nextMinCellSize = std::numeric_limits<double>::max();
  _nextMaxCellSize = -std::numeric_limits<double>::max(); // "-", min

  setNextGridUpdateRequested();
}

void exahype::solvers::ADERDGSolver::zeroTimeStepSizes() {
  _minCorrectorTimeStepSize = 0;
  _minPredictorTimeStepSize = 0;
}

void exahype::solvers::ADERDGSolver::reconstructStandardTimeSteppingData() {
  //  _previousMinCorrectorTimeStepSize = _minCorrectorTimeStepSize; // TODO(Dominic): Should not necessary.
  _minPredictorTimeStamp    = _minCorrectorTimeStamp+_minCorrectorTimeStepSize;
  _minCorrectorTimeStamp    = _minPredictorTimeStamp;
  _minCorrectorTimeStepSize = _minPredictorTimeStepSize;

  assertionEquals(_minCorrectorTimeStamp,_minPredictorTimeStamp);
  assertionEquals(_minCorrectorTimeStepSize,_minPredictorTimeStepSize);
}

void exahype::solvers::ADERDGSolver::reinitialiseTimeStepData() {
  switch (_timeStepping) {
    case TimeStepping::Global:
      _minPredictorTimeStepSize = _minNextPredictorTimeStepSize;
      _minPredictorTimeStamp    = _minCorrectorTimeStamp+_minNextPredictorTimeStepSize;
      break;
    case TimeStepping::GlobalFixed:
      //do nothing
      break;
  }
}

void exahype::solvers::ADERDGSolver::rollbackToPreviousTimeStep() {

  //  logInfo("rollbackToPreviousTimeStep(...)", // for  debugging purposes

//  logInfo("rollbackToPreviousTimeStep(...)",
//          "\n[pre]"
//          << "\n,previousPreviousMinCorrectorTimeStepSize="<<_previousPreviousMinCorrectorTimeStepSize
//          << "\n,previousMinCorrectorTimeStepSize="<<_previousMinCorrectorTimeStepSize
//          << "\n,previousMinCorrectorTimeStamp="<<_previousMinCorrectorTimeStamp
//          << "\n,minCorrectorTimeStamp="<<_minCorrectorTimeStamp
//          << "\n,minCorrectorTimeStepSize="<<_minCorrectorTimeStepSize
//          <<" \n,minPredictorTimeStamp="<<_minPredictorTimeStamp
//          << "\n,minPredictorTimeStepSize="<<_minPredictorTimeStepSize
//    );


  switch (_timeStepping) {
    case TimeStepping::Global:
      _minNextPredictorTimeStepSize             = std::numeric_limits<double>::max();

      _minPredictorTimeStamp                    = _minCorrectorTimeStamp;
      _minPredictorTimeStepSize                 = _minCorrectorTimeStepSize;

      _minCorrectorTimeStamp                    = _previousMinCorrectorTimeStamp;
      _minCorrectorTimeStepSize                 = _previousMinCorrectorTimeStepSize;

      _previousMinCorrectorTimeStamp            = std::numeric_limits<double>::max();
      _previousMinCorrectorTimeStepSize         = std::numeric_limits<double>::max();
      break;
    case TimeStepping::GlobalFixed:
      _minPredictorTimeStamp                    = _minCorrectorTimeStamp;
      _minPredictorTimeStepSize                 = _minCorrectorTimeStepSize;

      _minCorrectorTimeStamp                    = _previousMinCorrectorTimeStamp;
      _minCorrectorTimeStepSize                 = _previousMinCorrectorTimeStepSize;

      _previousMinCorrectorTimeStamp            = std::numeric_limits<double>::max();
      _previousMinCorrectorTimeStepSize         = std::numeric_limits<double>::max();
      break;
  }

//  logInfo("rollbackToPreviousTimeStep(...)", // for  debugging purposes
//          "\n[post]"
//          << "\n,previousPreviousMinCorrectorTimeStepSize="<<_previousPreviousMinCorrectorTimeStepSize
//          << "\n,previousMinCorrectorTimeStepSize="<<_previousMinCorrectorTimeStepSize
//          << "\n,previousMinCorrectorTimeStamp="<<_previousMinCorrectorTimeStamp
//          << "\n,minCorrectorTimeStamp="<<_minCorrectorTimeStamp
//          << "\n,minCorrectorTimeStepSize="<<_minCorrectorTimeStepSize
//          <<" \n,minPredictorTimeStamp="<<_minPredictorTimeStamp
//          << "\n,minPredictorTimeStepSize="<<_minPredictorTimeStepSize
//  );


//  _minCellSize     = _nextMinCellSize;
//  _maxCellSize     = _nextMaxCellSize; // TODO(Dominic): What to do with these fields?
  _nextMinCellSize = std::numeric_limits<double>::max();
  _nextMaxCellSize = -std::numeric_limits<double>::max(); // "-", min
}

void exahype::solvers::ADERDGSolver::reconstructStandardTimeSteppingDataAfterRollback() {
  _minPredictorTimeStamp            = _minCorrectorTimeStamp;    // corrector time stamp is now the previous corrector time stamp
  _minPredictorTimeStepSize         = _minCorrectorTimeStepSize; // corrector time step size is now the previous corrector time step size
  _previousMinCorrectorTimeStepSize = std::numeric_limits<double>::max();
}

void exahype::solvers::ADERDGSolver::updateMinNextPredictorTimeStepSize(
    const double& minNextPredictorTimeStepSize) {
  switch (_timeStepping) {
    case TimeStepping::Global:
      _minNextPredictorTimeStepSize =
          std::min(_minNextPredictorTimeStepSize, minNextPredictorTimeStepSize);
      break;
    case TimeStepping::GlobalFixed:
      _minNextPredictorTimeStepSize =
          _minNextPredictorTimeStepSize == std::numeric_limits<double>::max()
              ? std::min(_minNextPredictorTimeStepSize,
                         minNextPredictorTimeStepSize)
              : _minNextPredictorTimeStepSize;
      break;
  }
}

double exahype::solvers::ADERDGSolver::getMinNextPredictorTimeStepSize() const {
  return _minNextPredictorTimeStepSize;
}

void exahype::solvers::ADERDGSolver::setMinCorrectorTimeStamp(
    double minCorrectorTimeStamp) {
  _minCorrectorTimeStamp = minCorrectorTimeStamp;
}

double exahype::solvers::ADERDGSolver::getMinCorrectorTimeStamp() const {
  return _minCorrectorTimeStamp;
}

void exahype::solvers::ADERDGSolver::setMinPredictorTimeStamp(
    double minPredictorTimeStamp) {
  _minPredictorTimeStamp = minPredictorTimeStamp;
}

double exahype::solvers::ADERDGSolver::getMinPredictorTimeStamp() const {
  return _minPredictorTimeStamp;
}

void exahype::solvers::ADERDGSolver::setMinCorrectorTimeStepSize(
    double minCorrectorTimeStepSize) {
  _minCorrectorTimeStepSize = minCorrectorTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getMinCorrectorTimeStepSize() const {
  return _minCorrectorTimeStepSize;
}

void exahype::solvers::ADERDGSolver::setMinPredictorTimeStepSize(
    double minPredictorTimeStepSize) {
  _minPredictorTimeStepSize = minPredictorTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getMinPredictorTimeStepSize() const {
  return _minPredictorTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getPreviousMinCorrectorTimeStepSize() const {
  return _previousMinCorrectorTimeStepSize;
}

void exahype::solvers::ADERDGSolver::setPreviousMinCorrectorTimeStamp(double value) {
  _previousMinCorrectorTimeStamp = value;
}

double exahype::solvers::ADERDGSolver::getPreviousMinCorrectorTimeStamp() const {
  return _previousMinCorrectorTimeStamp;
}

void exahype::solvers::ADERDGSolver::setPreviousMinCorrectorTimeStepSize(double value) {
  _previousMinCorrectorTimeStepSize = value;
}

double exahype::solvers::ADERDGSolver::getMinTimeStamp() const {
  return getMinCorrectorTimeStamp();
}

double exahype::solvers::ADERDGSolver::getMinTimeStepSize() const {
  return getMinCorrectorTimeStepSize();
}

double exahype::solvers::ADERDGSolver::getMinNextTimeStepSize() const {
  return getMinNextPredictorTimeStepSize();
}

void exahype::solvers::ADERDGSolver::updateMinNextTimeStepSize( double value ) {
  updateMinNextPredictorTimeStepSize(value);
}

void exahype::solvers::ADERDGSolver::initSolver(const double timeStamp, const tarch::la::Vector<DIMENSIONS,double>& boundingBox) {
  _coarsestMeshLevel = exahype::solvers::Solver::computeMeshLevel(_maximumMeshSize,boundingBox[0]);

  setPreviousMinCorrectorTimeStepSize(0.0);
  setMinCorrectorTimeStepSize(0.0);
  setMinPredictorTimeStepSize(0.0);

  setPreviousMinCorrectorTimeStamp(timeStamp);
  setMinCorrectorTimeStamp(timeStamp);
  setMinPredictorTimeStamp(timeStamp);
}

void exahype::solvers::ADERDGSolver::initFusedSolverTimeStepSizes() {
  setPreviousMinCorrectorTimeStepSize(getMinPredictorTimeStepSize());
  setMinCorrectorTimeStepSize(getMinPredictorTimeStepSize());
  setMinPredictorTimeStepSize(getMinPredictorTimeStepSize());
}

bool exahype::solvers::ADERDGSolver::isValidCellDescriptionIndex(
      const int cellDescriptionsIndex) const {
    return Heap::getInstance().isValidIndex(cellDescriptionsIndex);
  }

int exahype::solvers::ADERDGSolver::tryGetElement(
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  if (Heap::getInstance().isValidIndex(cellDescriptionsIndex)) {
    int element=0;
    for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
      if (p.getSolverNumber()==solverNumber) {
        return element;
      }
      ++element;
    }
  }
  return NotFound;
}

exahype::solvers::Solver::SubcellPosition
exahype::solvers::ADERDGSolver::computeSubcellPositionOfCellOrAncestor(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& cellDescription =
      getCellDescription(cellDescriptionsIndex,element);

  return
      exahype::amr::computeSubcellPositionOfCellOrAncestor
      <CellDescription,Heap>(cellDescription);
}

///////////////////////////////////
// CELL-LOCAL MESH REFINEMENT
///////////////////////////////////
void exahype::solvers::ADERDGSolver::ensureConsistencyOfParentIndex(
    CellDescription& fineGridCellDescription,
    const int coarseGridCellDescriptionsIndex,
    const int solverNumber) {
  fineGridCellDescription.setParentIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
  int coarseGridElement = tryGetElement(coarseGridCellDescriptionsIndex,solverNumber);
  if (coarseGridElement!=exahype::solvers::Solver::NotFound) {
    fineGridCellDescription.setParentIndex(coarseGridCellDescriptionsIndex);
    // In this case, we are not at a master worker boundary only our parent is.
    #ifdef Parallel
    fineGridCellDescription.setHasToHoldDataForMasterWorkerCommunication(false);
    #endif
  }
}

bool exahype::solvers::ADERDGSolver::markForRefinement(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber) {
  // Fine grid cell based uniform mesh refinement.
  const int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);

  // Fine grid cell based adaptive mesh refinement operations.
  if (fineGridCellElement!=exahype::solvers::Solver::NotFound) {
    CellDescription& fineGridCellDescription = Heap::getInstance().getData(
        fineGridCell.getCellDescriptionsIndex())[fineGridCellElement];

    if (fineGridCellDescription.getLevel()<_coarsestMeshLevel+_maximumAdaptiveMeshDepth) {
      switch (determineLimiterStatus(fineGridCellDescription)) {
        // assertion: mergedLimiterStatus has been unified
        case CellDescription::LimiterStatus::NeighbourOfTroubled3:
        case CellDescription::LimiterStatus::NeighbourOfTroubled4:
        case CellDescription::LimiterStatus::Ok: {
          #ifdef Parallel
          ensureConsistencyOfParentIndex(fineGridCellDescription,coarseGridCell.getCellDescriptionsIndex(),solverNumber);
          #endif
          #if defined(Asserts) || defined(Debug)
          int coarseGridCellElement =
              tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
          #endif
          assertion3(coarseGridCellElement==exahype::solvers::Solver::NotFound ||
                     fineGridCellDescription.getParentIndex()==coarseGridCell.getCellDescriptionsIndex(),
                     fineGridCellDescription.toString(),fineGridCell.toString(),
                     coarseGridCell.toString()); // see mergeCellDescriptionsWithRemoteData.

          // marking for refinement
          return markForRefinement(fineGridCellDescription);
        }
        default:
          return false;
      }
    }
  }

  return false;
}

///////////////////////////////////
// CELL-LOCAL MESH REFINEMENT
///////////////////////////////////
bool exahype::solvers::ADERDGSolver::updateStateInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber) {
  bool refineFineGridCell = false;

  // Fine grid cell based uniform mesh refinement.
  int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (fineGridCellElement==exahype::solvers::Solver::NotFound &&
      tarch::la::allSmallerEquals(fineGridVerticesEnumerator.getCellSize(),getMaximumMeshSize()) &&
      tarch::la::allGreater(coarseGridVerticesEnumerator.getCellSize(),getMaximumMeshSize())) {
    addNewCell(fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
               multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex,
               solverNumber); // TODO(Dominic): Can directly refine if we directly evaluate initial conditions here.
  // Fine grid cell based adaptive mesh refinement operations.
  } else if (fineGridCellElement!=exahype::solvers::Solver::NotFound) {
    CellDescription& fineGridCellDescription = Heap::getInstance().getData(
        fineGridCell.getCellDescriptionsIndex())[fineGridCellElement];

    #ifdef Parallel
    ensureConsistencyOfParentIndex(fineGridCellDescription,coarseGridCell.getCellDescriptionsIndex(),solverNumber);
    #endif
    #if defined(Asserts) || defined(Debug)
    int coarseGridCellElement =
        tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
    #endif
    assertion3(coarseGridCellElement==exahype::solvers::Solver::NotFound ||
        fineGridCellDescription.getParentIndex()==coarseGridCell.getCellDescriptionsIndex(),
        fineGridCellDescription.toString(),fineGridCell.toString(),
        coarseGridCell.toString()); // see mergeCellDescriptionsWithRemoteData.

    // actions requiring adjacency info
    const tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int>&
      indicesAdjacentToFineGridVertices =
          VertexOperations::readCellDescriptionsIndex(
                          fineGridVerticesEnumerator,fineGridVertices);

    if (multiscalelinkedcell::adjacencyInformationIsConsistent(
        indicesAdjacentToFineGridVertices)) {
      const tarch::la::Vector<THREE_POWER_D, int> neighbourCellDescriptionIndices =
          multiscalelinkedcell::getIndicesAroundCell(
              indicesAdjacentToFineGridVertices);

      #ifdef Parallel
      fineGridCellDescription.setAdjacentToRemoteRank(
          exahype::Cell::isAdjacentToRemoteRank(fineGridVertices,fineGridVerticesEnumerator));
      #endif

      // marking for augmentation
      refineFineGridCell |= // TODO(Dominic): Change to the template version.
          markForAugmentation(
              fineGridCellDescription,
              neighbourCellDescriptionIndices,
              fineGridCell.isAssignedToRemoteRank());
    }
  }

  // Coarse grid cell based adaptive mesh refinement operations.
  // Add new cells to the grid and veto erasing or deaugmenting children
  // requests if there are cells on the fine level.
  int coarseGridCellElement =
      tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
  if (coarseGridCellElement!=exahype::solvers::Solver::NotFound) {
    CellDescription& coarseGridCellDescription = getCellDescription(
        coarseGridCell.getCellDescriptionsIndex(),coarseGridCellElement);

    vetoErasingOrDeaugmentingChildrenRequest(
        coarseGridCellDescription,
        fineGridCell.getCellDescriptionsIndex());

    // TODO(Dominic): Pass limiter status flag down to the new cell
    addNewDescendantIfAugmentingRequested(
            fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
            coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex());
    addNewCellIfRefinementRequested(
        fineGridCell,fineGridVertices,fineGridVerticesEnumerator,fineGridPositionOfCell,
        coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex(),
        initialGrid);
  }

  return refineFineGridCell;
}

bool exahype::solvers::ADERDGSolver::markForRefinement(
    CellDescription& fineGridCellDescription) {
  bool refineFineGridCell = false;
  bool vetoErasing        = true;

  double* solution = 0;
  exahype::solvers::Solver::RefinementControl refinementControl;

  switch (fineGridCellDescription.getType()) {
    case CellDescription::Cell:
      switch (fineGridCellDescription.getRefinementEvent()) {
        case CellDescription::RefiningRequested:
          refineFineGridCell = true;
          break;
        case CellDescription::None:
        case CellDescription::AugmentingRequested:
          solution = DataHeap::getInstance().getData(fineGridCellDescription.getSolution()).data();
          refinementControl =
              refinementCriterion(
                  solution,fineGridCellDescription.getOffset()+0.5*fineGridCellDescription.getSize(),
                  fineGridCellDescription.getSize(),
                  fineGridCellDescription.getCorrectorTimeStamp()+fineGridCellDescription.getCorrectorTimeStepSize(),
                  fineGridCellDescription.getLevel());

          switch (refinementControl) {
            case exahype::solvers::Solver::RefinementControl::Refine:
              fineGridCellDescription.setRefinementEvent(CellDescription::RefiningRequested);
              refineFineGridCell = true;
              break;
            case exahype::solvers::Solver::RefinementControl::Erase:
              vetoErasing = false;
              break;
            default:
              break;
          }
          break;
        default:
          break;
      }
      break;
    case CellDescription::Ancestor:
    case CellDescription::EmptyAncestor:
      if (
          fineGridCellDescription.getRefinementEvent()==CellDescription::None
          && !fineGridCellDescription.getNewlyCreated() // we do not erase newly created Ancestors
      ) {
        fineGridCellDescription.setRefinementEvent(CellDescription::ErasingChildrenRequested);

        /*  TODO(Dominic): Add to docu:
         * If this refinement event is set,
         * the parent Ancestor asks its
         * children if they want to be erased. If not,
         * the children change the RefinementEvent
         * of the parent to None. If so,
         * they leave the parent's RefinementEvent
         * unchanged.
         */
      }
      break;
    default:
      // do nothing
      break;
  }

  if (vetoErasing) {
    int coarseGridCellElement = tryGetElement(fineGridCellDescription.getParentIndex(),
                                              fineGridCellDescription.getSolverNumber());
    if (coarseGridCellElement!=exahype::solvers::Solver::NotFound) {
      auto& coarseGridCellDescription = getCellDescription(fineGridCellDescription.getParentIndex(),
                                                           coarseGridCellElement);

      if (coarseGridCellDescription.getRefinementEvent()==CellDescription::ErasingChildrenRequested ||
          coarseGridCellDescription.getRefinementEvent()==CellDescription::ChangeChildrenToDescendantsRequested) {
        coarseGridCellDescription.setRefinementEvent(CellDescription::None);

        assertion1(coarseGridCellDescription.getType()==CellDescription::EmptyAncestor ||
                   coarseGridCellDescription.getType()==CellDescription::Ancestor,
                   coarseGridCellDescription.toString());
      }
    }
  }

  return refineFineGridCell;
}

bool exahype::solvers::ADERDGSolver::markForAugmentation(
    CellDescription& fineGridCellDescription,
    const tarch::la::Vector<THREE_POWER_D, int>& neighbourCellDescriptionIndices,
    const bool assignedToRemoteRank) {
  int coarseGridElement = tryGetElement(
      fineGridCellDescription.getParentIndex(),
      fineGridCellDescription.getSolverNumber());
  // First check if we can set the deaugmenting children requested triggered event of the coarse grid cell
  // to a "real" deaugmenting children requested event.
  if (coarseGridElement!=exahype::solvers::Solver::NotFound) {
    CellDescription& coarseGridCellDescription =
        getCellDescription(fineGridCellDescription.getParentIndex(),coarseGridElement);
    if (coarseGridCellDescription.getRefinementEvent()==CellDescription::DeaugmentingChildrenRequestedTriggered) {
      coarseGridCellDescription.setRefinementEvent(CellDescription::DeaugmentingChildrenRequested);
    }
  }

  // Then do fine grid cell stuff
  bool refineFineGridCell = false;
  bool vetoDeaugmenting   = true;

  const exahype::solvers::Solver::AugmentationControl augmentationControl =
      exahype::amr::augmentationCriterion<CellDescription,Heap>(
          fineGridCellDescription.getSolverNumber(),fineGridCellDescription.getType(),
          fineGridCellDescription.getLevel(),
          neighbourCellDescriptionIndices);

  // Check if we have allocated enough memory.
  ensureOnlyNecessaryMemoryIsAllocated(
      fineGridCellDescription,augmentationControl,assignedToRemoteRank);

  // Further augment or deaugment cells and descendants if no other event
  // or an augmentation event has been triggered.
  switch (fineGridCellDescription.getRefinementEvent()) {
    case CellDescription::AugmentingRequested: // TODO(Dominic): Add to docu that the mergeWithNeighbourCall might set this.
      refineFineGridCell = true;
      break;
    case CellDescription::None:
      switch (fineGridCellDescription.getType()) {
        case CellDescription::Cell:
          fineGridCellDescription.setRefinementEvent(CellDescription::DeaugmentingChildrenRequestedTriggered);
          switch (augmentationControl) {
            case exahype::solvers::Solver::AugmentationControl::NextToAncestor:
            case exahype::solvers::Solver::AugmentationControl::NextToCellAndAncestor:
              fineGridCellDescription.setRefinementEvent(CellDescription::AugmentingRequested);
              refineFineGridCell = true;
              break;
            default:
              break;
          }
          break;
        case CellDescription::Descendant:
        case CellDescription::EmptyDescendant:
          fineGridCellDescription.setRefinementEvent(CellDescription::DeaugmentingChildrenRequestedTriggered);
          switch (augmentationControl) {
            case exahype::solvers::Solver::AugmentationControl::NextToAncestor:
            case exahype::solvers::Solver::AugmentationControl::NextToCellAndAncestor:
              fineGridCellDescription.setRefinementEvent(CellDescription::AugmentingRequested);
              refineFineGridCell = true;
              break;
            case exahype::solvers::Solver::AugmentationControl::Default:
              vetoDeaugmenting = false;
              break;
            default:
              break;
          }
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }

  // And check if we can reset the deaugmenting children request of the parent.
  if (vetoDeaugmenting) {
    int coarseGridCellElement = tryGetElement(fineGridCellDescription.getParentIndex(),
                                              fineGridCellDescription.getSolverNumber());
    if (coarseGridCellElement!=exahype::solvers::Solver::NotFound) {
      auto& coarseGridCellDescription = getCellDescription(fineGridCellDescription.getParentIndex(),
                                                           coarseGridCellElement);
      if (coarseGridCellDescription.getRefinementEvent()==CellDescription::DeaugmentingChildrenRequested) {
        coarseGridCellDescription.setRefinementEvent(CellDescription::None);

        assertion1(coarseGridCellDescription.getType()==CellDescription::Cell ||
                   coarseGridCellDescription.getType()==CellDescription::EmptyDescendant ||
                   coarseGridCellDescription.getType()==CellDescription::Descendant,
                   coarseGridCellDescription.toString());
      }
    }
  }

  return refineFineGridCell;
}

void exahype::solvers::ADERDGSolver::ensureOnlyNecessaryMemoryIsAllocated(
    CellDescription& fineGridCellDescription,
    const exahype::solvers::Solver::AugmentationControl& augmentationControl,
    const bool assignedToRemoteRank) {
  switch (fineGridCellDescription.getType()) {
    case CellDescription::Ancestor:
    case CellDescription::EmptyAncestor:
      fineGridCellDescription.setType(CellDescription::EmptyAncestor);
      #ifdef Parallel
      if (fineGridCellDescription.getHasToHoldDataForNeighbourCommunication()
          || fineGridCellDescription.getHasToHoldDataForMasterWorkerCommunication()) { // TODO(Dominic): Add to docu.
        fineGridCellDescription.setType(CellDescription::Ancestor);
      }
      #endif

      if (augmentationControl==exahype::solvers::Solver::AugmentationControl::NextToCell ||
          augmentationControl==exahype::solvers::Solver::AugmentationControl::NextToCellAndAncestor) {
        fineGridCellDescription.setType(CellDescription::Ancestor);
      }
      ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);
      ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
      break;
    case CellDescription::Descendant:
    case CellDescription::EmptyDescendant:
      fineGridCellDescription.setType(CellDescription::EmptyDescendant);

      #ifdef Parallel
      if (fineGridCellDescription.getHasToHoldDataForNeighbourCommunication()
          || fineGridCellDescription.getHasToHoldDataForMasterWorkerCommunication()) { // TODO(Dominic): Add to docu.
        fineGridCellDescription.setType(CellDescription::Descendant);
      }
      #endif

      if (augmentationControl==exahype::solvers::Solver::AugmentationControl::NextToCell ||
          augmentationControl==exahype::solvers::Solver::AugmentationControl::NextToCellAndAncestor) {
        fineGridCellDescription.setType(CellDescription::Descendant);
      }
      ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
      ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);
      break;
    default:
      break;
  }

  #ifdef Parallel
  fineGridCellDescription.setHasToHoldDataForNeighbourCommunication(false);
  #endif
}

void exahype::solvers::ADERDGSolver::vetoErasingOrDeaugmentingChildrenRequest(
    CellDescription& coarseGridCellDescription,
    const int fineGridCellDescriptionsIndex) {
  int coarseGridCellParentElement = tryGetElement(coarseGridCellDescription.getParentIndex(),
                                                  coarseGridCellDescription.getSolverNumber());
  int fineGridCellElement = tryGetElement(fineGridCellDescriptionsIndex,
                                          coarseGridCellDescription.getSolverNumber());
  if (fineGridCellElement!=exahype::solvers::Solver::NotFound &&
      coarseGridCellParentElement!=exahype::solvers::Solver::NotFound) {
    CellDescription& coarseGridCellDescriptionParent =
        getCellDescription(coarseGridCellDescription.getParentIndex(),coarseGridCellParentElement);

        switch (coarseGridCellDescriptionParent.getRefinementEvent()) {
          case CellDescription::DeaugmentingChildrenRequested:
            assertion1(coarseGridCellDescription.getType()==CellDescription::EmptyDescendant ||
                       coarseGridCellDescription.getType()==CellDescription::Descendant,coarseGridCellDescription.toString());
            coarseGridCellDescriptionParent.setRefinementEvent(CellDescription::None); // TODO(DOminic): Source of bug.
            break;
          case CellDescription::ErasingChildrenRequested:
            assertion1(coarseGridCellDescription.getType()==CellDescription::Cell,
                       coarseGridCellDescription.toString());

            coarseGridCellDescriptionParent.setRefinementEvent(
                CellDescription::ChangeChildrenToDescendantsRequested);
            break;
          default:
            break;
        }
  }
}

void exahype::solvers::ADERDGSolver::addNewCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const int coarseGridCellDescriptionsIndex,
    const int solverNumber) {
  fineGridCell.addNewCellDescription(
              solverNumber,
              CellDescription::Cell,
              CellDescription::None,
              fineGridVerticesEnumerator.getLevel(),
              coarseGridCellDescriptionsIndex,
              fineGridVerticesEnumerator.getCellSize(),
              fineGridVerticesEnumerator.getVertexPosition());
  int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  CellDescription& fineGridCellDescription =
      getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridCellElement);

  ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
  exahype::Cell::determineInsideAndOutsideFaces(
            fineGridCellDescription,
            fineGridVertices,
            fineGridVerticesEnumerator);
}

void exahype::solvers::ADERDGSolver::addNewDescendantIfAugmentingRequested(
     exahype::Cell& fineGridCell,
     exahype::Vertex* const fineGridVertices,
     const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
     CellDescription& coarseGridCellDescription,
     const int coarseGridCellDescriptionsIndex) {
  if (coarseGridCellDescription.getRefinementEvent()==CellDescription::AugmentingRequested ||
      coarseGridCellDescription.getRefinementEvent()==CellDescription::Augmenting) {
    assertion1(coarseGridCellDescription.getType()==CellDescription::Cell ||
               coarseGridCellDescription.getType()==CellDescription::EmptyDescendant ||
               coarseGridCellDescription.getType()==CellDescription::Descendant,
               coarseGridCellDescription.toString());
    int fineGridElement =
        tryGetElement(fineGridCell.getCellDescriptionsIndex(),
                      coarseGridCellDescription.getSolverNumber());

    coarseGridCellDescription.setRefinementEvent(CellDescription::None);
    if (fineGridElement==exahype::solvers::Solver::NotFound) {
      fineGridCell.addNewCellDescription( // (EmptyDescendant),None
          coarseGridCellDescription.getSolverNumber(),
          CellDescription::EmptyDescendant,
          CellDescription::None,
          fineGridVerticesEnumerator.getLevel(),
          coarseGridCellDescriptionsIndex,
          fineGridVerticesEnumerator.getCellSize(),
          fineGridVerticesEnumerator.getVertexPosition());
      int fineGridElement =
          tryGetElement(fineGridCell.getCellDescriptionsIndex(),
                        coarseGridCellDescription.getSolverNumber());
      CellDescription& fineGridCellDescription =
          getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridElement);
      exahype::Cell::determineInsideAndOutsideFaces(
          fineGridCellDescription,
          fineGridVertices,
          fineGridVerticesEnumerator);

      coarseGridCellDescription.setRefinementEvent(CellDescription::Augmenting);
    } else if (fineGridElement!=exahype::solvers::Solver::NotFound &&
               coarseGridCellDescription.getRefinementEvent()==CellDescription::AugmentingRequested){
      /**
       * Reset an augmentation request if the child cell does hold
       * a Descendant or EmptyDescendant cell description with
       * the same solver number.
       *
       * This scenario occurs if an augmentation request is triggered in
       * enterCell() or mergeWithNeighbourMetadata(...).
       *
       * A similar scenario can never occur for refinement requests
       * since only cell descriptions of type Cell can be refined.
       * Ancestors and EmptyAncestors can never request refinement.
       * TODO(Dominic): Add to docu.
       */
      coarseGridCellDescription.setRefinementEvent(CellDescription::None);

      #if defined(Debug) || defined(Asserts)
      CellDescription& fineGridCellDescription =
                getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridElement);
      #endif

      assertion1(fineGridCellDescription.getType()==CellDescription::EmptyDescendant ||
                 fineGridCellDescription.getType()==CellDescription::Descendant,
                 fineGridCellDescription.toString());
    }
  }
}

void exahype::solvers::ADERDGSolver::addNewCellIfRefinementRequested(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS,int>& fineGridPositionOfCell,
    CellDescription& coarseGridCellDescription,
    const int coarseGridCellDescriptionsIndex,
    const bool initialGrid) {
  if (coarseGridCellDescription.getRefinementEvent()==CellDescription::RefiningRequested ||
      coarseGridCellDescription.getRefinementEvent()==CellDescription::Refining) {
    assertion1(coarseGridCellDescription.getType()==CellDescription::Cell,
               coarseGridCellDescription.toString());
    int fineGridCellElement =
        tryGetElement(fineGridCell.getCellDescriptionsIndex(),
                      coarseGridCellDescription.getSolverNumber());

    if (fineGridCellElement==exahype::solvers::Solver::NotFound) {
      addNewCell(fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
                 coarseGridCellDescriptionsIndex,
                 coarseGridCellDescription.getSolverNumber());
      fineGridCellElement =
          tryGetElement(fineGridCell.getCellDescriptionsIndex(),
                        coarseGridCellDescription.getSolverNumber());
      CellDescription& fineGridCellDescription =
          getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridCellElement);
      prolongateVolumeData(
          fineGridCellDescription,coarseGridCellDescription,fineGridPositionOfCell,initialGrid);
    } else {
      CellDescription& fineGridCellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridCellElement);
      assertion2(fineGridCellDescription.getType()==CellDescription::EmptyDescendant
                 || fineGridCellDescription.getType()==CellDescription::Descendant,
                 fineGridCellDescription.toString(),coarseGridCellDescription.toString());
      assertion2(fineGridCellDescription.getParentIndex()==coarseGridCellDescriptionsIndex,
                 fineGridCellDescription.toString(),coarseGridCellDescriptionsIndex);
      fineGridCellDescription.setType(CellDescription::Cell);
      fineGridCellDescription.setRefinementEvent(CellDescription::None);
      ensureNecessaryMemoryIsAllocated(fineGridCellDescription);

      prolongateVolumeData(
          fineGridCellDescription,coarseGridCellDescription,fineGridPositionOfCell,initialGrid);
    }

    coarseGridCellDescription.setRefinementEvent(CellDescription::Refining);
  }
}

void exahype::solvers::ADERDGSolver::prolongateVolumeData(
    CellDescription&       fineGridCellDescription,
    const CellDescription& coarseGridCellDescription,
  const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
  const bool initialGrid) {
  const int levelFine = fineGridCellDescription.getLevel();
  const int levelCoarse = coarseGridCellDescription.getLevel();
  assertion(levelCoarse < levelFine);

  // current solution
  double* solutionFine   = DataHeap::getInstance().getData(
      fineGridCellDescription.getSolution()).data();
  double* solutionCoarse = DataHeap::getInstance().getData(
      coarseGridCellDescription.getSolution()).data();
  volumeUnknownsProlongation(
      solutionFine,solutionCoarse,
      levelCoarse,levelFine,
      subcellIndex);

  // previous solution
  assertion(DataHeap::getInstance().isValidIndex(fineGridCellDescription.getPreviousSolution()));
  double* previousSolutionFine   = DataHeap::getInstance().getData(
      fineGridCellDescription.getPreviousSolution()).data();
  double* previousSolutionCoarse = DataHeap::getInstance().getData(
      coarseGridCellDescription.getPreviousSolution()).data();
  volumeUnknownsProlongation(
      previousSolutionFine,previousSolutionCoarse,
      levelCoarse,levelFine,
      subcellIndex);

  fineGridCellDescription.setCorrectorTimeStamp(coarseGridCellDescription.getCorrectorTimeStamp());
  fineGridCellDescription.setPredictorTimeStamp(coarseGridCellDescription.getPredictorTimeStamp());
  fineGridCellDescription.setCorrectorTimeStepSize(coarseGridCellDescription.getCorrectorTimeStepSize());
  fineGridCellDescription.setPredictorTimeStepSize(coarseGridCellDescription.getPredictorTimeStepSize());

  // TODO Dominic: This is a little inconsistent since I orignially tried to hide
  // the limiting from the pure ADER-DG scheme
//  fineGridCellDescription.setPreviousLimiterStatus(coarseGridCellDescription.getPreviousLimiterStatus());
//  const CellDescription::LimiterStatus limiterStatus = determineLimiterStatus(fineGridCellDescription);
//
//  for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
//    fineGridCellDescription.setLimiterStatus(i,limiterStatus);
//  }

  fineGridCellDescription.setLimiterStatus(CellDescription::LimiterStatus::Ok);
  writeLimiterStatusOnBoundary(fineGridCellDescription);

  // TODO Dominic: During the inital mesh build where we only refine
  // according to the PAD, we don't want to have a too broad refined area.
  // We thus do not flag children cells with troubled
  if (!initialGrid) {
    if (coarseGridCellDescription.getLimiterStatus()==CellDescription::LimiterStatus::Troubled) {
      fineGridCellDescription.setLimiterStatus(CellDescription::LimiterStatus::Troubled);
      writeLimiterStatusOnBoundary(fineGridCellDescription);
    }
  }
}

bool exahype::solvers::ADERDGSolver::attainedStableState(
        exahype::Cell& fineGridCell,
        exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        const int solverNumber) const {
  assertion(Heap::getInstance().isValidIndex(fineGridCell.getCellDescriptionsIndex()));

  const int element = tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (element!=exahype::solvers::Solver::NotFound) {
    CellDescription& cellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),element);

      const tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int>&
        indicesAdjacentToFineGridVertices =
            VertexOperations::readCellDescriptionsIndex(
                fineGridVerticesEnumerator,fineGridVertices);

      return (cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::None
              && multiscalelinkedcell::adjacencyInformationIsConsistent(indicesAdjacentToFineGridVertices));
  }

  return true;
}

bool exahype::solvers::ADERDGSolver::updateStateInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
     const int solverNumber) {
  const int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (fineGridCellElement!=exahype::solvers::Solver::NotFound) {
    // TODO(Dominic): Only allow erasing if limiter status is Ok

    CellDescription& fineGridCellDescription = getCellDescription(
            fineGridCell.getCellDescriptionsIndex(),fineGridCellElement);

    // start or finish collective operations
    startOrFinishCollectiveRefinementOperations(fineGridCellDescription);

    // ensure again that we have allocated enough memory
    const tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int>&
        indicesAdjacentToFineGridVertices =
            VertexOperations::readCellDescriptionsIndex(
                            fineGridVerticesEnumerator,fineGridVertices);
    if (multiscalelinkedcell::adjacencyInformationIsConsistent(
        indicesAdjacentToFineGridVertices)) {
      const tarch::la::Vector<THREE_POWER_D, int> neighbourCellDescriptionIndices =
          multiscalelinkedcell::getIndicesAroundCell(
              indicesAdjacentToFineGridVertices);
      const exahype::solvers::Solver::AugmentationControl augmentationControl =
          exahype::amr::augmentationCriterion<CellDescription,Heap>(
              fineGridCellDescription.getSolverNumber(),fineGridCellDescription.getType(),
              fineGridCellDescription.getLevel(),
              neighbourCellDescriptionIndices);
      ensureOnlyNecessaryMemoryIsAllocated(
          fineGridCellDescription,augmentationControl,fineGridCell.isAssignedToRemoteRank());
    }

    const int coarseGridCellElement =
        tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
    if (coarseGridCellElement!=exahype::solvers::Solver::NotFound) {
      assertion3(fineGridCellDescription.getParentIndex()==coarseGridCell.getCellDescriptionsIndex(),
                 fineGridCellDescription.toString(),fineGridCell.toString(),
                 coarseGridCell.toString()); // see mergeCellDescriptionsWithRemoteData.

      CellDescription& coarseGridCellDescription = getCellDescription(
          fineGridCellDescription.getParentIndex(),coarseGridCellElement);
      assertion1(fineGridCellDescription.getSolverNumber()==
          coarseGridCellDescription.getSolverNumber(),
                     fineGridCellDescription.toString());

      bool eraseOfFineGridCellRequested =
          eraseCellDescriptionIfNecessary(
              fineGridCell.getCellDescriptionsIndex(),
              fineGridCellElement,
              fineGridPositionOfCell,
              coarseGridCellDescription);

      return eraseOfFineGridCellRequested;
    }
  }

  return false;
}

void exahype::solvers::ADERDGSolver::startOrFinishCollectiveRefinementOperations(
     CellDescription& fineGridCellDescription) {
  switch (fineGridCellDescription.getRefinementEvent()) {
    case CellDescription::Refining:
      assertion1(fineGridCellDescription.getType()==CellDescription::Cell,
                 fineGridCellDescription.toString());
      fineGridCellDescription.setType(CellDescription::EmptyAncestor);
      fineGridCellDescription.setRefinementEvent(CellDescription::None);
      fineGridCellDescription.setNewlyCreated(true);
      ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);

      // TODO(Dominic): Add to docu
      fineGridCellDescription.setLimiterStatus(
          std::max(0, fineGridCellDescription.getLimiterStatus()-2) );
      writeLimiterStatusOnBoundary(fineGridCellDescription);
      break;
    case CellDescription::ChangeChildrenToDescendantsRequested:
      fineGridCellDescription.setType(CellDescription::Cell);
      ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
      fineGridCellDescription.setRefinementEvent(CellDescription::ChangeChildrenToDescendants);
      break;
    case CellDescription::ChangeChildrenToDescendants:
      fineGridCellDescription.setRefinementEvent(CellDescription::None);
      break;
    case CellDescription::ErasingChildrenRequested:
      fineGridCellDescription.setType(CellDescription::Cell);
      ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
      fineGridCellDescription.setRefinementEvent(CellDescription::ErasingChildren);
      break;
    case CellDescription::ErasingChildren:
      fineGridCellDescription.setRefinementEvent(CellDescription::None);
      break;
    case CellDescription::Augmenting:
      fineGridCellDescription.setRefinementEvent(CellDescription::None);
      break;
    case CellDescription::DeaugmentingChildrenRequestedTriggered: // This means we are actually a leaf cell of this solver.
      fineGridCellDescription.setRefinementEvent(CellDescription::None);
      break;
    case CellDescription::DeaugmentingChildrenRequested:
      fineGridCellDescription.setRefinementEvent(CellDescription::DeaugmentingChildren);
      break;
    case CellDescription::DeaugmentingChildren:
      fineGridCellDescription.setRefinementEvent(CellDescription::None);
      break;
    default:
      break;
  }
}

bool exahype::solvers::ADERDGSolver::eraseCellDescriptionIfNecessary(
    const int cellDescriptionsIndex,
    const int fineGridCellElement,
    const tarch::la::Vector<DIMENSIONS,int>& fineGridPositionOfCell,
    CellDescription& coarseGridCellDescription) {
  if (coarseGridCellDescription.getRefinementEvent()==CellDescription::ErasingChildren) {
    CellDescription& fineGridCellDescription = getCellDescription(
          cellDescriptionsIndex,fineGridCellElement);
    // restrict values.
    restrictVolumeData(
        coarseGridCellDescription,
        fineGridCellDescription,
        fineGridPositionOfCell);
    // TODO(Dominic): Reconsider for anarchic time stepping.
    coarseGridCellDescription.setCorrectorTimeStamp(fineGridCellDescription.getCorrectorTimeStamp());
    coarseGridCellDescription.setPredictorTimeStamp(fineGridCellDescription.getPredictorTimeStamp());

    // erase cell description // or change to descendant
    fineGridCellDescription.setType(CellDescription::Erased);
    ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);

    getCellDescriptions(cellDescriptionsIndex).erase(
        getCellDescriptions(cellDescriptionsIndex).begin()+fineGridCellElement);

    return false;
  } else if (coarseGridCellDescription.getRefinementEvent()==CellDescription::ChangeChildrenToDescendants) {
    CellDescription& fineGridCellDescription = getCellDescription(
        cellDescriptionsIndex,fineGridCellElement);
    // restrict values.
    restrictVolumeData(
        coarseGridCellDescription,
        fineGridCellDescription,
        fineGridPositionOfCell);
    // TODO(Dominic): Reconsider for anarchic time stepping.
    coarseGridCellDescription.setCorrectorTimeStamp(fineGridCellDescription.getCorrectorTimeStamp());
    coarseGridCellDescription.setPredictorTimeStamp(fineGridCellDescription.getPredictorTimeStamp());

    // erase cell description // or change to descendant
    fineGridCellDescription.setType(CellDescription::EmptyDescendant);
    ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);

    return false;
  } else if (coarseGridCellDescription.getRefinementEvent()==CellDescription::DeaugmentingChildren) {
    CellDescription& fineGridCellDescription = getCellDescription(
          cellDescriptionsIndex,fineGridCellElement);

    fineGridCellDescription.setType(CellDescription::Erased);
    ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);

    getCellDescriptions(cellDescriptionsIndex).erase(
        getCellDescriptions(cellDescriptionsIndex).begin()+fineGridCellElement);
    return true;
  }

  return false;
}

void exahype::solvers::ADERDGSolver::restrictVolumeData(
    CellDescription&       coarseGridCellDescription,
    const CellDescription& fineGridCellDescription,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  const int levelFine = fineGridCellDescription.getLevel();
  const int levelCoarse = coarseGridCellDescription.getLevel();
  assertion(levelCoarse < levelFine);

  double* luhFine   = DataHeap::getInstance().getData(
      fineGridCellDescription.getSolution()).data();
  double* luhCoarse = DataHeap::getInstance().getData(
      coarseGridCellDescription.getSolution()).data();

  volumeUnknownsRestriction(
      luhCoarse,luhFine,
      levelCoarse,levelFine,
      subcellIndex);

  // TODO(Dominic): What to do in this case with the time step data for anarchic time stepping?
}

void exahype::solvers::ADERDGSolver::finaliseStateUpdates(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) {
  int element =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (element!=exahype::solvers::Solver::NotFound) {
    CellDescription& cellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),element);

    cellDescription.setNewlyCreated(false);
  }
}

////////////////////////////////////////
// CELL-LOCAL
////////////////////////////////////////
void exahype::solvers::ADERDGSolver::validateNoNansInADERDGSolver(
  const CellDescription& cellDescription,
  const std::string& methodTraceOfCaller
) {
  int dataPerCell             = 0;
  int unknownsPerCell         = 0;
  int dataPerCellBoundary     = 0;
  int unknownsPerCellBoundary = 0;

  #if defined(Debug) || defined(Asserts)
  double* luh = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
  double* lduh = DataHeap::getInstance().getData(cellDescription.getUpdate()).data();

  double* lQhbnd = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data();
  double* lFhbnd = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data();

  dataPerCell             = getDataPerCell();
  unknownsPerCell         = getUnknownsPerCell();

  dataPerCellBoundary     = getBndTotalSize();
  unknownsPerCellBoundary = getBndFluxTotalSize();
  #endif

  assertion1(getType()==exahype::solvers::Solver::Type::ADERDG,cellDescription.toString());

  assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()),cellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getUpdate()),cellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor()),cellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()),cellDescription.toString());

  for (int i=0; i<dataPerCell; i++) {
    assertion4(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(luh[i]),
        cellDescription.toString(),toString(),methodTraceOfCaller,i);
  }

  for (int i=0; i<unknownsPerCell; i++) {
    assertion4(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(lduh[i]),
          cellDescription.toString(),toString(),methodTraceOfCaller,i);
  }

  for (int i=0; i<dataPerCellBoundary; i++) {
    assertion4(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(lQhbnd[i]),
          cellDescription.toString(),toString(),methodTraceOfCaller,i);
  }

  for (int i=0; i<unknownsPerCellBoundary; i++) {
    assertion4(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(lFhbnd[i]),
        cellDescription.toString(),toString(),methodTraceOfCaller,i);
  } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
}

bool exahype::solvers::ADERDGSolver::evaluateRefinementCriterionAfterSolutionUpdate(
      const int cellDescriptionsIndex,
      const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Type::Cell) {
    const double* solution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
    RefinementControl refinementControl = refinementCriterion(
                      solution,cellDescription.getOffset()+0.5*cellDescription.getSize(),
                      cellDescription.getSize(),
                      cellDescription.getCorrectorTimeStamp()+cellDescription.getCorrectorTimeStepSize(),
                      cellDescription.getLevel());

    // TODO(Dominic): Set cell description refinement events? Yes or no?

    return (refinementControl==RefinementControl::Refine);
  }

  return false;
}

void exahype::solvers::ADERDGSolver::performPredictionAndVolumeIntegral(
    exahype::records::ADERDGCellDescription& cellDescription,
    double** tempSpaceTimeUnknowns,
    double** tempSpaceTimeFluxUnknowns,
    double*  tempUnknowns,
    double*  tempFluxUnknowns,
    double*  tempStateSizedVector,
    double*  tempPointForceSources) {
  assertion1(cellDescription.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,cellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()),cellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getUpdate()),cellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor()),cellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()),cellDescription.toString());

  assertion2(std::isfinite(cellDescription.getPredictorTimeStepSize()),
             cellDescription.toString(),toString());
  assertion3(cellDescription.getPredictorTimeStepSize()<
             std::numeric_limits<double>::max(),
             cellDescription.toString(),toString(),tarch::parallel::Node::getInstance().getRank());
  assertion2(cellDescription.getPredictorTimeStepSize()>0,
             cellDescription.toString(),toString());

  assertion2(std::isfinite(cellDescription.getPredictorTimeStamp()),
             cellDescription.toString(),toString());
  assertion2(cellDescription.getPredictorTimeStamp()<
             std::numeric_limits<double>::max(),
             cellDescription.toString(),toString());
  assertion2(cellDescription.getPredictorTimeStamp()>=0,
             cellDescription.toString(),toString());

  // persistent fields
  // volume DoF (basisSize**(DIMENSIONS))
  double* luh  = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
  double* lduh = DataHeap::getInstance().getData(cellDescription.getUpdate()).data();
  // face DoF (basisSize**(DIMENSIONS-1))
  double* lQhbnd = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data();
  double* lFhbnd = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data();

  for (int i=0; i<getUnknownsPerCell(); i++) { // cellDescription.getCorrectorTimeStepSize==0.0 is an initial condition
    assertion3(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(luh[i]),cellDescription.toString(),"performPredictionAndVolumeIntegral(...)",i);
  } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.

  if(usePointSource()) { //disable kernel if not needed
      pointSource(cellDescription.getCorrectorTimeStamp() , cellDescription.getCorrectorTimeStepSize(), cellDescription.getOffset()+0.5*cellDescription.getSize(), cellDescription.getSize(), tempPointForceSources); //TODO KD
      // luh, t, dt, cell cell center, cell size, data allocation for forceVect
    }
//TODO JMG move everything to inverseDx and use Peano to get it when Dominic implemente it
#ifdef OPT_KERNELS
  double* dx = &cellDescription.getSize()[0];
#if DIMENSIONS==2
  double inverseDx[2];
#else
  double inverseDx[3];
  inverseDx[2] = 1.0/dx[2];
#endif
  inverseDx[0] = 1.0/dx[0];
  inverseDx[1] = 1.0/dx[1];
  spaceTimePredictor(
      lQhbnd,
      lFhbnd,
      tempSpaceTimeUnknowns,
      tempSpaceTimeFluxUnknowns,
      tempUnknowns,
      tempFluxUnknowns,
      tempStateSizedVector,
      luh,
      &inverseDx[0], //TODO JMG use cellDescription.getInverseSize() when implemented
      cellDescription.getPredictorTimeStepSize(),
      tempPointForceSources);
      
  // TODO(Future Opt.)
  // Volume integral should be performed using the space time
  // flux unknowns. Something equivalent can also be done for
  // the extrpolated fluxes. Here, we can also perform the
  // time averaging on the fly.
  // Remove the tempFluxUnkowns and tempUnknowns.
  volumeIntegral(
      lduh,
      tempFluxUnknowns,
      &inverseDx[0]); //TODO JMG use cellDescription.getInverseSize() when implemented
#else 
  spaceTimePredictor(
      lQhbnd,
      lFhbnd,
      tempSpaceTimeUnknowns,
      tempSpaceTimeFluxUnknowns,
      tempUnknowns,
      tempFluxUnknowns,
      tempStateSizedVector,
      luh,
      cellDescription.getSize(),
      cellDescription.getPredictorTimeStepSize(),
      tempPointForceSources);

  // TODO(Future Opt.)
  // Volume integral should be performed using the space time
  // flux unknowns. Something equivalent can also be done for
  // the extrpolated fluxes. Here, we can also perform the
  // time averaging on the fly.
  // Remove the tempFluxUnkowns and tempUnknowns.
  volumeIntegral(
      lduh,
      tempFluxUnknowns,
      cellDescription.getSize());
#endif

  for (int i=0; i<getTempSpaceTimeUnknownsSize(); i++) { // cellDescription.getCorrectorTimeStepSize==0.0 is an initial condition
    assertion3(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(tempSpaceTimeUnknowns[0][i]),cellDescription.toString(),"performPredictionAndVolumeIntegral(...)",i);
  } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
  for (int i=0; i<getSpaceTimeFluxUnknownsPerCell(); i++) {
    assertion3(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(tempSpaceTimeFluxUnknowns[0][i]), cellDescription.toString(),"performPredictionAndVolumeIntegral",i);
  } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.

  #if defined(Debug) || defined(Asserts)
  if(usePaddedData_nVar()) {
    //TODO JMG add assert ignoring padding
  } else {
    for (int i=0; i<getDataPerCell(); i++) {
    assertion3(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(tempUnknowns[i]),cellDescription.toString(),"performPredictionAndVolumeIntegral(...)",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
    for (int i=0; i<getFluxUnknownsPerCell(); i++) {
      assertion3(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(tempFluxUnknowns[i]),cellDescription.toString(),"performPredictionAndVolumeIntegral(...)",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
  }
  #endif
}

double exahype::solvers::ADERDGSolver::startNewTimeStep(
    const int cellDescriptionsIndex,
    const int element,
    double*   tempEigenvalues) {
  CellDescription& cellDescription =
      exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex)[element];

  if (cellDescription.getType()==exahype::records::ADERDGCellDescription::Cell) {
    assertion1(cellDescription.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,cellDescription.toString());
    const double* luh = exahype::DataHeap::getInstance().getData(cellDescription.getSolution()).data();

    validateNoNansInADERDGSolver(cellDescription,"startNewTimeStep(...)");
//TODO JMG move everything to inverseDx and use Peano to get it when Dominic implemente it
#ifdef OPT_KERNELS
    double* dx = &cellDescription.getSize()[0];
#if DIMENSIONS==2
    double inverseDx[2];
#else
    double inverseDx[3];
    inverseDx[2] = 1.0/dx[2];
#endif
    inverseDx[0] = 1.0/dx[0];
    inverseDx[1] = 1.0/dx[1];
    double admissibleTimeStepSize =
        stableTimeStepSize(luh,tempEigenvalues,&inverseDx[0]); //TODO JMG use cellDescription.getInverseSize() when implemented
#else
    double admissibleTimeStepSize =
        stableTimeStepSize(luh,tempEigenvalues,cellDescription.getSize());
#endif
    assertion2(admissibleTimeStepSize>0,admissibleTimeStepSize,cellDescription.toString());
    assertion3(admissibleTimeStepSize<std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),admissibleTimeStepSize,cellDescription.toString());
    assertion2(std::isfinite(admissibleTimeStepSize),admissibleTimeStepSize,cellDescription.toString());

    // Direct update of the cell description time steps
    // Note that these local quantities might
    // be overwritten again by the synchronisation
    // happening in the next time step.
    // n-1
    cellDescription.setPreviousCorrectorTimeStamp(cellDescription.getCorrectorTimeStamp());
    cellDescription.setPreviousCorrectorTimeStepSize(cellDescription.getCorrectorTimeStepSize());

    // n
    cellDescription.setCorrectorTimeStamp(cellDescription.getPredictorTimeStamp());
    cellDescription.setCorrectorTimeStepSize(cellDescription.getPredictorTimeStepSize());

    // n+1
    cellDescription.setPredictorTimeStamp(cellDescription.getPredictorTimeStamp() + cellDescription.getPredictorTimeStepSize());
    cellDescription.setPredictorTimeStepSize(admissibleTimeStepSize);

    return admissibleTimeStepSize;
  }

  return std::numeric_limits<double>::max();
}

void exahype::solvers::ADERDGSolver::zeroTimeStepSizes(const int cellDescriptionsIndex, const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Cell) {
    cellDescription.setCorrectorTimeStepSize(0.0);
    cellDescription.setPredictorTimeStepSize(0.0);

    cellDescription.setPredictorTimeStamp(
        cellDescription.getCorrectorTimeStamp());
  }
}

void exahype::solvers::ADERDGSolver::reconstructStandardTimeSteppingData(const int cellDescriptionsIndex,int element) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Cell) {
    //  cellDescription.setPreviousCorrectorTimeStepSize(cellDescription.getCorrectorTimeStepSize()); TODO(Dominic): Should not be necessary: Prove by induction
//    logInfo("reconstructStandardTimeSteppingData(...)","cellDescription.getCorrectorTimeStamp()="<<cellDescription.getCorrectorTimeStamp());
//    logInfo("reconstructStandardTimeSteppingData(...)","cellDescription.getCorrectorTimeStepSize()="<<cellDescription.getCorrectorTimeStepSize()); TODO(Dominic): remove

    cellDescription.setPredictorTimeStamp(cellDescription.getCorrectorTimeStamp()+cellDescription.getCorrectorTimeStepSize());
    cellDescription.setCorrectorTimeStamp(cellDescription.getPredictorTimeStamp());
    cellDescription.setCorrectorTimeStepSize(cellDescription.getPredictorTimeStepSize());

    assertionEquals(cellDescription.getCorrectorTimeStamp(),cellDescription.getPredictorTimeStamp());
    assertionEquals(cellDescription.getCorrectorTimeStepSize(),cellDescription.getPredictorTimeStepSize());
  }
}


void exahype::solvers::ADERDGSolver::rollbackToPreviousTimeStep(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  // n+1
  cellDescription.setPredictorTimeStamp(cellDescription.getCorrectorTimeStamp());
  cellDescription.setPredictorTimeStepSize(cellDescription.getCorrectorTimeStepSize());

  // n
  cellDescription.setCorrectorTimeStamp(cellDescription.getPreviousCorrectorTimeStamp());
  cellDescription.setCorrectorTimeStepSize(cellDescription.getPreviousCorrectorTimeStepSize());

  // n-1
  cellDescription.setPreviousCorrectorTimeStamp(std::numeric_limits<double>::max());
  cellDescription.setPreviousCorrectorTimeStepSize(std::numeric_limits<double>::max()); // TODO(Dominic): get rid of the last time level.
}

void exahype::solvers::ADERDGSolver::reconstructStandardTimeSteppingDataAfterRollback(
    const int cellDescriptionsIndex,
    const int element) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  cellDescription.setPredictorTimeStamp(cellDescription.getCorrectorTimeStamp());
  cellDescription.setPredictorTimeStepSize(cellDescription.getCorrectorTimeStepSize());

  cellDescription.setPreviousCorrectorTimeStepSize(std::numeric_limits<double>::max());
}

void exahype::solvers::ADERDGSolver::setInitialConditions(
    const int cellDescriptionsIndex,
    const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  // reset helper variables
  CellDescription& cellDescription  = getCellDescription(cellDescriptionsIndex,element);
  exahype::Cell::resetNeighbourMergeHelperVariables(
          cellDescription,fineGridVertices,fineGridVerticesEnumerator);

  // initial conditions
  if (cellDescription.getType()==exahype::records::ADERDGCellDescription::Cell &&
      cellDescription.getRefinementEvent()==exahype::records::ADERDGCellDescription::None) {
    double* luh = exahype::DataHeap::getInstance().getData(cellDescription.getSolution()).data();

    if (
      useAdjustSolution(
        cellDescription.getOffset()+0.5*cellDescription.getSize(),
        cellDescription.getSize(),
        cellDescription.getCorrectorTimeStamp(),
        cellDescription.getCorrectorTimeStepSize()
      )
      !=AdjustSolutionValue::No
    ) {
        adjustSolution(
            luh,
            cellDescription.getOffset()+0.5*cellDescription.getSize(),
            cellDescription.getSize(),
            cellDescription.getCorrectorTimeStamp(),
            cellDescription.getCorrectorTimeStepSize());
    }

    for (int i=0; i<getUnknownsPerCell(); i++) {
      assertion3(std::isfinite(luh[i]),cellDescription.toString(),"setInitialConditions(...)",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
  }
}

void exahype::solvers::ADERDGSolver::updateSolution(
    const int cellDescriptionsIndex,
    const int element,
    double** tempStateSizedArrays,
    double** tempUnknowns,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  // reset helper variables
  CellDescription& cellDescription  = getCellDescription(cellDescriptionsIndex,element);
  exahype::Cell::resetNeighbourMergeHelperVariables(
        cellDescription,fineGridVertices,fineGridVerticesEnumerator);

  if (cellDescription.getType()==exahype::records::ADERDGCellDescription::Cell &&
      cellDescription.getRefinementEvent()==exahype::records::ADERDGCellDescription::None) {
    double* solution    = DataHeap::getInstance().getData(cellDescription.getPreviousSolution()).data();
    double* newSolution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
    std::copy(newSolution,newSolution+_dofPerCell,solution); // Copy (current solution) in old solution field.

    double* lduh   = exahype::DataHeap::getInstance().getData(cellDescription.getUpdate()).data();
    double* lFhbnd = exahype::DataHeap::getInstance().getData(cellDescription.getFluctuation()).data();

    for (int i=0; i<getDataPerCell(); i++) { // cellDescription.getCorrectorTimeStepSize()==0.0 is an initial condition
      assertion3(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0)  || std::isfinite(solution[i]),cellDescription.toString(),"updateSolution(...)",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.

    for (int i=0; i<getUnknownsPerCell(); i++) {
      assertion3(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0)  || std::isfinite(lduh[i]),cellDescription.toString(),"updateSolution",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.

    for (int i=0; i<getBndFluxTotalSize(); i++) {
      assertion3(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0)  || tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(lFhbnd[i]),cellDescription.toString(),"updateSolution",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
//TODO JMG move everything to inverseDx and use Peano to get it when Dominic implemente it
#ifdef OPT_KERNELS
    double* dx = &cellDescription.getSize()[0];
#if DIMENSIONS==2
    double inverseDx[2];
#else
    double inverseDx[3];
    inverseDx[2] = 1.0/dx[2];
#endif
    inverseDx[0] = 1.0/dx[0];
    inverseDx[1] = 1.0/dx[1];
    surfaceIntegral(lduh,lFhbnd,&inverseDx[0]); //TODO JMG use cellDescription.getInverseSize() when implemented
#else 
    surfaceIntegral(lduh,lFhbnd,cellDescription.getSize());
#endif

    for (int i=0; i<getUnknownsPerCell(); i++) {
      assertion3(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0)  || std::isfinite(lduh[i]),cellDescription.toString(),"updateSolution(...)",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.

    solutionUpdate(newSolution,lduh,cellDescription.getCorrectorTimeStepSize());
    
    if (
      useAdjustSolution(
        cellDescription.getOffset()+0.5*cellDescription.getSize(),
        cellDescription.getSize(),
        cellDescription.getCorrectorTimeStamp()+cellDescription.getCorrectorTimeStepSize(),
        cellDescription.getCorrectorTimeStamp()
      )
      !=AdjustSolutionValue::No
    ) {
      adjustSolution(
          newSolution,
          cellDescription.getOffset()+0.5*cellDescription.getSize(),
          cellDescription.getSize(),
          cellDescription.getCorrectorTimeStamp()+cellDescription.getCorrectorTimeStepSize(),
          cellDescription.getCorrectorTimeStepSize());
    }

    for (int i=0; i<getDataPerCell(); i++) {
      assertion4(std::isfinite(newSolution[i]),cellDescriptionsIndex,cellDescription.toString(),"updateSolution(...)",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
  }
  assertion(cellDescription.getRefinementEvent()==exahype::records::ADERDGCellDescription::None);
}

void exahype::solvers::ADERDGSolver::rollbackSolution(
    const int cellDescriptionsIndex,
    const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  // reset helper variables
  CellDescription& cellDescription  = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==exahype::records::ADERDGCellDescription::Cell &&
      cellDescription.getRefinementEvent()==exahype::records::ADERDGCellDescription::None) {
    swapSolutionAndPreviousSolution(cellDescriptionsIndex,element);
  }
  assertion(cellDescription.getRefinementEvent()==exahype::records::ADERDGCellDescription::None);
}

void exahype::solvers::ADERDGSolver::swapSolutionAndPreviousSolution(
    const int cellDescriptionsIndex,
    const int element) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  assertion(cellDescription.getType()==exahype::records::ADERDGCellDescription::Cell);
  assertion(cellDescription.getRefinementEvent()==exahype::records::ADERDGCellDescription::None);

  // Simply swap the heap indices
  const int previousSolution = cellDescription.getPreviousSolution();
  cellDescription.setPreviousSolution(cellDescription.getSolution());
  cellDescription.setSolution(previousSolution);
}

void exahype::solvers::ADERDGSolver::preProcess(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (
    cellDescription.getType()==CellDescription::Type::Cell
    #ifdef Parallel
    &&
    !cellDescription.getAdjacentToRemoteRank()
    #endif
  ) {
    uncompress(cellDescription);
  }
}

void exahype::solvers::ADERDGSolver::postProcess(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (
      cellDescription.getType()==CellDescription::Type::Cell
      #ifdef Parallel
      &&
      !cellDescription.getAdjacentToRemoteRank()
      #endif
    ) {
    compress(cellDescription);
  }
}

void exahype::solvers::ADERDGSolver::prepareFaceDataOfAncestor(CellDescription& cellDescription) {
  logDebug("prepareFaceDataOfAncestor(...)","cell="<<cellDescription.getOffset()+0.5*cellDescription.getSize() <<
          ", level=" << cellDescription.getLevel());

  assertion1(cellDescription.getType()==CellDescription::Ancestor,cellDescription.toString());
  std::fill_n(DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).begin(),
              getBndTotalSize(), 0.0);
  std::fill_n(DataHeap::getInstance().getData(cellDescription.getFluctuation()).begin(),
              getBndFluxTotalSize(), 0.0);

  #if defined(Debug) || defined(Asserts)
  double* Q = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data();
  double* F = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data();
  #endif

  for(int i=0; i<getBndTotalSize(); ++i) {
    assertion2(tarch::la::equals(Q[i],0.0),i,Q[i]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.


  for(int i=0; i<getBndFluxTotalSize(); ++i) {
    assertion2(tarch::la::equals(F[i],0.0),i,F[i]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.
}

void exahype::solvers::ADERDGSolver::prolongateFaceDataToDescendant(
    CellDescription& cellDescription,
    SubcellPosition& subcellPosition) {
  assertion2(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(
      subcellPosition.parentCellDescriptionsIndex),
      subcellPosition.parentCellDescriptionsIndex,cellDescription.toString());

  CellDescription& cellDescriptionParent = Heap::getInstance().getData(
      subcellPosition.parentCellDescriptionsIndex)[subcellPosition.parentElement];

  assertion(cellDescriptionParent.getSolverNumber() == cellDescription.getSolverNumber());
  assertion(cellDescriptionParent.getType() == exahype::records::ADERDGCellDescription::Cell ||
            cellDescriptionParent.getType() == exahype::records::ADERDGCellDescription::Descendant);

  const int levelFine   = cellDescription.getLevel();
  const int levelCoarse = cellDescriptionParent.getLevel();
  assertion(levelCoarse < levelFine);
  const int levelDelta = levelFine - levelCoarse;

  for (int d = 0; d < DIMENSIONS; ++d) {
    // Check if cell is at "left" or "right" d face of parent
    if (subcellPosition.subcellIndex[d]==0 ||
        subcellPosition.subcellIndex[d]==tarch::la::aPowI(levelDelta,3)-1) {
      const int faceIndex = 2*d + ((subcellPosition.subcellIndex[d]==0) ? 0 : 1); // Do not remove brackets.

      const int numberOfFaceDof = getBndFaceSize();
      const int numberOfFluxDof = getBndFluxSize();

      double* lQhbndFine = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
          (faceIndex * numberOfFaceDof);
      const double* lQhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getExtrapolatedPredictor()).data() +
          (faceIndex * numberOfFaceDof);

      double* lFhbndFine = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
          (faceIndex * numberOfFluxDof);
      const double* lFhbndCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getFluctuation()).data() +
          (faceIndex * numberOfFluxDof);

      faceUnknownsProlongation(lQhbndFine,lFhbndFine,lQhbndCoarse,
                               lFhbndCoarse, levelCoarse, levelFine,
                               exahype::amr::getSubfaceIndex(subcellPosition.subcellIndex,d));
    }
  }
}

void exahype::solvers::ADERDGSolver::prolongateDataAndPrepareDataRestriction(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& cellDescription =
      exahype::solvers::ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Type::Ancestor) {
    prepareFaceDataOfAncestor(cellDescription);
  } else if (cellDescription.getType()==CellDescription::Type::Descendant &&
      isValidCellDescriptionIndex(cellDescription.getParentIndex())) {
    exahype::solvers::Solver::SubcellPosition
    subcellPosition =
        exahype::amr::computeSubcellPositionOfDescendant<CellDescription,Heap,false>(
            cellDescription);

    prolongateFaceDataToDescendant(cellDescription,subcellPosition);
  }
  assertion1(
      cellDescription.getType()!=CellDescription::Type::Descendant ||
      (isValidCellDescriptionIndex(cellDescription.getParentIndex()) ||
          cellDescription.getParentIndex()==multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex),
          cellDescription.toString());
}

void exahype::solvers::ADERDGSolver::restrictData(
                  const int cellDescriptionsIndex,
                  const int element,
                  const int parentCellDescriptionsIndex,
                  const int parentElement,
                  const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  CellDescription& cellDescription =
      getCellDescription(cellDescriptionsIndex,element);
  CellDescription& parentCellDescription =
      getCellDescription(parentCellDescriptionsIndex,parentElement);
  assertion(parentCellDescription.getSolverNumber()==cellDescription.getSolverNumber());
  assertion1(parentCellDescription.getType()==CellDescription::Type::Ancestor,
            parentCellDescription.toString());

  const int levelFine   = cellDescription.getLevel();
  const int levelCoarse = parentCellDescription.getLevel();
  assertion(levelCoarse < levelFine);
  const int levelDelta  = levelFine - levelCoarse;

  for (int d = 0; d < DIMENSIONS; d++) {
    if (subcellIndex[d]==0 ||
        subcellIndex[d]==tarch::la::aPowI(levelDelta,3)-1) {
      const int faceIndex = 2*d + ((subcellIndex[d]==0) ? 0 : 1); // Do not remove brackets.

      logDebug("restrictData(...)","cell=" << cellDescription.getOffset()+0.5*cellDescription.getSize() <<
               ",level=" << cellDescription.getLevel() <<
               ",d=" << d <<
               ",face=" << faceIndex << ",subcellIndex" << subcellIndex.toString() <<
               " to " <<
               " cell="<<parentCellDescription.getOffset()+0.5*parentCellDescription.getSize()<<
               " level="<<parentCellDescription.getLevel());

      #ifdef Parallel
      logDebug("restrictData(...)","forNeighbourComm="<<cellDescription.getHasToHoldDataForNeighbourCommunication() <<
               ",forMasterWorkerComm="<<cellDescription.getHasToHoldDataForMasterWorkerCommunication());
      #endif

      const int numberOfFaceDof = getBndFaceSize();
      const int numberOfFluxDof = getBndFluxSize();

      const double* lQhbndFine = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
          (faceIndex * numberOfFaceDof);
      double* lQhbndCoarse = DataHeap::getInstance().getData(parentCellDescription.getExtrapolatedPredictor()).data() +
          (faceIndex * numberOfFaceDof);

      const double* lFhbndFine = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
          (faceIndex * numberOfFluxDof);
      double* lFhbndCoarse = DataHeap::getInstance().getData(parentCellDescription.getFluctuation()).data() +
          (faceIndex * numberOfFluxDof);

      faceUnknownsRestriction(lQhbndCoarse,lFhbndCoarse,lQhbndFine,lFhbndFine,
                              levelCoarse, levelFine,
                              exahype::amr::getSubfaceIndex(subcellIndex,d));
    }
  }
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::ADERDGSolver::mergeNeighbours(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));
  // !!! In Riemann solve we consider "left" face of "right" cell and
  // "right" face of "left" cell. !!!
  const int normalDirection = tarch::la::equalsReturnIndex(pos1, pos2);
  assertion(normalDirection >= 0 && normalDirection < DIMENSIONS);
  const int faceIndex1 = 2 * normalDirection +
      (pos2(normalDirection) > pos1(normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
  const int faceIndex2 = 2 * normalDirection +
      (pos1(normalDirection) > pos2(normalDirection) ? 1 : 0);   // !!! Be aware of the ">" !!!

  int cellDescriptionsIndexLeft  = cellDescriptionsIndex1;
  int elementLeft                = element1;
  int faceIndexLeft              = faceIndex1;

  int cellDescriptionsIndexRight = cellDescriptionsIndex2;
  int elementRight               = element2;
  int faceIndexRight             = faceIndex2;

  if (pos1(normalDirection) > pos2(normalDirection)) {
    cellDescriptionsIndexLeft  = cellDescriptionsIndex2;
    elementLeft                = element2;
    faceIndexLeft              = faceIndex2;

    cellDescriptionsIndexRight = cellDescriptionsIndex1;
    elementRight               = element1;
    faceIndexRight             = faceIndex1;
  }

  CellDescription& pLeft  = getCellDescription(cellDescriptionsIndexLeft,elementLeft);
  CellDescription& pRight = getCellDescription(cellDescriptionsIndexRight,elementRight);

  peano::datatraversal::TaskSet uncompression(
    [&] () -> void {
      uncompress(pLeft);
    },
    [&] () -> void {
      uncompress(pRight);
    },
    true
  );

  solveRiemannProblemAtInterface(
      pLeft,pRight,faceIndexLeft,faceIndexRight,
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
}

void exahype::solvers::ADERDGSolver::solveRiemannProblemAtInterface(
    CellDescription& pLeft,
    CellDescription& pRight,
    const int faceIndexLeft,
    const int faceIndexRight,
    double**  tempFaceUnknowns,
    double**  tempStateSizedVectors,
    double**  tempStateSizedSquareMatrices) {
  if (pLeft.getType()==CellDescription::Cell ||
      pRight.getType()==CellDescription::Cell) {
    assertion1(DataHeap::getInstance().isValidIndex(pLeft.getExtrapolatedPredictor()),pLeft.toString());
    assertion1(DataHeap::getInstance().isValidIndex(pLeft.getFluctuation()),pLeft.toString());
    assertion1(DataHeap::getInstance().isValidIndex(pRight.getExtrapolatedPredictor()),pRight.toString());
    assertion1(DataHeap::getInstance().isValidIndex(pRight.getFluctuation()),pRight.toString());
    assertion1(holdsFaceData(pLeft.getType()),pLeft.toString());
    assertion1(holdsFaceData(pRight.getType()),pRight.toString());
    assertion1(pLeft.getRefinementEvent()==CellDescription::None,pLeft.toString());
    assertion1(pRight.getRefinementEvent()==CellDescription::None,pRight.toString());
    assertionEquals4(pLeft.getRiemannSolvePerformed(faceIndexLeft),pRight.getRiemannSolvePerformed(faceIndexRight),faceIndexLeft,faceIndexRight,pLeft.toString(),pRight.toString());
    assertion4(std::abs(faceIndexLeft-faceIndexRight)==1,faceIndexLeft,faceIndexRight,pLeft.toString(),pRight.toString());

    const int dataPerFace = getBndFaceSize();
    const int dofPerFace  = getBndFluxSize();

    double* QL = DataHeap::getInstance() .getData(pLeft.getExtrapolatedPredictor()).data() + /// !!! Be aware of the dataPerFace, Left, Right
        (faceIndexLeft * dataPerFace);
    double* QR = DataHeap::getInstance().getData(pRight.getExtrapolatedPredictor()).data() +
        (faceIndexRight * dataPerFace);

    double* FL = DataHeap::getInstance().getData(pLeft.getFluctuation()).data() + /// !!! Be aware of the dofPerFace, Left, Right
        (faceIndexLeft * dofPerFace);
    double* FR = DataHeap::getInstance().getData(pRight.getFluctuation()).data() +
        (faceIndexRight * dofPerFace);

    // todo Time step must be interpolated in local time stepping case
    // both time step sizes are the same, so the min has no effect here.
    assertion1(faceIndexLeft>=0,faceIndexLeft);
    assertion1(faceIndexRight>=0,faceIndexRight);
    assertion1(faceIndexLeft<DIMENSIONS_TIMES_TWO,faceIndexLeft);
    assertion1(faceIndexRight<DIMENSIONS_TIMES_TWO,faceIndexRight);
    assertion1(faceIndexRight%2==0,faceIndexRight);
    const int normalDirection = (faceIndexRight - (faceIndexRight %2))/2;
    assertion3(normalDirection==(faceIndexLeft - (faceIndexLeft %2))/2,normalDirection,faceIndexLeft,faceIndexRight);
    assertion3(normalDirection<DIMENSIONS,normalDirection,faceIndexLeft,faceIndexRight);

    // Synchronise time stepping.
    synchroniseTimeStepping(pLeft);
    synchroniseTimeStepping(pRight);

    assertion3(std::isfinite(pLeft.getCorrectorTimeStepSize()),pLeft.toString(),faceIndexLeft,normalDirection);
    assertion3(std::isfinite(pRight.getCorrectorTimeStepSize()),pRight.toString(),faceIndexRight,normalDirection);
    assertion3(pLeft.getCorrectorTimeStepSize()>=0.0,pLeft.toString(),faceIndexLeft,normalDirection);
    assertion3(pRight.getCorrectorTimeStepSize()>=0.0,pRight.toString(),faceIndexRight,normalDirection);

    for(int i=0; i<dataPerFace; ++i) {
      assertion5(tarch::la::equals(pLeft.getCorrectorTimeStepSize(),0.0) || std::isfinite(QL[i]),pLeft.toString(),faceIndexLeft,normalDirection,i,QL[i]);
      assertion5(tarch::la::equals(pRight.getCorrectorTimeStepSize(),0.0) || std::isfinite(QR[i]),pRight.toString(),faceIndexRight,normalDirection,i,QR[i]);
    }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

    for(int i=0; i<dofPerFace; ++i) {
      assertion5(tarch::la::equals(pLeft.getCorrectorTimeStepSize(),0.0) || std::isfinite(FL[i]),pLeft.toString(),faceIndexLeft,normalDirection,i,FL[i]);
      assertion5(tarch::la::equals(pRight.getCorrectorTimeStepSize(),0.0) || std::isfinite(FR[i]),pRight.toString(),faceIndexRight,normalDirection,i,FR[i]);
    }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

    riemannSolver(
        FL,FR,QL,QR,
        tempFaceUnknowns[0],tempStateSizedVectors,tempStateSizedSquareMatrices,
        std::min(pLeft.getCorrectorTimeStepSize(),
            pRight.getCorrectorTimeStepSize()),
            normalDirection);

    for(int i=0; i<dofPerFace; ++i) {
      assertion8(tarch::la::equals(pLeft.getCorrectorTimeStepSize(),0.0) || (std::isfinite(FL[i]) && std::isfinite(FR[i])),
                 pLeft.toString(),faceIndexLeft,pRight.toString(),faceIndexRight,normalDirection,i,FL[i],FR[i]);
    }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

    assertion1(DataHeap::getInstance().isValidIndex(pLeft.getExtrapolatedPredictor()),pLeft.toString());
    assertion1(DataHeap::getInstance().isValidIndex(pLeft.getFluctuation()),pLeft.toString());
    assertion1(DataHeap::getInstance().isValidIndex(pRight.getExtrapolatedPredictor()),pRight.toString());
    assertion1(DataHeap::getInstance().isValidIndex(pRight.getFluctuation()),pRight.toString());
  }
}

void exahype::solvers::ADERDGSolver::mergeWithBoundaryData(
    const int                                 cellDescriptionsIndex,
    const int                                 element,
    const tarch::la::Vector<DIMENSIONS, int>& posCell,
    const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  if (tarch::la::countEqualEntries(posCell,posBoundary)!=(DIMENSIONS-1)) {
    return; // We only consider faces; no corners.
  }

  CellDescription& cellDescription =
      getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Cell) {
    // !!! Left face of right cell.
    const int normalOfExchangedFace = tarch::la::equalsReturnIndex(posCell, posBoundary);
    assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
    const int faceIndex = 2 * normalOfExchangedFace +
        (posCell(normalOfExchangedFace) < posBoundary(normalOfExchangedFace) ? 1 : 0);

    uncompress(cellDescription);

    applyBoundaryConditions(
        cellDescription,faceIndex,
        tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
  }
}

void exahype::solvers::ADERDGSolver::applyBoundaryConditions(
    CellDescription& p,
    const int faceIndex,
    double**  tempFaceUnknowns,
    double**  tempStateSizedVectors,
    double**  tempStateSizedSquareMatrices) {
  assertion1(p.getType()==CellDescription::Cell,p.toString());
  assertion1(p.getRefinementEvent()==CellDescription::None,p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getExtrapolatedPredictor()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getFluctuation()),p.toString());

  const int dataPerFace = getBndFaceSize();
  const int dofPerFace  = getBndFluxSize();

  double* QIn = DataHeap::getInstance().getData(p.getExtrapolatedPredictor()).data() +
      (faceIndex * dataPerFace);
  double* FIn = DataHeap::getInstance().getData(p.getFluctuation()).data() +
      (faceIndex * dofPerFace);

  const int normalDirection = (faceIndex - faceIndex % 2)/2;
  assertion2(normalDirection<DIMENSIONS,faceIndex,normalDirection);
  
  for(int i=0; i<dataPerFace; ++i) {
    assertion5(std::isfinite(QIn[i]), p.toString(),
        faceIndex, normalDirection, i, QIn[i]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  for(int i=0; i<dofPerFace; ++i) {
    assertion5(std::isfinite(FIn[i]), p.toString(),
        faceIndex, normalDirection, i, FIn[i]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  // Synchronise time stepping.
  synchroniseTimeStepping(p);

  double* QOut = tempFaceUnknowns[1];
  double* FOut  = tempFaceUnknowns[2];

  // TODO(Dominic): Hand in space-time volume data. Time integrate it afterwards
  boundaryConditions(FOut,QOut,
      FIn,QIn,
      p.getOffset() + 0.5*p.getSize(), // centre
      p.getSize(),
      p.getCorrectorTimeStamp(),
      p.getCorrectorTimeStepSize(),
      faceIndex,
      normalDirection);

  assertion4(std::isfinite(p.getCorrectorTimeStamp()),p.toString(),faceIndex,normalDirection,p.getCorrectorTimeStamp());
  assertion4(std::isfinite(p.getCorrectorTimeStepSize()),p.toString(),faceIndex,normalDirection,p.getCorrectorTimeStepSize());
  assertion4(p.getCorrectorTimeStepSize()>=0.0, p.toString(),faceIndex, normalDirection,p.getCorrectorTimeStepSize());
  for(int i=0; i<dataPerFace; ++i) {
    assertion5(tarch::la::equals(p.getCorrectorTimeStepSize(),0.0) || std::isfinite(QIn[i]),p.toString(),faceIndex,normalDirection,i,QIn[i]);
    assertion5(tarch::la::equals(p.getCorrectorTimeStepSize(),0.0) || std::isfinite(QOut[i]),p.toString(),faceIndex,normalDirection,i,QOut[i]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.
  for(int i=0; i<dofPerFace; ++i) {
    assertion5(tarch::la::equals(p.getCorrectorTimeStepSize(),0.0) || std::isfinite(FIn[i]),p.toString(),faceIndex,normalDirection,i,FIn[i]);
    assertion5(tarch::la::equals(p.getCorrectorTimeStepSize(),0.0) || std::isfinite(FOut[i]),p.toString(),faceIndex,normalDirection,i,FOut[i]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  // @todo(Dominic): Add to docu why we need this. Left or right inpclut
  if (faceIndex % 2 == 0) {
    riemannSolver(FOut, FIn, QOut, QIn,
        tempFaceUnknowns[0],tempStateSizedVectors,tempStateSizedSquareMatrices,
        p.getCorrectorTimeStepSize(),
        normalDirection);
  } else {
    riemannSolver(FIn, FOut, QIn, QOut,
        tempFaceUnknowns[0],tempStateSizedVectors,tempStateSizedSquareMatrices,
        p.getCorrectorTimeStepSize(),
        normalDirection);
  }

  for(int i=0; i<dofPerFace; ++i) {
    assertion5(tarch::la::equals(p.getCorrectorTimeStepSize(),0.0) || std::isfinite(FIn[i]),p.toString(),faceIndex,normalDirection,i,FIn[i]);
    assertion5(tarch::la::equals(p.getCorrectorTimeStepSize(),0.0) || std::isfinite(FOut[i]),p.toString(),faceIndex,normalDirection,i,FOut[i]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.
}

void exahype::solvers::ADERDGSolver::prepareNextNeighbourMerging(
    const int cellDescriptionsIndex,const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  exahype::Cell::resetNeighbourMergeHelperVariables(
      cellDescription,fineGridVertices,fineGridVerticesEnumerator);
}

#ifdef Parallel
const int exahype::solvers::ADERDGSolver::DataMessagesPerNeighbourCommunication    = 2;
const int exahype::solvers::ADERDGSolver::DataMessagesPerForkOrJoinCommunication   = 1;
const int exahype::solvers::ADERDGSolver::DataMessagesPerMasterWorkerCommunication = 2;

/**
 * After the forking the master's cell descriptions
 * are not accessed by enterCell(...) on the master
 * anymore. However we still use them as buffer for saving data.
 */
void exahype::solvers::ADERDGSolver::sendCellDescriptions(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),
      cellDescriptionsIndex);

  for (CellDescription& cellDescription : Heap::getInstance().getData(cellDescriptionsIndex)) {
    if (cellDescription.getType()==CellDescription::Type::EmptyAncestor
        || cellDescription.getType()==CellDescription::Type::Ancestor) {
      Solver::SubcellPosition subcellPosition =
          exahype::amr::computeSubcellPositionOfCellOrAncestorOrEmptyAncestor
          <CellDescription,Heap>(cellDescription);

      if (subcellPosition.parentElement!=NotFound) {
        cellDescription.setType(CellDescription::Type::Ancestor);
        cellDescription.setHasToHoldDataForMasterWorkerCommunication(true);

        auto* solver = exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()];

        switch (solver->getType()) {
          case exahype::solvers::Solver::Type::ADERDG:
            static_cast<exahype::solvers::ADERDGSolver*>(solver)->ensureNecessaryMemoryIsAllocated(cellDescription);
            break;
          case exahype::solvers::Solver::Type::LimitingADERDG:
            static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->ensureNecessaryMemoryIsAllocated(cellDescription);
            break;
          case exahype::solvers::Solver::Type::FiniteVolumes:
            assertionMsg(false,"Solver type not supported!");
            break;
        }
      }
    } else if (cellDescription.getType()==CellDescription::Type::EmptyDescendant
        || cellDescription.getType()==CellDescription::Type::Descendant) {
      cellDescription.setType(CellDescription::Type::Descendant);
      cellDescription.setHasToHoldDataForMasterWorkerCommunication(true);

      auto* solver = exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()];

      switch (solver->getType()) {
        case exahype::solvers::Solver::Type::ADERDG:
          static_cast<exahype::solvers::ADERDGSolver*>(solver)->ensureNecessaryMemoryIsAllocated(cellDescription);
          break;
        case exahype::solvers::Solver::Type::LimitingADERDG:
          static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->ensureNecessaryMemoryIsAllocated(cellDescription);
          break;
        case exahype::solvers::Solver::Type::FiniteVolumes:
          assertionMsg(false,"Solver type not supported!");
          break;
      }
    }
  }

  Heap::getInstance().sendData(cellDescriptionsIndex,
                               toRank,x,level,messageType);
}

void exahype::solvers::ADERDGSolver::sendEmptyCellDescriptions(
    const int                                     toRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::HeapEntries emptyMessage(0);
  Heap::getInstance().sendData(emptyMessage,
      toRank,x,level,messageType);
}

void exahype::solvers::ADERDGSolver::mergeCellDescriptionsWithRemoteData(
    const int                                    fromRank,
    exahype::Cell&                               localCell,
    const peano::heap::MessageType&              messageType,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  waitUntilAllBackgroundTasksHaveTerminated();
  tarch::multicore::Lock lock(_heapSemaphore);

  int receivedCellDescriptionsIndex =
      Heap::getInstance().createData(0,exahype::solvers::RegisteredSolvers.size());
  Heap::getInstance().receiveData(
      receivedCellDescriptionsIndex,fromRank,x,level,messageType);

  if (!Heap::getInstance().getData(receivedCellDescriptionsIndex).empty()) {
    // TODO(Dominic): We reset the parent and heap indices of a received cell
    // to -1 and RemoteAdjacencyIndex, respectively.
    // If we receive parent and children cells during a fork event,
    //
    // We use the information of a parentIndex of a fine grid cell description
    // set to RemoteAdjacencyIndex to update the parent index with
    // the index of the coarse grid cell description in enterCell(..).
    resetDataHeapIndices(receivedCellDescriptionsIndex,
                         multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex);

    localCell.setupMetaData();
    assertion1(Heap::getInstance().isValidIndex(localCell.getCellDescriptionsIndex()),
               localCell.getCellDescriptionsIndex());
    Heap::getInstance().getData(localCell.getCellDescriptionsIndex()).reserve(
        std::max(Heap::getInstance().getData(localCell.getCellDescriptionsIndex()).size(),
                 Heap::getInstance().getData(receivedCellDescriptionsIndex).size()));

    logDebug("mergeCellDescriptionsWithRemoteData(...)","Received " <<
            Heap::getInstance().getData(receivedCellDescriptionsIndex).size() <<
            " cell descriptions for cell ("
            "offset="<< x.toString() <<
            "level="<< level << ")");

    for (auto& pReceived : Heap::getInstance().getData(receivedCellDescriptionsIndex)) {
      logDebug("mergeCellDescriptionsWithRemoteData(...)","Received " <<
              " cell description for cell ("
              "offset="<< x.toString() <<
              ",level="<< level <<
              ",isRoot="<< localCell.isRoot() <<
              ",isAssignedToRemoteRank="<< localCell.isAssignedToRemoteRank() <<
              ") with type="<< pReceived.getType());

      bool found = false;
      for (auto& pLocal : Heap::getInstance().getData(localCell.getCellDescriptionsIndex())) {
        if (pReceived.getSolverNumber()==pLocal.getSolverNumber()) {
          found = true;

          assertion(pReceived.getType()==pLocal.getType());
          if (pLocal.getType()==CellDescription::Type::Cell ||
              pLocal.getType()==CellDescription::Type::EmptyAncestor ||
              pLocal.getType()==CellDescription::Type::Ancestor ||
              pLocal.getType()==CellDescription::Type::Descendant
          ) {
            assertionNumericalEquals2(pLocal.getCorrectorTimeStamp(),pReceived.getCorrectorTimeStamp(),
                                      pLocal.toString(),pReceived.toString());
            assertionNumericalEquals2(pLocal.getCorrectorTimeStepSize(),pReceived.getCorrectorTimeStepSize(),
                                      pLocal.toString(),pReceived.toString());
            assertionNumericalEquals2(pLocal.getPredictorTimeStamp(),pReceived.getPredictorTimeStamp(),
                                      pLocal.toString(),pReceived.toString());
            assertionNumericalEquals2(pLocal.getPredictorTimeStepSize(),pReceived.getPredictorTimeStepSize(),
                                      pLocal.toString(),pReceived.toString());
          }
        }
      }

      if (!found) {
        auto* solver = exahype::solvers::RegisteredSolvers[pReceived.getSolverNumber()];

        switch (solver->getType()) {
          case exahype::solvers::Solver::Type::ADERDG:
            static_cast<exahype::solvers::ADERDGSolver*>(solver)->ensureNecessaryMemoryIsAllocated(pReceived);
            break;
          case exahype::solvers::Solver::Type::LimitingADERDG:
            static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->ensureNecessaryMemoryIsAllocated(pReceived);
            break;
          case exahype::solvers::Solver::Type::FiniteVolumes:
            assertionMsg(false,"Solver type not supported!");
            break;
        }


        Heap::getInstance().getData(localCell.getCellDescriptionsIndex()).
            push_back(pReceived);
      }
    }
  }
}

void exahype::solvers::ADERDGSolver::resetDataHeapIndices(
    const int cellDescriptionsIndex,
    const int parentIndex) {
  for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
    p.setParentIndex(parentIndex);

    // Default field data indices
    p.setSolution(-1);
    p.setUpdate(-1);
    p.setExtrapolatedPredictor(-1);
    p.setFluctuation(-1);

    // Limiter meta data (oscillations identificator)
    p.setSolutionMin(-1);
    p.setSolutionMax(-1);
  }
}

/**
 * Drop cell descriptions received from \p fromRank.
 */
void exahype::solvers::ADERDGSolver::dropCellDescriptions(
    const int                                     fromRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::getInstance().receiveData(fromRank,x,level,messageType);
}

///////////////////////////////////
// FORK OR JOIN
///////////////////////////////////

void exahype::solvers::ADERDGSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  CellDescription& p = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  if (p.getType()==CellDescription::Cell) {
    double* solution = DataHeap::getInstance().getData(p.getSolution()).data();

    logDebug("sendDataToWorkerOrMasterDueToForkOrJoin(...)","solution of solver " << p.getSolverNumber() << " sent to rank "<<toRank<<
             ", cell: "<< x << ", level: " << level);

    DataHeap::getInstance().sendData(
        solution, getUnknownsPerCell(), toRank, x, level,
        peano::heap::MessageType::ForkOrJoinCommunication);
  }
}


void exahype::solvers::ADERDGSolver::sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerForkOrJoinCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, toRank, x, level,
        peano::heap::MessageType::ForkOrJoinCommunication);
}


void exahype::solvers::ADERDGSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  auto& p = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex)[element];
  assertion4(tarch::la::equals(x,p.getOffset()+0.5*p.getSize()),x,p.getOffset()+0.5*p.getSize(),level,p.getLevel());
  assertion2(p.getLevel()==level,p.getLevel(),level);

  if (p.getType()==CellDescription::Cell) {
    logDebug("mergeWithRemoteDataDueToForkOrJoin(...)","[solution] receive from rank "<<fromRank<<
             ", cell: "<< x << ", level: " << level);

    DataHeap::getInstance().getData(p.getSolution()).clear();
    DataHeap::getInstance().receiveData(
        p.getSolution(),fromRank,x,level,
        peano::heap::MessageType::ForkOrJoinCommunication);
  }
}

void exahype::solvers::ADERDGSolver::dropWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  for(int receives=0; receives<DataMessagesPerForkOrJoinCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        fromRank, x, level,
        peano::heap::MessageType::ForkOrJoinCommunication);
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::ADERDGSolver::mergeWithNeighbourMetadata(
    const exahype::MetadataHeap::HeapEntries& neighbourMetadata,
    const tarch::la::Vector<DIMENSIONS, int>& src,
    const tarch::la::Vector<DIMENSIONS, int>& dest,
    const int                                 cellDescriptionsIndex,
    const int                                 element) {
  if (tarch::la::countEqualEntries(src,dest)!=DIMENSIONS-1) { // only consider faces
    return;
  }

  CellDescription& p = getCellDescription(cellDescriptionsIndex,element);

  CellDescription::Type neighbourType =
      static_cast<CellDescription::Type>(neighbourMetadata[exahype::MetadataCellType].getU());
  switch(p.getType()) {
    case CellDescription::Cell:
      if (p.getRefinementEvent()==CellDescription::None &&
          (neighbourType==CellDescription::Ancestor ||
              neighbourType==CellDescription::EmptyAncestor)) {
        p.setRefinementEvent(CellDescription::AugmentingRequested);
      }
      break;
    case CellDescription::Descendant:
    case CellDescription::EmptyDescendant:
      // TODO(Dominic): Add to docu what we do here.
      if (neighbourType==CellDescription::Cell) {
        p.setHasToHoldDataForNeighbourCommunication(true);
      }

      // 2. Request further augmentation if necessary (this might get reset if the traversal
      // is able to descend and finds existing descendants).
      if (p.getRefinementEvent()==CellDescription::None &&
          (neighbourType==CellDescription::Ancestor ||
              neighbourType==CellDescription::EmptyAncestor)) {
        p.setRefinementEvent(CellDescription::AugmentingRequested);
      }
      break;
    case CellDescription::Ancestor:
    case CellDescription::EmptyAncestor:
      // TODO(Dominic): Add to docu what we do here.
      if (neighbourType==CellDescription::Cell) {
        p.setHasToHoldDataForNeighbourCommunication(true);
      }
      break;
    default:
      assertionMsg(false,"Should never be entered in static AMR scenarios!");
      break;
  }
}


void exahype::solvers::ADERDGSolver::sendDataToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertion( tarch::la::countEqualEntries(src,dest)==(DIMENSIONS-1) );
/*
  if (tarch::la::countEqualEntries(src,dest)!=(DIMENSIONS-1)) {
    return; // We only consider faces; no corners.
  }
*/

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!

  CellDescription& cellDescription = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  if (holdsFaceData(cellDescription.getType())) {
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor()));
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()));

    const int numberOfFaceDof = getBndFaceSize();
    const double* lQhbnd = DataHeap::getInstance().getData(
        cellDescription.getExtrapolatedPredictor()).data() +
        (faceIndex * numberOfFaceDof);

    const int numberOfFluxDof = getBndFluxSize();
    const double* lFhbnd = DataHeap::getInstance().getData(
        cellDescription.getFluctuation()).data() +
        (faceIndex * numberOfFluxDof);

    logDebug(
        "sendDataToNeighbour(...)",
        "send "<<DataMessagesPerNeighbourCommunication<<" arrays to rank " <<
        toRank << " for cell="<<cellDescription.getOffset()<< " and face=" << faceIndex << " from vertex x=" << x << ", level=" << level <<
        ", src type=" << multiscalelinkedcell::indexToString(cellDescriptionsIndex) <<
        ", src=" << src << ", dest=" << dest <<
        ", counter=" << cellDescription.getFaceDataExchangeCounter(faceIndex)
    );

    // TODO(Dominic): Send predictor time step data separately
//    // We append all the max values to the min values.
//    // We further append the predictor time stamp and
//    // and the predictor time step size
//    std::vector<double> sentMinMax( 2*getNumberOfVariables() + 2 );
//    for (int i=0; i<getNumberOfVariables(); i++) {
//      sentMinMax[i]                        = DataHeap::getInstance().getData( cellDescription.getSolutionMin() )[faceIndex*getNumberOfVariables()+i];
//      sentMinMax[i+getNumberOfVariables()] = DataHeap::getInstance().getData( cellDescription.getSolutionMax() )[faceIndex*getNumberOfVariables()+i];
//    }
//    sentMinMax[2*getNumberOfVariables()] = cellDescription.getPredictorTimeStamp();
//    sentMinMax[2*getNumberOfVariables()+1] = cellDescription.getPredictorTimeStepSize();
//    assertionEquals(sentMinMax.size(),2*static_cast<unsigned int>(getNumberOfVariables())+2);

    // Send order: lQhbnd,lFhbnd
    // Receive order: lFhbnd,lQhbnd
    DataHeap::getInstance().sendData(
        lQhbnd, numberOfFaceDof, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
    DataHeap::getInstance().sendData(
        lFhbnd, numberOfFluxDof, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
    // TODO(Dominic): If anarchic time stepping send the time step over too.
  } else {
    std::vector<double> emptyArray(0,0);

    for(int sends=0; sends<DataMessagesPerNeighbourCommunication; ++sends) {
      DataHeap::getInstance().sendData(
          emptyArray, toRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
    }
  }
}

void exahype::solvers::ADERDGSolver::sendEmptyDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerNeighbourCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
}

void exahype::solvers::ADERDGSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const MetadataHeap::HeapEntries&             neighbourMetadata,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    double**                                     tempFaceUnknowns,
    double**                                     tempStateSizedVectors,
    double**                                     tempStateSizedSquareMatrices,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  if (tarch::la::countEqualEntries(src,dest)!=(DIMENSIONS-1)) {
    return; // We only consider faces; no corners.
  }

  waitUntilAllBackgroundTasksHaveTerminated();
  // @todo Das kommt mir sehr seltsam vor, weil unten ja nochmal ein Lock kommt.
  //tarch::multicore::Lock lock(_heapSemaphore);

  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  CellDescription::Type neighbourType =
      static_cast<CellDescription::Type>(neighbourMetadata[exahype::MetadataCellType].getU());

  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) > dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

  // TODO(Dominic): Add to docu: We only perform a Riemann solve if a Cell is involved.
  // Solving Riemann problems at a Ancestor Ancestor boundary might lead to problems
  // if one Ancestor is just used for restriction.
  if(neighbourType==CellDescription::Type::Cell || cellDescription.getType()==CellDescription::Type::Cell){
    tarch::multicore::Lock lock(_heapSemaphore);

    assertion1(holdsFaceData(neighbourType),neighbourType);
    assertion1(holdsFaceData(cellDescription.getType()),cellDescription.toString());

    assertion4(!cellDescription.getRiemannSolvePerformed(faceIndex),
        faceIndex,cellDescriptionsIndex,cellDescription.getOffset().toString(),cellDescription.getLevel());
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor()));
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()));

    logDebug(
        "mergeWithNeighbourData(...)", "receive "<<DataMessagesPerNeighbourCommunication<<" arrays from rank " <<
        fromRank << " for vertex x=" << x << ", level=" << level <<
        ", src type=" << cellDescription.getType() <<
        ", src=" << src << ", dest=" << dest <<
        ", counter=" << cellDescription.getFaceDataExchangeCounter(faceIndex)
    );

    const int dataPerFace = getBndFaceSize();
    const int dofPerFace  = getBndFluxSize();
    int receivedlQhbndIndex   = DataHeap::getInstance().createData(dataPerFace, dataPerFace);
    int receivedlFhbndIndex   = DataHeap::getInstance().createData(dofPerFace, dofPerFace);
    assertion(!DataHeap::getInstance().getData(receivedlQhbndIndex).empty());
    assertion(!DataHeap::getInstance().getData(receivedlFhbndIndex).empty());


    // Send order: lQhbnd,lFhbnd
    // Receive order: lFhbnd,lQhbnd // TODO change to double variant
    DataHeap::getInstance().receiveData(
        DataHeap::getInstance().getData(receivedlFhbndIndex).data(),dataPerFace,
        fromRank, x, level,peano::heap::MessageType::NeighbourCommunication);
    DataHeap::getInstance().receiveData(
        DataHeap::getInstance().getData(receivedlQhbndIndex).data(),dataPerFace,
        fromRank, x, level,peano::heap::MessageType::NeighbourCommunication);
    logDebug(
        "mergeWithNeighbourData(...)", "[pre] solve Riemann problem with received data." <<
        " cellDescription=" << cellDescription.toString() <<
        ",faceIndexForCell=" << faceIndex <<
        ",normalOfExchangedFac=" << normalOfExchangedFace <<
        ",x=" << x.toString() << ", level=" << level <<
        ", counter=" << cellDescription.getFaceDataExchangeCounter(faceIndex)
    );

    // TODO collect neighbour time step data here out of receivedMinMax and
    // use it for solveRiemannProblemAtInterface.

    solveRiemannProblemAtInterface(
        cellDescription,
        faceIndex,
        receivedlQhbndIndex,
        receivedlFhbndIndex,
        tempFaceUnknowns,
        tempStateSizedVectors,
        tempStateSizedSquareMatrices);

    // TODO(Dominic): If anarchic time stepping, receive the time step too.

    DataHeap::getInstance().deleteData(receivedlQhbndIndex,true);
    DataHeap::getInstance().deleteData(receivedlFhbndIndex,true);
  } else  {
    logDebug(
        "mergeWithNeighbourData(...)", "drop three arrays from rank " <<
        fromRank << " for vertex x=" << x << ", level=" << level <<
        ", src type=" << multiscalelinkedcell::indexToString(cellDescriptionsIndex) <<
        ", src=" << src << ", dest=" << dest <<
        ", counter=" << cellDescription.getFaceDataExchangeCounter(faceIndex)
    );

    dropNeighbourData(fromRank,src,dest,x,level);
  }
}

void exahype::solvers::ADERDGSolver::solveRiemannProblemAtInterface(
    records::ADERDGCellDescription& cellDescription,
    const int faceIndex,
    const int indexOfQValues,
    const int indexOfFValues,
    double**  tempFaceUnknowns,
    double**  tempStateSizedVectors,
    double**  tempStateSizedSquareMatrices) {
  cellDescription.setRiemannSolvePerformed(faceIndex, true);

  const int dataPerFace = getBndFaceSize();
  const int dofPerFace  = getBndFluxSize();

  logDebug("solveRiemannProblemAtInterface(...)",
      "cell-description=" << cellDescription.toString());

  double* QL = 0;
  double* QR = 0;
  double* FL = 0;
  double* FR = 0;

  assertionEquals(DataHeap::getInstance().getData(indexOfQValues).size(),
      static_cast<unsigned int>(dataPerFace));
  assertionEquals(DataHeap::getInstance().getData(indexOfFValues).size(),
      static_cast<unsigned int>(dofPerFace));

  // @todo Doku im Header warum wir das hier brauchen,
  if (faceIndex % 2 == 0) {
    QL = DataHeap::getInstance().getData(indexOfQValues).data();
    QR = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
        (faceIndex * dataPerFace);
    FL = DataHeap::getInstance().getData(indexOfFValues).data();
    FR = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
        (faceIndex * dofPerFace);
  } else {
    QR = DataHeap::getInstance().getData(indexOfQValues).data();
    QL = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
        (faceIndex * dataPerFace);
    FR = DataHeap::getInstance().getData(indexOfFValues).data();
    FL = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
        (faceIndex * dofPerFace);
  }

  // Synchronise time stepping.
  synchroniseTimeStepping(cellDescription);

  const int normalDirection = (faceIndex - faceIndex%2)/2; // faceIndex=2*normalNonZero+f, f=0,1
  assertion2(normalDirection<DIMENSIONS,faceIndex,normalDirection);
  riemannSolver(FL, FR, QL, QR,
      tempFaceUnknowns[0],tempStateSizedVectors,tempStateSizedSquareMatrices,
      cellDescription.getCorrectorTimeStepSize(),
      normalDirection);

  for (int ii = 0; ii<dataPerFace; ii++) {
    assertion10(std::isfinite(QR[ii]), cellDescription.toString(),
        faceIndex, normalDirection, indexOfQValues, indexOfFValues,
        ii, QR[ii], QL[ii], FR[ii], FL[ii]);
    assertion10(std::isfinite(QL[ii]), cellDescription.toString(),
        faceIndex, normalDirection, indexOfQValues, indexOfFValues,
        ii, QR[ii], QL[ii], FR[ii], FL[ii]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  for (int ii = 0; ii<dofPerFace; ii++) {
    assertion10(std::isfinite(FR[ii]), cellDescription.toString(),
        faceIndex, normalDirection, indexOfQValues, indexOfFValues,
        ii, QR[ii], QL[ii], FR[ii], FL[ii]);
    assertion10(std::isfinite(FL[ii]), cellDescription.toString(),
        faceIndex, normalDirection, indexOfQValues, indexOfFValues,
        ii, QR[ii], QL[ii], FR[ii], FL[ii]);
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.
}

void exahype::solvers::ADERDGSolver::dropNeighbourData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  for(int receives=0; receives<DataMessagesPerNeighbourCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////

/*
 * At the time of sending data to the master,
 * we have already performed a time step update locally
 * on the rank. We thus need to communicate the
 * current min predictor time step size to the master.
 * The next min predictor time step size is
 * already reset locally to the maximum double value.
 *
 * However on the master's side, we need to
 * merge the received time step size with
 * the next min predictor time step size since
 * the master has not yet performed a time step update
 * (i.e. called TimeStepSizeComputation::endIteration()).
 */
void exahype::solvers::ADERDGSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level){
  std::vector<double> timeStepDataToReduce(0,4);
  timeStepDataToReduce.push_back(_minPredictorTimeStepSize);
  timeStepDataToReduce.push_back(_gridUpdateRequested ? 1.0 : -1.0); // TODO(Dominic): ugly
  timeStepDataToReduce.push_back(_minCellSize);
  timeStepDataToReduce.push_back(_maxCellSize);

  assertion1(timeStepDataToReduce.size()==4,timeStepDataToReduce.size());
  assertion1(std::isfinite(timeStepDataToReduce[0]),timeStepDataToReduce[0]);
  if (_timeStepping==TimeStepping::Global) {
    assertionNumericalEquals1(_minNextPredictorTimeStepSize,std::numeric_limits<double>::max(),
                                tarch::parallel::Node::getInstance().getRank());
  }

  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToMaster(...)","Sending time step data: " <<
             "data[0]=" << timeStepDataToReduce[0] <<
             ",data[1]=" << timeStepDataToReduce[1] <<
             ",data[2]=" << timeStepDataToReduce[2] <<
             ",data[3]=" << timeStepDataToReduce[3]);
  }

  DataHeap::getInstance().sendData(
      timeStepDataToReduce.data(), timeStepDataToReduce.size(),
      masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

/**
 * At the time of the merging,
 * the workers and the master have already performed
 * at local update of the next predictor time step size
 * and of the predictor time stamp.
 * We thus need to minimise over both quantities.
 */
void exahype::solvers::ADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  std::vector<double> receivedTimeStepData(4);

  if (true || tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","Receiving time step data [pre] from rank " << workerRank);
  }

  DataHeap::getInstance().receiveData(
      receivedTimeStepData.data(),receivedTimeStepData.size(),workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);

  assertion1(receivedTimeStepData.size()==4,receivedTimeStepData.size());
  assertion1(receivedTimeStepData[0]>=0,receivedTimeStepData[0]);
  assertion1(std::isfinite(receivedTimeStepData[0]),receivedTimeStepData[0]);
  // The master solver has not yet updated its minNextPredictorTimeStepSize.
  // Thus it does not equal MAX_DOUBLE.

  int index=0;
  _minNextPredictorTimeStepSize  = std::min( _minNextPredictorTimeStepSize, receivedTimeStepData[index++] );
  _nextGridUpdateRequested      |= (receivedTimeStepData[index++]) > 0 ? true : false;
  _nextMinCellSize               = std::min( _nextMinCellSize, receivedTimeStepData[index++] );
  _nextMaxCellSize               = std::max( _nextMaxCellSize, receivedTimeStepData[index++] );

  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","Receiving time step data: " <<
             "data[0]=" << receivedTimeStepData[0] <<
             ",data[1]=" << receivedTimeStepData[1] <<
             ",data[2]=" << receivedTimeStepData[2] <<
             ",data[3]=" << receivedTimeStepData[3] );

    logDebug("mergeWithWorkerData(...)","Updated time step fields: " <<
             "_minNextPredictorTimeStepSize=" << _minNextPredictorTimeStepSize <<
             "_nextGridUpdateRequested=" << _nextGridUpdateRequested <<
             ",_nextMinCellSize=" << _nextMinCellSize <<
             ",_nextMaxCellSize=" << _nextMaxCellSize);
  }
}

void exahype::solvers::ADERDGSolver::sendEmptyDataToMaster(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  logDebug("sendEmptyDataToMaster(...)","Empty data for solver sent to rank "<<masterRank<<
           ", cell: "<< x << ", level: " << level);

  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerMasterWorkerCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::ADERDGSolver::sendDataToMaster(
    const int                                     masterRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  CellDescription& cellDescription = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  if (cellDescription.getType()==CellDescription::Ancestor) {
    double* extrapolatedPredictor = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data();
    double* fluctuations          = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data();

    logDebug("sendDataToMaster(...)","Face data of solver " << cellDescription.getSolverNumber() << " sent to rank "<<masterRank<<
             ", cell: "<< x << ", level: " << level);

    // No inverted message order since we do synchronous data exchange.
    // Order: extrapolatedPredictor,fluctuations.
    DataHeap::getInstance().sendData(
        extrapolatedPredictor, getBndTotalSize(), masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
    DataHeap::getInstance().sendData(
        fluctuations, getBndFluxTotalSize(), masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
  } else {
    sendEmptyDataToMaster(masterRank,x,level);
  }
}

void exahype::solvers::ADERDGSolver::mergeWithWorkerData(
    const int                                     workerRank,
    const MetadataHeap::HeapEntries&              workerMetadata,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  logDebug("mergeWithWorkerData(...)","Merge with worker data from rank "<<workerRank<<
             ", cell: "<< x << ", level: " << level);

  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  CellDescription& cellDescription = Heap::getInstance().getData(cellDescriptionsIndex)[element];
  const CellDescription::Type workerType =
      static_cast<CellDescription::Type>(workerMetadata[exahype::MetadataCellType].getU());

  // The following two assertions assert that cell descriptions on both ranks are together of
  // type Cell, Descendant, EmptyAncestor, or Ancestor.
  // Pairwise differing EmptyAncestor-Ancestor configurations as well as EmptyDescendants are not allowed.
  assertion4(cellDescription.getType()==CellDescription::Type::Cell
             || cellDescription.getType()==CellDescription::Type::Ancestor
             || cellDescription.getType()==CellDescription::Type::EmptyAncestor
             || cellDescription.getType()==CellDescription::Type::Descendant,
             cellDescription.toString(),CellDescription::toString(workerType),tarch::parallel::Node::getInstance().getRank(),workerRank);
  assertion4(cellDescription.getType()!=CellDescription::Type::EmptyDescendant
               && workerType!=CellDescription::Type::EmptyDescendant
               && !(workerType==CellDescription::Type::EmptyAncestor
                   && cellDescription.getType()==CellDescription::Type::Ancestor)
                   && !(workerType==CellDescription::Type::Ancestor &&
                       cellDescription.getType()==CellDescription::Type::EmptyAncestor),
                       cellDescription.toString(),CellDescription::toString(workerType),
                       tarch::parallel::Node::getInstance().getRank(),workerRank);

  if (cellDescription.getType()==CellDescription::Type::Ancestor) {
    logDebug("mergeWithWorkerData(...)","Received face data for solver " <<
             cellDescription.getSolverNumber() << " from Rank "<<workerRank<<
             ", cell: "<< x << ", level: " << level);

    // No inverted message order since we do synchronous data exchange.
    // Order: extrapolatedPredictor,fluctuations.
    // Make sure you clear the arrays before you append(!) data via receive!
    DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).clear();
    DataHeap::getInstance().getData(cellDescription.getFluctuation()).clear();

    DataHeap::getInstance().receiveData(
        cellDescription.getExtrapolatedPredictor(), workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
    DataHeap::getInstance().receiveData(
        cellDescription.getFluctuation(), workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);

    exahype::solvers::Solver::SubcellPosition subcellPosition =
        computeSubcellPositionOfCellOrAncestor(cellDescriptionsIndex,element);

    // TODO(Dominic): Add to docu. I can be the top most Ancestor too.
    if (subcellPosition.parentElement!=exahype::solvers::Solver::NotFound) {
      #if defined(Debug)
      CellDescription& parentCellDescription =
          getCellDescription(subcellPosition.parentCellDescriptionsIndex,subcellPosition.parentElement);
      #endif

      logDebug("mergeWithWorkerData(...)","Restricting face data for solver " <<
               cellDescription.getSolverNumber() << " from Rank "<<workerRank<<
               " from cell="<< x << ", level=" << level <<
               " to cell="<<parentCellDescription.getOffset()+0.5*parentCellDescription.getSize() <<
               " level="<<parentCellDescription.getLevel());

      restrictData(
          cellDescriptionsIndex,element,
          subcellPosition.parentCellDescriptionsIndex,subcellPosition.parentElement,
          subcellPosition.subcellIndex);
    }
  } else  {
    dropWorkerData(workerRank,x,level);
  }
}

void exahype::solvers::ADERDGSolver::dropWorkerData(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  logDebug("dropWorkerData(...)","Dropping worker data from rank "<<workerRank<<
             ", cell: "<< x << ", level: " << level);

  for(int receives=0; receives<DataMessagesPerMasterWorkerCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
void exahype::solvers::ADERDGSolver::sendDataToWorker(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  std::vector<double> timeStepDataToSend(0,6);
  timeStepDataToSend.push_back(_minCorrectorTimeStamp);  // TODO(Dominic): Append previous time step size
  timeStepDataToSend.push_back(_minCorrectorTimeStepSize);
  timeStepDataToSend.push_back(_minPredictorTimeStamp);
  timeStepDataToSend.push_back(_minPredictorTimeStepSize);

  timeStepDataToSend.push_back(_minCellSize);
  timeStepDataToSend.push_back(_maxCellSize);

  assertion1(timeStepDataToSend.size()==6,timeStepDataToSend.size());
  assertion1(std::isfinite(timeStepDataToSend[0]),timeStepDataToSend[0]);
  assertion1(std::isfinite(timeStepDataToSend[1]),timeStepDataToSend[1]);
  assertion1(std::isfinite(timeStepDataToSend[2]),timeStepDataToSend[2]);
  assertion1(std::isfinite(timeStepDataToSend[3]),timeStepDataToSend[3]);

  if (_timeStepping==TimeStepping::Global) {
    assertionEquals1(_minNextPredictorTimeStepSize,std::numeric_limits<double>::max(),
                                   tarch::parallel::Node::getInstance().getRank());
  }

  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToWorker(...)","Broadcasting time step data: " <<
            " data[0]=" << timeStepDataToSend[0] <<
            ",data[1]=" << timeStepDataToSend[1] <<
            ",data[2]=" << timeStepDataToSend[2] <<
            ",data[3]=" << timeStepDataToSend[3] <<
            ",data[4]=" << timeStepDataToSend[4] <<
            ",data[5]=" << timeStepDataToSend[5]);
    logDebug("sendDataWorker(...)","_minNextPredictorTimeStepSize="<<_minNextPredictorTimeStepSize);
  }

  DataHeap::getInstance().sendData(
      timeStepDataToSend.data(), timeStepDataToSend.size(),
      workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::ADERDGSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  std::vector<double> receivedTimeStepData(6);
  DataHeap::getInstance().receiveData(
      receivedTimeStepData.data(),receivedTimeStepData.size(),masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
  assertion1(receivedTimeStepData.size()==6,receivedTimeStepData.size());

  if (_timeStepping==TimeStepping::Global) {
    assertionNumericalEquals1(_minNextPredictorTimeStepSize,std::numeric_limits<double>::max(),
                                  _minNextPredictorTimeStepSize);
  }

  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithMasterData(...)","Received time step data: " <<
            "data[0]="  << receivedTimeStepData[0] <<
            ",data[1]=" << receivedTimeStepData[1] <<
            ",data[2]=" << receivedTimeStepData[2] <<
            ",data[3]=" << receivedTimeStepData[3] <<
            ",data[4]=" << receivedTimeStepData[4] <<
            ",data[5]=" << receivedTimeStepData[5]);
  }

  _minCorrectorTimeStamp    = receivedTimeStepData[0];
  _minCorrectorTimeStepSize = receivedTimeStepData[1];
  _minPredictorTimeStamp    = receivedTimeStepData[2];
  _minPredictorTimeStepSize = receivedTimeStepData[3];

  _minCellSize              = receivedTimeStepData[4];
  _maxCellSize              = receivedTimeStepData[5];
}

bool exahype::solvers::ADERDGSolver::hasToSendDataToMaster(
    const int cellDescriptionsIndex,
    const int element) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  CellDescription& cellDescription = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  if (cellDescription.getType()==CellDescription::Ancestor) {
    return true;
  } else if (cellDescription.getType()==CellDescription::EmptyAncestor) {
    #if defined(Debug) || defined(Asserts)
    exahype::solvers::Solver::SubcellPosition subcellPosition =
        computeSubcellPositionOfCellOrAncestor(cellDescriptionsIndex,element);
    assertion(subcellPosition.parentElement==exahype::solvers::Solver::NotFound);
    #endif
  }

  return false;
}

void exahype::solvers::ADERDGSolver::sendDataToWorker(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  CellDescription& cellDescription = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  if (cellDescription.getType()==CellDescription::Descendant) {
    exahype::solvers::Solver::SubcellPosition subcellPosition =
        exahype::amr::computeSubcellPositionOfDescendant<CellDescription,Heap,false>(cellDescription);
    prolongateFaceDataToDescendant(cellDescription,subcellPosition);

    double* extrapolatedPredictor = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data();
    double* fluctuations          = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data();

    // No inverted message order since we do synchronous data exchange.
    // Order: extraplolatedPredictor,fluctuations.
    DataHeap::getInstance().sendData(
        extrapolatedPredictor, getBndTotalSize(), workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
    DataHeap::getInstance().sendData(
        fluctuations, getBndFluxTotalSize(), workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);

    logDebug("sendDataToWorker(...)","Sent face data of solver " <<
             cellDescription.getSolverNumber() << " to rank "<< workerRank <<
             ", cell: "<< x << ", level: " << level);
  } else {
    sendEmptyDataToWorker(workerRank,x,level);
  }
}

void exahype::solvers::ADERDGSolver::sendEmptyDataToWorker(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerMasterWorkerCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, workerRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::ADERDGSolver::mergeWithMasterData(
    const int                                     masterRank,
    const MetadataHeap::HeapEntries&              masterMetadata,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  CellDescription& cellDescription = Heap::getInstance().getData(cellDescriptionsIndex)[element];
  const CellDescription::Type masterType =
      static_cast<CellDescription::Type>(masterMetadata[exahype::MetadataCellType].getU());

  // The following two assertions assert that cell descriptions on both ranks are together of
  // type Cell, Descendant, EmptyAncestor, or Ancestor.
  // Pairwise differing EmptyAncestor-Ancestor configurations as well as EmptyDescendants are not allowed.
  assertion4(cellDescription.getType()==CellDescription::Type::Cell
             || cellDescription.getType()==CellDescription::Type::Ancestor
             || cellDescription.getType()==CellDescription::Type::EmptyAncestor
             || cellDescription.getType()==CellDescription::Type::Descendant,
             cellDescription.toString(),CellDescription::toString(masterType),
             tarch::parallel::Node::getInstance().getRank(),masterRank);
  assertion4(cellDescription.getType()!=CellDescription::Type::EmptyDescendant
             && masterType!=CellDescription::Type::EmptyDescendant
             && !(masterType==CellDescription::Type::EmptyAncestor
                 && cellDescription.getType()==CellDescription::Type::Ancestor)
                 && !(masterType==CellDescription::Type::Ancestor &&
                     cellDescription.getType()==CellDescription::Type::EmptyAncestor),
                     cellDescription.toString(),CellDescription::toString(masterType),
                     tarch::parallel::Node::getInstance().getRank(),masterRank);

  if (cellDescription.getType()==CellDescription::Descendant) {
    logDebug("mergeWithMasterData(...)","Received face data for solver " <<
             cellDescription.getSolverNumber() << " from rank "<<masterRank<<
             ", cell: "<< x << ", level: " << level);

    // No inverted send and receives order since we do synchronous data exchange.
    // Order: extraplolatedPredictor,fluctuations
    DataHeap::getInstance().receiveData(
        cellDescription.getExtrapolatedPredictor(), masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
    DataHeap::getInstance().receiveData(
        cellDescription.getFluctuation(), masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
  } else {
    dropMasterData(masterRank,x,level);
  }
}

void exahype::solvers::ADERDGSolver::dropMasterData(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
        const int                                     level) {
  for(int receives=0; receives<DataMessagesPerMasterWorkerCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        masterRank, x, level,
        peano::heap::MessageType::MasterWorkerCommunication);
}
#endif

std::string exahype::solvers::ADERDGSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::ADERDGSolver::toString (std::ostream& out) const {
  out << "(";
  out << "_identifier:" << _identifier;
  out << ",";
  out << "_type:" << exahype::solvers::Solver::toString(_type);
  out << ",";
  out << "_numberOfVariables:" << _numberOfVariables;
  out << ",";
  out << "_numberOfParameters:" << _numberOfParameters;
  out << ",";
  out << "_nodesPerCoordinateAxis:" << _nodesPerCoordinateAxis;
  out << ",";
  out << "_maximumMeshSize:" << _maximumMeshSize;
  out << ",";
  out << "_timeStepping:" << exahype::solvers::Solver::toString(_timeStepping); // only solver attributes
  out << ",";
  out << "_unknownsPerFace:" << _dofPerFace;
  out << ",";
  out << "_unknownsPerCellBoundary:" << _dofPerCellBoundary;
  out << ",";
  out << "_unknownsPerCell:" << _dofPerCell;
  out << ",";
  out << "_fluxUnknownsPerCell:" << _fluxDofPerCell;
  out << ",";
  out << "_spaceTimeUnknownsPerCell:" << _spaceTimeDofPerCell;
  out << ",";
  out << "_spaceTimeFluxUnknownsPerCell:" << _spaceTimeFluxDofPerCell;
  out << ",";
  out << "_previousMinCorrectorTimeStepSize:" << _previousMinCorrectorTimeStepSize;
  out << ",";
  out << "_minCorrectorTimeStamp:" << _minCorrectorTimeStamp;
  out << ",";
  out << "_minCorrectorTimeStepSize:" << _minCorrectorTimeStepSize;
  out << ",";
  out << "_minPredictorTimeStepSize:" << _minPredictorTimeStepSize;
  out << ",";
  out << "_minNextPredictorTimeStepSize:" << _minNextPredictorTimeStepSize;
  out <<  ")";
}


exahype::solvers::ADERDGSolver::CompressionTask::CompressionTask(
  ADERDGSolver&                             solver,
  exahype::records::ADERDGCellDescription&  cellDescription
):
  _solver(solver),
  _cellDescription(cellDescription) {
}


void exahype::solvers::ADERDGSolver::CompressionTask::operator()() {
  _solver.determineUnknownAverages(_cellDescription);
  _solver.computeHierarchicalTransform(_cellDescription,-1.0);
  _solver.putUnknownsIntoByteStream(_cellDescription);

  tarch::multicore::Lock lock(_heapSemaphore);
  _cellDescription.setCompressionState(exahype::records::ADERDGCellDescription::Compressed);
  _NumberOfTriggeredTasks--;
  assertion( _NumberOfTriggeredTasks>=0 );
}


void exahype::solvers::ADERDGSolver::compress(exahype::records::ADERDGCellDescription& cellDescription) {
  assertion1( cellDescription.getCompressionState() ==  exahype::records::ADERDGCellDescription::Uncompressed, cellDescription.toString() );
  if (CompressionAccuracy>0.0) {
    if (SpawnCompressionAsBackgroundThread) {
      cellDescription.setCompressionState(exahype::records::ADERDGCellDescription::CurrentlyProcessed);

      tarch::multicore::Lock lock(_heapSemaphore);
      _NumberOfTriggeredTasks++;
      lock.free();

      CompressionTask myTask( *this, cellDescription );
      peano::datatraversal::TaskSet spawnedSet( myTask );
    }
    else {
      determineUnknownAverages(cellDescription);
      computeHierarchicalTransform(cellDescription,-1.0);
      putUnknownsIntoByteStream(cellDescription);
      cellDescription.setCompressionState(exahype::records::ADERDGCellDescription::Compressed);
    }
  }
}


void exahype::solvers::ADERDGSolver::uncompress(exahype::records::ADERDGCellDescription& cellDescription) {
  #ifdef SharedMemoryParallelisation
  bool madeDecision = CompressionAccuracy<=0.0;
  bool uncompress   = false;

  while (!madeDecision) {
    tarch::multicore::Lock lock(_heapSemaphore);
    madeDecision = cellDescription.getCompressionState() != exahype::records::ADERDGCellDescription::CurrentlyProcessed;
    uncompress   = cellDescription.getCompressionState() == exahype::records::ADERDGCellDescription::Compressed;
    if (uncompress) {
      cellDescription.setCompressionState( exahype::records::ADERDGCellDescription::CurrentlyProcessed );
    }
    lock.free();

    tarch::multicore::BooleanSemaphore::sendTaskToBack();
  }
  #else
  bool uncompress = CompressionAccuracy>0.0
      && cellDescription.getCompressionState() == exahype::records::ADERDGCellDescription::Compressed;
  #endif

/*
  #ifdef Parallel
  assertion1(!cellDescription.getAdjacentToRemoteRank() || cellDescription.getCompressionState() == exahype::records::ADERDGCellDescription::Compressed,
             cellDescription.toString());
  #endif
*/

  if (uncompress) {
    pullUnknownsFromByteStream(cellDescription);
    computeHierarchicalTransform(cellDescription,1.0);

    tarch::multicore::Lock lock(_heapSemaphore);
    cellDescription.setCompressionState(exahype::records::ADERDGCellDescription::Uncompressed);
  }
}


void exahype::solvers::ADERDGSolver::determineUnknownAverages(
  exahype::records::ADERDGCellDescription& cellDescription
) {
  for (int variableNumber=0; variableNumber<getNumberOfVariables()+getNumberOfParameters(); variableNumber++) {
    double solutionAverage         = 0.0;
    double previousSolutionAverage = 0.0;
    double updateAverage           = 0.0;

    int    numberOfDoFsPerVariable  = power(getNodesPerCoordinateAxis(), DIMENSIONS);
    for (int i=0; i<numberOfDoFsPerVariable; i++) {
      solutionAverage         += DataHeap::getInstance().getData( cellDescription.getSolution()         )[i + variableNumber * numberOfDoFsPerVariable];
      previousSolutionAverage += DataHeap::getInstance().getData( cellDescription.getPreviousSolution() )[i + variableNumber * numberOfDoFsPerVariable];
      updateAverage           += DataHeap::getInstance().getData( cellDescription.getUpdate()           )[i + variableNumber * numberOfDoFsPerVariable];
    }
    DataHeap::getInstance().getData( cellDescription.getSolutionAverages()         )[variableNumber] = solutionAverage         / numberOfDoFsPerVariable;
    DataHeap::getInstance().getData( cellDescription.getPreviousSolutionAverages() )[variableNumber] = previousSolutionAverage / numberOfDoFsPerVariable;
    DataHeap::getInstance().getData( cellDescription.getUpdateAverages()           )[variableNumber] = updateAverage           / numberOfDoFsPerVariable;

    for (int face=0; face<2*DIMENSIONS; face++) {
      double extrapolatedPredictorAverage = 0.0;
      double boundaryFluxAverage          = 0.0;
      int    numberOfDoFsPerVariable  = power(getNodesPerCoordinateAxis(), DIMENSIONS-1);
      for (int i=0; i<numberOfDoFsPerVariable; i++) {
        extrapolatedPredictorAverage += DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictor() )[i + variableNumber * numberOfDoFsPerVariable + getNumberOfVariables() * numberOfDoFsPerVariable * face ];
        boundaryFluxAverage          += DataHeap::getInstance().getData( cellDescription.getFluctuation()           )[i + variableNumber * numberOfDoFsPerVariable + getNumberOfVariables() * numberOfDoFsPerVariable * face ];
      }
      DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictorAverages() )[variableNumber+getNumberOfVariables() * face] = extrapolatedPredictorAverage / numberOfDoFsPerVariable;
      DataHeap::getInstance().getData( cellDescription.getFluctuationAverages()           )[variableNumber+getNumberOfVariables() * face] = boundaryFluxAverage   / numberOfDoFsPerVariable;
    }
  }
}


void exahype::solvers::ADERDGSolver::computeHierarchicalTransform(exahype::records::ADERDGCellDescription& cellDescription, double sign) {
  for (int variableNumber=0; variableNumber<getNumberOfVariables()+getNumberOfParameters(); variableNumber++) {
    int    numberOfDoFsPerVariable  = power(getNodesPerCoordinateAxis(), DIMENSIONS);
    for (int i=0; i<numberOfDoFsPerVariable; i++) {
      DataHeap::getInstance().getData( cellDescription.getSolution()         )[i + variableNumber * numberOfDoFsPerVariable] += sign * DataHeap::getInstance().getData( cellDescription.getSolutionAverages() )[variableNumber];
      DataHeap::getInstance().getData( cellDescription.getPreviousSolution() )[i + variableNumber * numberOfDoFsPerVariable] += sign * DataHeap::getInstance().getData( cellDescription.getPreviousSolutionAverages() )[variableNumber];
      DataHeap::getInstance().getData( cellDescription.getUpdate()           )[i + variableNumber * numberOfDoFsPerVariable] += sign * DataHeap::getInstance().getData( cellDescription.getUpdateAverages()   )[variableNumber];
    }

    for (int face=0; face<2*DIMENSIONS; face++) {
      int    numberOfDoFsPerVariable  = power(getNodesPerCoordinateAxis(), DIMENSIONS-1);
      for (int i=0; i<numberOfDoFsPerVariable; i++) {
        DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictor() )[i + variableNumber * numberOfDoFsPerVariable + getNumberOfVariables() * numberOfDoFsPerVariable * face ] += sign * DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictorAverages() )[variableNumber+getNumberOfVariables() * face];
        DataHeap::getInstance().getData( cellDescription.getFluctuation()           )[i + variableNumber * numberOfDoFsPerVariable + getNumberOfVariables() * numberOfDoFsPerVariable * face ] += sign * DataHeap::getInstance().getData( cellDescription.getFluctuationAverages()           )[variableNumber+getNumberOfVariables() * face];
      }
    }
  }
}


void exahype::solvers::ADERDGSolver::tearApart(int numberOfEntries, int normalHeapIndex, int compressedHeapIndex, int bytesForMantissa) {
  char exponent;
  long int mantissa;
  char* pMantissa = reinterpret_cast<char*>( &(mantissa) );

  assertion( DataHeap::getInstance().isValidIndex(normalHeapIndex) );
  assertion( CompressedDataHeap::getInstance().isValidIndex(compressedHeapIndex) );
  assertion2( static_cast<int>(DataHeap::getInstance().getData(normalHeapIndex).size())==numberOfEntries, DataHeap::getInstance().getData(normalHeapIndex).size(), numberOfEntries );
  assertion( CompressedDataHeap::getInstance().getData(compressedHeapIndex).empty() );

  CompressedDataHeap::getInstance().getData( compressedHeapIndex ).resize(numberOfEntries * (bytesForMantissa+1));

  int compressedDataHeapIndex = 0;
  for (int i=0; i<numberOfEntries; i++) {
    peano::heap::decompose(
      DataHeap::getInstance().getData( normalHeapIndex )[i],
      exponent, mantissa, bytesForMantissa
    );
    CompressedDataHeap::getInstance().getData( compressedHeapIndex )[compressedDataHeapIndex]._persistentRecords._u = exponent;
    compressedDataHeapIndex++;
    for (int j=0; j<bytesForMantissa; j++) {
      CompressedDataHeap::getInstance().getData( compressedHeapIndex )[compressedDataHeapIndex]._persistentRecords._u = pMantissa[j];
      compressedDataHeapIndex++;
    }
  }
}


void exahype::solvers::ADERDGSolver::glueTogether(int numberOfEntries, int normalHeapIndex, int compressedHeapIndex, int bytesForMantissa) {
  char exponent  = 0;
  long int mantissa;
  char* pMantissa = reinterpret_cast<char*>( &(mantissa) );

  assertion( DataHeap::getInstance().isValidIndex(normalHeapIndex) );
  assertion( CompressedDataHeap::getInstance().isValidIndex(compressedHeapIndex) );
  assertion5(
    static_cast<int>(CompressedDataHeap::getInstance().getData(compressedHeapIndex).size())==numberOfEntries * (bytesForMantissa+1),
    CompressedDataHeap::getInstance().getData(compressedHeapIndex).size(), numberOfEntries * (bytesForMantissa+1),
    numberOfEntries, compressedHeapIndex, bytesForMantissa
  );

  #ifdef ValidateCompressedVsUncompressedData
  assertion( static_cast<int>(DataHeap::getInstance().getData(normalHeapIndex).size())==numberOfEntries );
  #else
  DataHeap::getInstance().getData(normalHeapIndex).resize(numberOfEntries);
  #endif

  int compressedDataHeapIndex = numberOfEntries * (bytesForMantissa+1)-1;
  for (int i=numberOfEntries-1; i>=0; i--) {
    mantissa = 0;
    for (int j=bytesForMantissa-1; j>=0; j--) {
      pMantissa[j] = CompressedDataHeap::getInstance().getData( compressedHeapIndex )[compressedDataHeapIndex]._persistentRecords._u; // TODO(Dominic):This line fails
      compressedDataHeapIndex--;
    }
    exponent = CompressedDataHeap::getInstance().getData( compressedHeapIndex )[compressedDataHeapIndex]._persistentRecords._u;
    compressedDataHeapIndex--;
    double reconstructedValue = peano::heap::compose(
      exponent, mantissa, bytesForMantissa
    );
    #ifdef ValidateCompressedVsUncompressedData
    assertion7(
      tarch::la::equals( DataHeap::getInstance().getData(normalHeapIndex)[i], reconstructedValue, CompressionAccuracy ),
      DataHeap::getInstance().getData(normalHeapIndex)[i], reconstructedValue, DataHeap::getInstance().getData(normalHeapIndex)[i] - reconstructedValue,
      CompressionAccuracy, bytesForMantissa, numberOfEntries, normalHeapIndex
    );
    #else
    DataHeap::getInstance().getData(normalHeapIndex)[i] = reconstructedValue;
    #endif
  }
}



void exahype::solvers::ADERDGSolver::putUnknownsIntoByteStream(exahype::records::ADERDGCellDescription& cellDescription) {
  assertion(CompressionAccuracy>0.0);

  assertion( cellDescription.getPreviousSolutionCompressed()==-1 );
  assertion( cellDescription.getSolutionCompressed()==-1 );
  assertion( cellDescription.getUpdateCompressed()==-1 );
  assertion( cellDescription.getExtrapolatedPredictorCompressed()==-1 );
  assertion( cellDescription.getFluctuationCompressed()==-1 );

  int compressionOfPreviousSolution;
  int compressionOfSolution;
  int compressionOfUpdate;
  int compressionOfExtrapolatedPredictor;
  int compressionOfFluctuation;

  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getUpdate() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictor() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getFluctuation() ));

  peano::datatraversal::TaskSet compressionFactorIdentification(
    [&]() -> void  { compressionOfPreviousSolution = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).data(),
      getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS),
      CompressionAccuracy
      );},
    [&] () -> void  { compressionOfSolution = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getSolution() ).data(),
      getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS),
      CompressionAccuracy
      );},
    [&]() -> void  { compressionOfUpdate = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getUpdate() ).data(),
      getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS),
      CompressionAccuracy
      );},
    [&]() -> void  { compressionOfExtrapolatedPredictor = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictor() ).data(),
      getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS-1) * 2 * DIMENSIONS,
      CompressionAccuracy
      );},
    [&]() -> void  { compressionOfFluctuation = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getFluctuation() ).data(),
      getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS-1) * 2 * DIMENSIONS,
      CompressionAccuracy
      );},
      true
  );

  assertion(1<=compressionOfPreviousSolution);
  assertion(1<=compressionOfSolution);
  assertion(1<=compressionOfUpdate);
  assertion(1<=compressionOfExtrapolatedPredictor);
  assertion(1<=compressionOfFluctuation);

  assertion(compressionOfPreviousSolution<=7);
  assertion(compressionOfSolution<=7);
  assertion(compressionOfUpdate<=7);
  assertion(compressionOfExtrapolatedPredictor<=7);
  assertion(compressionOfFluctuation<=7);

  peano::datatraversal::TaskSet runParallelTasks(
    [&]() -> void {
      cellDescription.setBytesPerDoFInPreviousSolution(compressionOfPreviousSolution);
      if (compressionOfPreviousSolution<7) {
        tarch::multicore::Lock lock(_heapSemaphore);
        cellDescription.setPreviousSolutionCompressed( CompressedDataHeap::getInstance().createData(0,0,CompressedDataHeap::Allocation::UseOnlyRecycledEntries) );
        assertion( cellDescription.getPreviousSolutionCompressed()>=0 );
        lock.free();

        const int numberOfEntries = getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS);
        tearApart(numberOfEntries, cellDescription.getPreviousSolution(), cellDescription.getPreviousSolutionCompressed(), compressionOfPreviousSolution);

        #if defined(Asserts)
        lock.lock();
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).size() * 8.0;
        PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getPreviousSolutionCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
        DataHeap::getInstance().deleteData( cellDescription.getPreviousSolution(), true );
        cellDescription.setPreviousSolution( -1 );
        #endif
      }
      else {
        #if defined(Asserts)
        tarch::multicore::Lock lock(_heapSemaphore);
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).size() * 8.0;
        PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).size() * 8.0;
        #endif
      }
    },
    [&]() -> void {
      cellDescription.setBytesPerDoFInSolution(compressionOfSolution);
      if (compressionOfSolution<7) {
        tarch::multicore::Lock lock(_heapSemaphore);
        cellDescription.setSolutionCompressed(CompressedDataHeap::getInstance().createData(0,0,CompressedDataHeap::Allocation::UseOnlyRecycledEntries));
        assertion1( cellDescription.getSolutionCompressed()>=0, cellDescription.getSolutionCompressed() );
        lock.free();

        const int numberOfEntries = getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS);

        tearApart(numberOfEntries, cellDescription.getSolution(), cellDescription.getSolutionCompressed(), compressionOfSolution);

        #if defined(Asserts)
        lock.lock();
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getSolution() ).size() * 8.0;
        PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getSolutionCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
        DataHeap::getInstance().deleteData( cellDescription.getSolution(), true );
        cellDescription.setSolution( -1 );
        #endif
      }
      else {
        #if defined(Asserts)
        tarch::multicore::Lock lock(_heapSemaphore);
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getSolution() ).size() * 8.0;
        PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getSolution() ).size() * 8.0;
        #endif
      }
    },
    [&]() -> void {
      cellDescription.setBytesPerDoFInUpdate(compressionOfUpdate);
      if (compressionOfUpdate<7) {
        tarch::multicore::Lock lock(_heapSemaphore);
        cellDescription.setUpdateCompressed( CompressedDataHeap::getInstance().createData(0,0,CompressedDataHeap::Allocation::UseOnlyRecycledEntries) );
        assertion( cellDescription.getUpdateCompressed()>=0 );
        lock.free();

        const int numberOfEntries = getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS);
        tearApart(numberOfEntries, cellDescription.getUpdate(), cellDescription.getUpdateCompressed(), compressionOfUpdate);

        #if defined(Asserts)
        lock.lock();
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getUpdate() ).size() * 8.0;
        PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getUpdateCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
        DataHeap::getInstance().deleteData( cellDescription.getUpdate(), true );
        cellDescription.setUpdate( -1 );
        #endif
      }
      else {
        #if defined(Asserts)
        tarch::multicore::Lock lock(_heapSemaphore);
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getUpdate() ).size() * 8.0;
        PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getUpdate() ).size() * 8.0;
        #endif
      }
    },
    [&]() -> void {
      cellDescription.setBytesPerDoFInExtrapolatedPredictor(compressionOfExtrapolatedPredictor);
      if (compressionOfExtrapolatedPredictor<7) {
        tarch::multicore::Lock lock(_heapSemaphore);
        cellDescription.setExtrapolatedPredictorCompressed( CompressedDataHeap::getInstance().createData(0,0,CompressedDataHeap::Allocation::UseOnlyRecycledEntries) );
        assertion( cellDescription.getExtrapolatedPredictorCompressed()>=0 );
        lock.free();

        const int numberOfEntries = getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS-1) * 2 * DIMENSIONS;
        tearApart(numberOfEntries, cellDescription.getExtrapolatedPredictor(), cellDescription.getExtrapolatedPredictorCompressed(), compressionOfExtrapolatedPredictor);

        #if defined(Asserts)
        lock.lock();
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictor() ).size() * 8.0;
        PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictorCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
        DataHeap::getInstance().deleteData( cellDescription.getExtrapolatedPredictor(), true );
        cellDescription.setExtrapolatedPredictor( -1 );
        #endif
      }
      else {
        #if defined(Asserts)
        tarch::multicore::Lock lock(_heapSemaphore);
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictor() ).size() * 8.0;
        PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictor() ).size() * 8.0;
        #endif
      }
    },
    [&]() -> void {
      cellDescription.setBytesPerDoFInFluctuation(compressionOfFluctuation);
      if (compressionOfFluctuation<7) {
        tarch::multicore::Lock lock(_heapSemaphore);
        cellDescription.setFluctuationCompressed( CompressedDataHeap::getInstance().createData(0,0,CompressedDataHeap::Allocation::UseOnlyRecycledEntries) );
        assertion( cellDescription.getFluctuationCompressed()>=0 );
        lock.free();

        const int numberOfEntries = getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS-1) * 2 * DIMENSIONS;
        tearApart(numberOfEntries, cellDescription.getFluctuation(), cellDescription.getFluctuationCompressed(), compressionOfFluctuation);

        #if defined(Asserts)
        lock.lock();
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getFluctuation() ).size() * 8.0;
        PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getFluctuationCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
        DataHeap::getInstance().deleteData( cellDescription.getFluctuation(), true );
        cellDescription.setFluctuation( -1 );
        #endif
      }
      else {
        #if defined(Asserts)
        tarch::multicore::Lock lock(_heapSemaphore);
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getFluctuation() ).size() * 8.0;
        PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getFluctuation() ).size() * 8.0;
        #endif
      }
    },
    true
  );
}


void exahype::solvers::ADERDGSolver::pullUnknownsFromByteStream(exahype::records::ADERDGCellDescription& cellDescription) {
  assertion(CompressionAccuracy>0.0);

  #if !defined(ValidateCompressedVsUncompressedData)
  const int unknownsPerCell         = getUnknownsPerCell();
  const int dataPointsPerCell       = getDataPerCell();
  const int unknownsPerCellBoundary = getUnknownsPerCellBoundary();

  {
    tarch::multicore::Lock lock(_heapSemaphore);
    cellDescription.setPreviousSolution( DataHeap::getInstance().createData( dataPointsPerCell,         dataPointsPerCell,         DataHeap::Allocation::UseOnlyRecycledEntries) );
    cellDescription.setSolution( DataHeap::getInstance().createData(         dataPointsPerCell,         dataPointsPerCell,         DataHeap::Allocation::UseOnlyRecycledEntries) );
    cellDescription.setUpdate( DataHeap::getInstance().createData(           unknownsPerCell,         unknownsPerCell,         DataHeap::Allocation::UseOnlyRecycledEntries) );

    cellDescription.setExtrapolatedPredictor( DataHeap::getInstance().createData( unknownsPerCellBoundary, unknownsPerCellBoundary, DataHeap::Allocation::UseOnlyRecycledEntries) );
    cellDescription.setFluctuation( DataHeap::getInstance().createData(           unknownsPerCellBoundary, unknownsPerCellBoundary, DataHeap::Allocation::UseOnlyRecycledEntries) );
    lock.free();

    if (cellDescription.getPreviousSolution()==-1) {
      waitUntilAllBackgroundTasksHaveTerminated();
      lock.lock();
      cellDescription.setPreviousSolution( DataHeap::getInstance().createData( dataPointsPerCell, dataPointsPerCell, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired) );
      lock.free();
    }
    if (cellDescription.getSolution()==-1) {
      waitUntilAllBackgroundTasksHaveTerminated();
      lock.lock();
      cellDescription.setSolution( DataHeap::getInstance().createData( dataPointsPerCell, dataPointsPerCell, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired) );
      lock.free();
    }
    if (cellDescription.getUpdate()==-1) {
      waitUntilAllBackgroundTasksHaveTerminated();
      lock.lock();
      cellDescription.setUpdate( DataHeap::getInstance().createData( unknownsPerCell, unknownsPerCell, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired) );
      lock.free();
    }
    if (cellDescription.getExtrapolatedPredictor()==-1) {
      waitUntilAllBackgroundTasksHaveTerminated();
      lock.lock();
      cellDescription.setExtrapolatedPredictor( DataHeap::getInstance().createData(unknownsPerCellBoundary, unknownsPerCellBoundary, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired) );
      lock.free();
    }
    if (cellDescription.getFluctuation()==-1) {
      waitUntilAllBackgroundTasksHaveTerminated();
      lock.lock();
      cellDescription.setFluctuation( DataHeap::getInstance().createData( unknownsPerCellBoundary, unknownsPerCellBoundary, DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired) );
      lock.free();
    }
  }
  #else
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getUpdate() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictor() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getFluctuation() ));
  #endif

  assertion1(
      CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressed() ),
      cellDescription.getPreviousSolutionCompressed()
    );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressed() ),
    cellDescription.getSolutionCompressed()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getUpdateCompressed() ),
    cellDescription.getUpdateCompressed()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictorCompressed() ),
    cellDescription.getExtrapolatedPredictorCompressed()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getFluctuationCompressed() ),
    cellDescription.getFluctuationCompressed()
  );

  peano::datatraversal::TaskSet glueTasks(
    [&]() -> void {
      if (cellDescription.getBytesPerDoFInPreviousSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolution() ), cellDescription.getPreviousSolution());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressed() ));
        const int numberOfEntries = ( getNumberOfVariables()+getNumberOfParameters() ) * power(getNodesPerCoordinateAxis(), DIMENSIONS);
        glueTogether(numberOfEntries, cellDescription.getPreviousSolution(), cellDescription.getPreviousSolutionCompressed(), cellDescription.getBytesPerDoFInPreviousSolution());
        tarch::multicore::Lock lock(_heapSemaphore);
        CompressedDataHeap::getInstance().deleteData( cellDescription.getPreviousSolutionCompressed(), true );
        cellDescription.setPreviousSolutionCompressed( -1 );
      }
    },
    [&]() -> void {
      if (cellDescription.getBytesPerDoFInSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getSolution() ), cellDescription.getSolution() );
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressed() ));
        const int numberOfEntries = ( getNumberOfVariables()+getNumberOfParameters() ) * power(getNodesPerCoordinateAxis(), DIMENSIONS);
        glueTogether(numberOfEntries, cellDescription.getSolution(), cellDescription.getSolutionCompressed(), cellDescription.getBytesPerDoFInSolution());
        tarch::multicore::Lock lock(_heapSemaphore);
        CompressedDataHeap::getInstance().deleteData( cellDescription.getSolutionCompressed(), true );
        cellDescription.setSolutionCompressed( -1 );
      }
    },
    [&]() -> void {
      if (cellDescription.getBytesPerDoFInUpdate()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getUpdate() ), cellDescription.getUpdate());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getUpdateCompressed() ));
        const int numberOfEntries = getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS);
        glueTogether(numberOfEntries, cellDescription.getUpdate(), cellDescription.getUpdateCompressed(), cellDescription.getBytesPerDoFInUpdate());
        tarch::multicore::Lock lock(_heapSemaphore);
        CompressedDataHeap::getInstance().deleteData( cellDescription.getUpdateCompressed(), true );
        cellDescription.setUpdateCompressed( -1 );
      }
    },
    [&]() -> void {
      if (cellDescription.getBytesPerDoFInExtrapolatedPredictor()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictor() ), cellDescription.getExtrapolatedPredictor());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictorCompressed() ));
        const int numberOfEntries = getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS-1) * 2 * DIMENSIONS;
        glueTogether(numberOfEntries, cellDescription.getExtrapolatedPredictor(), cellDescription.getExtrapolatedPredictorCompressed(), cellDescription.getBytesPerDoFInExtrapolatedPredictor());
        tarch::multicore::Lock lock(_heapSemaphore);
        CompressedDataHeap::getInstance().deleteData( cellDescription.getExtrapolatedPredictorCompressed(), true );
        cellDescription.setExtrapolatedPredictorCompressed( -1 );
      }
    },
    [&]() -> void {
      if (cellDescription.getBytesPerDoFInFluctuation()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getFluctuation() ), cellDescription.getFluctuation());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getFluctuationCompressed() ));
        const int numberOfEntries = getNumberOfVariables() * power(getNodesPerCoordinateAxis(), DIMENSIONS-1) * 2 * DIMENSIONS;
        glueTogether(numberOfEntries, cellDescription.getFluctuation(), cellDescription.getFluctuationCompressed(), cellDescription.getBytesPerDoFInFluctuation());
        tarch::multicore::Lock lock(_heapSemaphore);
        CompressedDataHeap::getInstance().deleteData( cellDescription.getFluctuationCompressed(), true );
        cellDescription.setFluctuationCompressed( -1 );
      }
    },
    true
  );
}

void exahype::solvers::ADERDGSolver::pointSource(
    const double t,
    const double dt, 
    const tarch::la::Vector<DIMENSIONS,double>& center,
    const tarch::la::Vector<DIMENSIONS,double>& dx, 
    double* tempPointForceSources) {}
