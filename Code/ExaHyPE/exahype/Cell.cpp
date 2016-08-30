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
 
#include "exahype/Cell.h"
#include "exahype/State.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "tarch/la/ScalarOperations.h"

#include "kernels/KernelCalls.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "exahype/records/ADERDGCellDescription.h"

tarch::logging::Log exahype::Cell::_log("exahype::Cell");

exahype::Cell::Cell() : Base() {
  // We initialise cells which are not touched by the
  // createCell(...) events of Peano's spacetree traversal automaton
  // with default ("do-nothing") values.
  _cellData.setCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

exahype::Cell::Cell(const Base::DoNotCallStandardConstructor& value)
    : Base(value) {
}

exahype::Cell::Cell(const Base::PersistentCell& argument) : Base(argument) {
  // This constructor is used to create a cell from persistent data.
  // Do not use it. This would overwrite persistent data.
}

void exahype::Cell::setupMetaData() {
  assertion1(!ADERDGCellDescriptionHeap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),toString());

  const int cellDescriptionIndex = ADERDGCellDescriptionHeap::getInstance().createData(0, 0);
  FiniteVolumesCellDescriptionHeap::getInstance().createDataForIndex(cellDescriptionIndex,0,0);

  _cellData.setCellDescriptionsIndex(cellDescriptionIndex);
}

void exahype::Cell::shutdownMetaData() {
  assertion1(
    ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      _cellData.getCellDescriptionsIndex()
    ),
    toString());

  ADERDGCellDescriptionHeap::getInstance().deleteData(_cellData.getCellDescriptionsIndex());
  FiniteVolumesCellDescriptionHeap::getInstance().deleteData(_cellData.getCellDescriptionsIndex());

  _cellData.setCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

#ifdef Parallel
exahype::MetadataHeap::HeapEntries exahype::Cell::encodeMetadata(int cellDescriptionsIndex) {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  const int nADERDG = ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex).size();
  const int nFV     = FiniteVolumesCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex).size();
  const int length  = 2 + 2*(nADERDG+nFV);

  exahype::MetadataHeap::HeapEntries encodedMetaData(0,length);
  // ADER-DG
  encodedMetaData.push_back(nADERDG); // Implicit conversion.
  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)) {
    encodedMetaData.push_back(p.getSolverNumber());           // Implicit conversion.
    encodedMetaData.push_back(static_cast<int>(p.getType())); // Implicit conversion.
  }
  // FV
  encodedMetaData.push_back(nFV);
  for (auto& p : FiniteVolumesCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)) {
    encodedMetaData.push_back(p.getSolverNumber());           // Implicit conversion.
    encodedMetaData.push_back(static_cast<int>(p.getType())); // Implicit conversion.
  }
  return encodedMetaData;
}

exahype::MetadataHeap::HeapEntries exahype::Cell::createEncodedMetadataSequenceForInvalidCellDescriptionsIndex() {
  exahype::MetadataHeap::HeapEntries metadata(0,2);
  metadata.push_back(0); // ADER-DG
  metadata.push_back(0); // FV
  return metadata;
}

bool exahype::Cell::isEncodedMetadataSequenceForInvalidCellDescriptionsIndex(exahype::MetadataHeap::HeapEntries& sequence) {
  return sequence.size()==2    &&
         sequence[0].getU()==0 && // ADER-DG
         sequence[1].getU()==0;   // FV
}

#endif

bool exahype::Cell::isInitialised() const {
  if (_cellData.getCellDescriptionsIndex()!=multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) {
    assertion1( ADERDGCellDescriptionHeap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),
            _cellData.getCellDescriptionsIndex());
    assertion1( FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),
        _cellData.getCellDescriptionsIndex());
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  return _cellData.getCellDescriptionsIndex()!=multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex;
}

int exahype::Cell::getCellDescriptionsIndex() const {
  return _cellData.getCellDescriptionsIndex();
}

void exahype::Cell::setCellDescriptionsIndex(int cellDescriptionsIndex) {
  _cellData.setCellDescriptionsIndex(cellDescriptionsIndex);
}

void exahype::Cell::addNewCellDescription(
    const int solverNumber,
    const exahype::records::FiniteVolumesCellDescription::Type cellType,
/*
    const exahype::records::FiniteVolumesCellDescription::RefinementEvent
        refinementEvent,
*/
    const int level, const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, double>& size,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre) {
  if (_cellData.getCellDescriptionsIndex() == multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) {
    setupMetaData();
  }

  assertion1(
    FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(
      _cellData.getCellDescriptionsIndex()
    ),
    toString());

  assertion2(static_cast<unsigned int>(solverNumber) <
                 solvers::RegisteredSolvers.size(),
             solverNumber, exahype::solvers::RegisteredSolvers.size());

//  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
//  assertion(solver->getType()==exahype::solvers::Solver::Type::FiniteVolumes);

  exahype::records::FiniteVolumesCellDescription newCellDescription;
  newCellDescription.setSolverNumber(solverNumber);

  // Default AMR settings
  newCellDescription.setType(cellType);
  newCellDescription.setLevel(level);
  //newCellDescription.setRefinementEvent(refinementEvent);

  // Pass geometry information to the cellDescription description
  newCellDescription.setSize(size);
  newCellDescription.setOffset(cellCentre);

  // Initialise helper variables
//  // TODO(Dominic):
//  newCellDescription.setHelperCellNeedsToStoreFaceData(false);
//  newCellDescription.setHelperCellNeedsToStoreFaceData(false);

  // Default field data indices
  newCellDescription.setSolution(-1);

  FiniteVolumesCellDescriptionHeap::getInstance()
      .getData(_cellData.getCellDescriptionsIndex())
      .push_back(newCellDescription);

}


void exahype::Cell::addNewCellDescription(
  const int                                     solverNumber,
  const exahype::records::ADERDGCellDescription::Type cellType,
  const exahype::records::ADERDGCellDescription::RefinementEvent refinementEvent,
  const int                                     level,
  const int                                     parentIndex,
  const tarch::la::Vector<DIMENSIONS, double>&  size,
  const tarch::la::Vector<DIMENSIONS, double>&  cellCentre) {
  if (_cellData.getCellDescriptionsIndex() == multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) {
    setupMetaData();
  }

  assertion1(
    ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      _cellData.getCellDescriptionsIndex()
    ),
    toString());

  assertion2(parentIndex == -1 ||
             parentIndex != _cellData.getCellDescriptionsIndex(),
             parentIndex, _cellData.getCellDescriptionsIndex());

  assertion2(parentIndex != _cellData.getCellDescriptionsIndex(),
             parentIndex, _cellData.getCellDescriptionsIndex());

  assertion2(static_cast<unsigned int>(solverNumber) <
                 solvers::RegisteredSolvers.size(),
             solverNumber, exahype::solvers::RegisteredSolvers.size());

//  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
//  assertion(solver->getType()==exahype::solvers::Solver::Type::ADER_DG);

  exahype::records::ADERDGCellDescription newCellDescription;
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
  newCellDescription.setSize(size);
  newCellDescription.setOffset(cellCentre);

  // Initialise helper variables
  #ifdef Parallel
  newCellDescription.setHelperCellNeedsToStoreFaceData(false);
  for (int i = 0; i < DIMENSIONS_TIMES_TWO; i++) {
    newCellDescription.setFaceDataExchangeCounter(i,TWO_POWER_D);
  }
  #endif

  // Default field data indices
  newCellDescription.setSolution(-1);
  newCellDescription.setUpdate(-1);
  newCellDescription.setExtrapolatedPredictor(-1);
  newCellDescription.setFluctuation(-1);

  // Limiter meta data (oscillations identificator)
  newCellDescription.setSolutionMin(-1);
  newCellDescription.setSolutionMax(-1);

  ADERDGCellDescriptionHeap::getInstance()
      .getData(_cellData.getCellDescriptionsIndex())
      .push_back(newCellDescription);
}

void exahype::Cell::ensureNecessaryMemoryIsAllocated(exahype::records::ADERDGCellDescription& cellDescription) {
  const solvers::Solver* solver = solvers::RegisteredSolvers[cellDescription.getSolverNumber()];
  assertion1(solver->getType()==exahype::solvers::Solver::Type::ADER_DG,cellDescription.toString());

  switch (cellDescription.getType()) {
    case exahype::records::ADERDGCellDescription::Cell:
      if (!DataHeap::getInstance().isValidIndex(cellDescription.getSolution())) {
        assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getUpdate()));
        // Allocate volume DoF for limiter
        const int unknownsPerCell = static_cast<const exahype::solvers::ADERDGSolver*>(solver)->getUnknownsPerCell();
        cellDescription.setUpdate(DataHeap::getInstance().createData(
            unknownsPerCell, unknownsPerCell));
        cellDescription.setSolution(DataHeap::getInstance().createData(
            unknownsPerCell, unknownsPerCell));
      }
      break;
    default:
      break;
  }

  // Allocate face DoF
  switch (cellDescription.getType()) {
    case exahype::records::ADERDGCellDescription::Cell:
    case exahype::records::ADERDGCellDescription::Ancestor:
    case exahype::records::ADERDGCellDescription::Descendant:
      if (!DataHeap::getInstance().isValidIndex(
          cellDescription.getExtrapolatedPredictor())) {
        assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()));

        // Allocate face DoF
        const int unknownsPerCellBoundary = static_cast<const exahype::solvers::ADERDGSolver*>(solver)->getUnknownsPerCellBoundary();
        cellDescription.setExtrapolatedPredictor(DataHeap::getInstance().createData(
            unknownsPerCellBoundary, unknownsPerCellBoundary));
        cellDescription.setFluctuation(DataHeap::getInstance().createData(
            unknownsPerCellBoundary, unknownsPerCellBoundary));

        // Allocate volume DoF for limiter
        const int unknownsPerCell = static_cast<const exahype::solvers::ADERDGSolver*>(solver)->getUnknownsPerCell();
        cellDescription.setSolutionMin(DataHeap::getInstance().createData(
            unknownsPerCell * 2 * DIMENSIONS, unknownsPerCell * 2 * DIMENSIONS));
        cellDescription.setSolutionMax(DataHeap::getInstance().createData(
            unknownsPerCell * 2 * DIMENSIONS, unknownsPerCell * 2 * DIMENSIONS));
        for (int i=0; i<unknownsPerCell * 2 * DIMENSIONS; i++) {
          DataHeap::getInstance().getData( cellDescription.getSolutionMax() )[i] = std::numeric_limits<double>::max();
          DataHeap::getInstance().getData( cellDescription.getSolutionMin() )[i] = std::numeric_limits<double>::min();
        }
      }
      break;
    default:
      break;
  }
}

void exahype::Cell::ensureNoUnnecessaryMemoryIsAllocated(exahype::records::ADERDGCellDescription& cellDescription) {
  #ifdef Asserts
  const solvers::Solver* solver = solvers::RegisteredSolvers[cellDescription.getSolverNumber()];
  assertion1(solver->getType()==exahype::solvers::Solver::Type::ADER_DG,cellDescription.toString());
  #endif

  if (DataHeap::getInstance().isValidIndex(cellDescription.getSolution())) {
    switch (cellDescription.getType()) {
      case exahype::records::ADERDGCellDescription::Erased:
      case exahype::records::ADERDGCellDescription::EmptyAncestor:
      case exahype::records::ADERDGCellDescription::EmptyDescendant:
      case exahype::records::ADERDGCellDescription::Ancestor:
      case exahype::records::ADERDGCellDescription::Descendant:
        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()));
        assertion(DataHeap::getInstance().isValidIndex(cellDescription.getUpdate()));

        DataHeap::getInstance().deleteData(cellDescription.getUpdate());
        DataHeap::getInstance().deleteData(cellDescription.getSolution());

        cellDescription.setSolution(-1);
        cellDescription.setUpdate(-1);
        break;
      default:
        break;
    }
  }

  if (DataHeap::getInstance().isValidIndex(
      cellDescription.getExtrapolatedPredictor())) {
    switch (cellDescription.getType()) {
      case exahype::records::ADERDGCellDescription::Erased:
      case exahype::records::ADERDGCellDescription::EmptyAncestor:
      case exahype::records::ADERDGCellDescription::EmptyDescendant:
        assertion(
            DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()));

        DataHeap::getInstance().deleteData(
            cellDescription.getExtrapolatedPredictor());
        DataHeap::getInstance().deleteData(cellDescription.getFluctuation());

        cellDescription.setExtrapolatedPredictor(-1);
        cellDescription.setFluctuation(-1);
        break;
      default:
        break;
    }
  }
}

void exahype::Cell::ensureNecessaryMemoryIsAllocated(const int solverNumber) {
  // Ensure that -1 value is invalid index (cf. addNewCellDescription)
  assertion(!DataHeap::getInstance().isValidIndex(-1));

  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 _cellData.getCellDescriptionsIndex()),
             toString());

  assertion3(static_cast<unsigned int>(solverNumber) <
                 solvers::RegisteredSolvers.size(),
             solverNumber, solvers::RegisteredSolvers.size(), toString());

  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];

  switch (solver->getType()) {
    case exahype::solvers::Solver::Type::ADER_DG:
      for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(
               _cellData.getCellDescriptionsIndex())) {
        if (solverNumber == p.getSolverNumber()) {
          ensureNecessaryMemoryIsAllocated(p);
        }
      }
      break;
    case exahype::solvers::Solver::Type::FiniteVolumes:
      for (auto& p : FiniteVolumesCellDescriptionHeap::getInstance().getData(_cellData.getCellDescriptionsIndex())) {
        if (solverNumber == p.getSolverNumber()) {
          switch (p.getType()) {
            case exahype::records::FiniteVolumesCellDescription::Cell:
              const int unknownsPerCell =
                  static_cast<const exahype::solvers::FiniteVolumesSolver*>(solver)->getUnknownsPerCell();
              assertion(unknownsPerCell>0);
              p.setSolution(DataHeap::getInstance().createData(unknownsPerCell, unknownsPerCell));
              logDebug( "initialiseCellDescription(...)", "allocated " << unknownsPerCell << " records for a cell" );
              break;
          }
        }
      }
      break;
    default:
      logDebug("initialiseCellDescription(...)",
               "solver is not associated with any cell descriptions of this "
               "cell. cell="
                   << toString());
      break;
  }
}



void exahype::Cell::ensureNoUnnecessaryMemoryIsAllocated(const int solverNumber) {
  // Ensure that -1 value is invalid index (cf. addNewCellDescription)
  assertion(!DataHeap::getInstance().isValidIndex(-1));

  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 _cellData.getCellDescriptionsIndex()),
             toString());

  assertion3(static_cast<unsigned int>(solverNumber) <
                 solvers::RegisteredSolvers.size(),
             solverNumber, solvers::RegisteredSolvers.size(), toString());

  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];

  switch (solver->getType()) {
    case exahype::solvers::Solver::Type::ADER_DG:
      for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(
               _cellData.getCellDescriptionsIndex())) {
        if (solverNumber == p.getSolverNumber()) {
          ensureNoUnnecessaryMemoryIsAllocated(p);
        }
      }
      break;
    case exahype::solvers::Solver::Type::FiniteVolumes:
      // TODO(Dominic): Implement!
      assertionMsg(false,"Not implemented yet!");
      break;
    default:
      logDebug("cleanCellDescription(...)",
               "solver is not associated with any cell descriptions of this "
               "cell. cell="
                   << toString());
      break;
  }
}


void exahype::Cell::validateNoNansInADERDGSolver(
  int                                  number,
  const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
  const std::string&                   methodTraceOfCaller
) {
  int unknownsPerCell              = 0;
  int unknownsPerCellBoundary      = 0;

#if defined(Debug) || defined(Asserts)
  auto& p = getADERDGCellDescription(number);

  double* luh = DataHeap::getInstance().getData(p.getSolution()).data();
  double* lduh = DataHeap::getInstance().getData(p.getUpdate()).data();

  double* lQhbnd = DataHeap::getInstance().getData(p.getExtrapolatedPredictor()).data();
  double* lFhbnd = DataHeap::getInstance().getData(p.getFluctuation()).data();

  exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
      exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);
  unknownsPerCell              = solver->getUnknownsPerCell();
  unknownsPerCellBoundary      = solver->getUnknownsPerCellBoundary();
  #endif


  assertion1(solver->getType()==exahype::solvers::Solver::Type::ADER_DG,p.toString());

  assertion1(DataHeap::getInstance().isValidIndex(p.getSolution()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getUpdate()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getExtrapolatedPredictor()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getFluctuation()),p.toString());

  assertionEquals4(p.getPredictorTimeStepSize(),p.getPredictorTimeStepSize(),
                   fineGridVerticesEnumerator.toString(),
                   p.toString(),toString(),methodTraceOfCaller);

  for (int i=0; i<unknownsPerCell; i++) {
    assertion5(std::isfinite(luh[i]), fineGridVerticesEnumerator.toString(),
        p.toString(),toString(),methodTraceOfCaller,i);
    assertion5(std::isfinite(lduh[i]), fineGridVerticesEnumerator.toString(),
        p.toString(),toString(),methodTraceOfCaller,i);
  }

  for (int i=0; i<unknownsPerCellBoundary; i++) {
    assertion5(std::isfinite(lQhbnd[i]), fineGridVerticesEnumerator.toString(),
        p.toString(),toString(),methodTraceOfCaller,i);
    assertion5(std::isfinite(lFhbnd[i]), fineGridVerticesEnumerator.toString(),
        p.toString(),toString(),methodTraceOfCaller,i);
  } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
}


exahype::Cell::SubcellPosition
exahype::Cell::computeSubcellPositionOfCellOrAncestor(
    const exahype::records::ADERDGCellDescription& pChild) const {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 _cellData.getCellDescriptionsIndex()),
             toString());
  assertion1(solvers::RegisteredSolvers[pChild.getSolverNumber()]->getType() ==
                 exahype::solvers::Solver::Type::ADER_DG,
             toString());
  assertion1(
      pChild.getType() == exahype::records::ADERDGCellDescription::Cell ||
          pChild.getType() == exahype::records::ADERDGCellDescription::Ancestor,
      toString());

  exahype::Cell::SubcellPosition subcellPosition;
  // Initialisation.
  subcellPosition.parentIndex = pChild.getParentIndex();
  exahype::records::ADERDGCellDescription* pParent = 0;
  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(
           subcellPosition.parentIndex)) {  // Loop over cell descriptions
    if (p.getSolverNumber() == pChild.getSolverNumber()) {
      pParent = &p;
      assertion(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
          subcellPosition.parentIndex));
    }
  }

  if (pParent != 0) {
    // Iterative determining of the top most parent that might hold data.
    while (pParent->getType() == exahype::records::ADERDGCellDescription::EmptyAncestor &&
           ADERDGCellDescriptionHeap::getInstance().isValidIndex(pParent->getParentIndex())) {
      const int currentParentIndex =
          pParent->getParentIndex();  // Value must be fixed. We update pParent
                                      // within the loop.

      for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(
               currentParentIndex)) {  // Loop over cell descriptions
        if (p.getSolverNumber() == pChild.getSolverNumber()) {
          subcellPosition.parentIndex = pParent->getParentIndex();
          pParent = &p;
        }
      }
    }
    assertion(pParent->getType() ==
                  exahype::records::ADERDGCellDescription::Ancestor ||
              exahype::records::ADERDGCellDescription::EmptyAncestor);

    // compute subcell index
    for (int xi = 0; xi < DIMENSIONS; ++xi) {
      assertion((pChild.getOffset(xi) - pParent->getOffset(xi)) >= 0);
      subcellPosition.subcellIndex[xi] = tarch::la::round(
          (pChild.getOffset(xi) - pParent->getOffset(xi))/pChild.getSize(xi));
    }
  }

  return subcellPosition;
}


exahype::Cell::SubcellPosition
exahype::Cell::computeSubcellPositionOfDescendant(
    const exahype::records::ADERDGCellDescription& pChild) const {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),toString());
  assertion1(solvers::RegisteredSolvers[pChild.getSolverNumber()]->getType() == exahype::solvers::Solver::Type::ADER_DG,toString());
  assertion1(pChild.getType() == exahype::records::ADERDGCellDescription::Descendant,toString());

  exahype::Cell::SubcellPosition subcellPosition;

  // Initialisation.
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 pChild.getParentIndex()),
             pChild.getParentIndex());
  subcellPosition.parentIndex = pChild.getParentIndex();
  exahype::records::ADERDGCellDescription* pParent = 0;

  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(
           pChild.getParentIndex())) {
    if (p.getSolverNumber() == pChild.getSolverNumber()) {
      pParent = &p;
    }
  }

  if (pParent != 0) {
    // recursion
    while (pParent->getType() == exahype::records::ADERDGCellDescription::EmptyDescendant) {
      const int currentParentIndex =
          pParent->getParentIndex();  // Value must be fixed. We update pParent
                                      // within the loop.
      assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                     currentParentIndex),
                 currentParentIndex);
      for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(
               currentParentIndex)) {  // Loop over cell descriptions
        if (p.getSolverNumber() == pChild.getSolverNumber()) {
          subcellPosition.parentIndex = pParent->getParentIndex();
          pParent = &p;
        }
      }
    }

    assertion(pParent->getType() == exahype::records::ADERDGCellDescription::Descendant ||
              pParent->getType() == exahype::records::ADERDGCellDescription::Cell);

    // compute subcell index
    for (int xi = 0; xi < DIMENSIONS; ++xi) {
      assertion((pChild.getOffset(xi) - pParent->getOffset(xi)) >= 0);
      subcellPosition.subcellIndex[xi] = tarch::la::round(
          (pChild.getOffset(xi) - pParent->getOffset(xi))/pChild.getSize(0));
    }
  } else {
    std::cerr << "exahype::Cell::computeSubcellPositionOfDescendant: parent of "
                 "descendant could not be found!"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  return subcellPosition;
}


int exahype::Cell::getNumberOfADERDGCellDescriptions() const {
  return ADERDGCellDescriptionHeap::getInstance().getData(
      getCellDescriptionsIndex()).size();
}


int exahype::Cell::getNumberOfFiniteVolumeCellDescriptions() const {
  return FiniteVolumesCellDescriptionHeap::getInstance().getData(
      getCellDescriptionsIndex()).size();
}


#ifdef Parallel
void exahype::Cell::clearLoadBalancingWorkloads() {
  if (isRefined()) {
    _cellData.setLocalWorkload(0.0);
    _cellData.setGlobalWorkload(0.0);
  }
  else {
    // @todo really insert here the number of real solvers and weight them accordingly
    _cellData.setLocalWorkload(1.0);
    _cellData.setGlobalWorkload(1.0);
  }
}


void exahype::Cell::restrictLoadBalancingWorkloads(const Cell& childCell, bool isRemote) {
  if (isRemote) {
    // _cellData.setLocalWorkload(  _cellData.getLocalWorkload()  + childCell._cellData.getLocalWorkload() );
    _cellData.setGlobalWorkload(
      std::max(_cellData.getLocalWorkload(), childCell._cellData.getGlobalWorkload())
    );
  }
  else {
    _cellData.setLocalWorkload(  _cellData.getLocalWorkload()  + childCell._cellData.getLocalWorkload() );
    _cellData.setGlobalWorkload( _cellData.getGlobalWorkload() + childCell._cellData.getGlobalWorkload() );
  }

  if ( ADERDGCellDescriptionHeap::getInstance().isValidIndex(getCellDescriptionsIndex()) ) {
    const double embeddedADERDGCells = getNumberOfADERDGCellDescriptions();
    const double embeddedFVPatches   = getNumberOfFiniteVolumeCellDescriptions();

    // @todo this will require further tuning and it might become necessary to
    //       take the order or patch size into account.
    double additionalWeight = 4.0 * embeddedADERDGCells + 1.0 * embeddedFVPatches;

    assertion(additionalWeight>=0.0);

    _cellData.setLocalWorkload(  _cellData.getLocalWorkload()  + additionalWeight );
    _cellData.setGlobalWorkload( _cellData.getGlobalWorkload() + additionalWeight );
  }
}


double exahype::Cell::getLocalWorkload() const {
  return _cellData.getLocalWorkload();
}


double exahype::Cell::getGlobalWorkload() const {
  return _cellData.getGlobalWorkload();
}
#endif


bool exahype::Cell::setSolutionMinMaxAndAnalyseValidity( double* min, double* max, int solverIndex ) {
  assertion( ADERDGCellDescriptionHeap::getInstance().isValidIndex( getCellDescriptionsIndex() ) ) ;
  assertion( ADERDGCellDescriptionHeap::getInstance().getData( getCellDescriptionsIndex() ).size()>static_cast<unsigned int>(solverIndex) ) ;
  assertion( max>=min );

  exahype::records::ADERDGCellDescription& cellDescription = ADERDGCellDescriptionHeap::getInstance().getData(getCellDescriptionsIndex())[solverIndex];

  assertion( exahype::solvers::RegisteredSolvers[ cellDescription.getSolverNumber() ]->getType()==exahype::solvers::Solver::Type::ADER_DG );
  const int numberOfVariables = exahype::solvers::RegisteredSolvers[ cellDescription.getSolverNumber() ]->getNumberOfVariables();

  std::vector<double> cellMinimum( numberOfVariables , std::numeric_limits<double>::max() );
  std::vector<double> cellMaximum( numberOfVariables , std::numeric_limits<double>::min() );

  for (int faceNumber = 0; faceNumber<DIMENSIONS_TIMES_TWO; faceNumber++ )
  for (int i=0; i<numberOfVariables; i++) {
    cellMinimum[i] = std::min(cellMinimum[i], DataHeap::getInstance().getData( cellDescription.getSolutionMin() )[i+faceNumber*numberOfVariables] );
    cellMaximum[i] = std::max(cellMinimum[i], DataHeap::getInstance().getData( cellDescription.getSolutionMax() )[i+faceNumber*numberOfVariables] );
  }

  bool isValidNewCombination = true;
  for (int i=0; i<numberOfVariables; i++) {
    isValidNewCombination &= tarch::la::greater(min[i],cellMinimum[i]);
    isValidNewCombination &= tarch::la::greater(cellMaximum[i],max[i]);
  }

  for (int faceNumber = 0; faceNumber<DIMENSIONS_TIMES_TWO; faceNumber++ )
  for (int i=0; i<numberOfVariables; i++) {
    DataHeap::getInstance().getData( cellDescription.getSolutionMin() )[i+faceNumber*numberOfVariables] = min[i];
    DataHeap::getInstance().getData( cellDescription.getSolutionMax() )[i+faceNumber*numberOfVariables] = max[i];
  }

  return isValidNewCombination;
}


void exahype::Cell::mergeSolutionMinMaxOnFace(
  const int cellDescriptionsIndexOfLeftCell,
  const int cellDescriptionsIndexOfRightCell,
  const int faceIndexForLeftCell,
  const int faceIndexForRightCell
) {
  if (
    ADERDGCellDescriptionHeap::getInstance().isValidIndex( cellDescriptionsIndexOfLeftCell )
    &&
    ADERDGCellDescriptionHeap::getInstance().isValidIndex( cellDescriptionsIndexOfRightCell )
  ) {
    for (auto& leftCellDescription: ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndexOfLeftCell))
    for (auto& rightCellDescription: ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndexOfRightCell)) {
      if (
        (leftCellDescription.getType() == exahype::records::ADERDGCellDescription::Cell ||
        leftCellDescription.getType() == exahype::records::ADERDGCellDescription::Ancestor ||
        leftCellDescription.getType() == exahype::records::ADERDGCellDescription::Descendant)
        &&
        (rightCellDescription.getType() == exahype::records::ADERDGCellDescription::Cell ||
         rightCellDescription.getType() == exahype::records::ADERDGCellDescription::Ancestor ||
         rightCellDescription.getType() == exahype::records::ADERDGCellDescription::Descendant)
      ) {
        assertion( leftCellDescription.getSolverNumber() == rightCellDescription.getSolverNumber() );
        assertion( exahype::solvers::RegisteredSolvers[ leftCellDescription.getSolverNumber() ]->getType()==exahype::solvers::Solver::Type::ADER_DG );
        const int numberOfVariables = exahype::solvers::RegisteredSolvers[ leftCellDescription.getSolverNumber() ]->getNumberOfVariables();
        for (int i=0; i<numberOfVariables; i++) {
          double min = std::min(
            DataHeap::getInstance().getData( leftCellDescription.getSolutionMin() )[i+faceIndexForLeftCell*numberOfVariables],
            DataHeap::getInstance().getData( leftCellDescription.getSolutionMin() )[i+faceIndexForRightCell*numberOfVariables]
          );
          double max = std::min(
            DataHeap::getInstance().getData( leftCellDescription.getSolutionMax() )[i+faceIndexForLeftCell*numberOfVariables],
            DataHeap::getInstance().getData( leftCellDescription.getSolutionMax() )[i+faceIndexForRightCell*numberOfVariables]
	  );

          DataHeap::getInstance().getData( leftCellDescription.getSolutionMin()  )[i+faceIndexForLeftCell*numberOfVariables]  = min;
          DataHeap::getInstance().getData( rightCellDescription.getSolutionMin() )[i+faceIndexForRightCell*numberOfVariables] = min;

          DataHeap::getInstance().getData( leftCellDescription.getSolutionMax()  )[i+faceIndexForLeftCell*numberOfVariables]  = max;
          DataHeap::getInstance().getData( rightCellDescription.getSolutionMax() )[i+faceIndexForRightCell*numberOfVariables] = max;
        }
      }
    }
  }
}


void exahype::Cell::mergeSolutionMinMaxOnFace(
  records::ADERDGCellDescription&  cellDescription,
  int                              faceNumber,
  double* min, double* max
) {
  if (cellDescription.getType() == exahype::records::ADERDGCellDescription::Cell ||
      cellDescription.getType() == exahype::records::ADERDGCellDescription::Ancestor ||
      cellDescription.getType() == exahype::records::ADERDGCellDescription::Descendant
      ) {
    assertion( exahype::solvers::RegisteredSolvers[ cellDescription.getSolverNumber() ]->getType()==exahype::solvers::Solver::Type::ADER_DG );
    const int numberOfVariables = exahype::solvers::RegisteredSolvers[ cellDescription.getSolverNumber() ]->getNumberOfVariables();

    for (int i=0; i<numberOfVariables; i++) {
      DataHeap::getInstance().getData( cellDescription.getSolutionMin()  )[i+faceNumber*numberOfVariables]  = min[i];
      DataHeap::getInstance().getData( cellDescription.getSolutionMax()  )[i+faceNumber*numberOfVariables]  = max[i];
    }
  }
  else {
    assertionMsg(false, "Dominic, please have a look whether we have to do something different here" );
  }
}
