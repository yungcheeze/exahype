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
 
#ifndef _EXAHYPE_CELL_H_
#define _EXAHYPE_CELL_H_

#include "exahype/State.h"

#include "exahype/records/Cell.h"
#include "peano/grid/Cell.h"
#include "peano/grid/VertexEnumerator.h"

#include "exahype/records/ADERDGCellDescription.h"
#include "exahype/records/FiniteVolumesCellDescription.h"
#include "peano/heap/DoubleHeap.h"

namespace exahype {
  class Cell;

  typedef peano::heap::PlainHeap<exahype::records::ADERDGCellDescription>         ADERDGCellDescriptionHeap;
  typedef peano::heap::PlainHeap<exahype::records::FiniteVolumesCellDescription>  FiniteVolumesCellDescriptionHeap;
  typedef peano::heap::PlainDoubleHeap DataHeap;
}

/**
 * @todo Dominic We should add some proper descriptions here one day.
 */
class exahype::Cell : public peano::grid::Cell<exahype::records::Cell> {
 private:
  typedef class peano::grid::Cell<exahype::records::Cell> Base;

  static tarch::logging::Log _log;

  /**
   * Each cell points to a series of cell descriptions. The array holding the
   * series has to be stored on the heap and, consequently, initialised
   * properly. This is done by create() while destroy() cleans up. Please note
   * that you have to invoke create() once before you do anything with the cell
   * at all. You should destroy() in return in the very end.
   *
   * The operation shows that each cell in the tree can theoretically hold a
   * solver though only few do.
   *
   * This operation is used by addNewCellDescription().
   */
  void setupMetaData();

 public:
  typedef struct {
    int parentIndex;
    tarch::la::Vector<DIMENSIONS, int> subcellIndex;
  } SubcellPosition;

  /**
   * Default Constructor
   *
   * This constructor is required by the framework's data container. Do not
   * remove it.
   */
  Cell();

  /**
   * This constructor should not set any attributes. It is used by the
   * traversal algorithm whenever it allocates an array whose elements
   * will be overwritten later anyway.
   */
  Cell(const Base::DoNotCallStandardConstructor&);

  /**
   * Constructor
   *
   * This constructor is required by the framework's data container. Do not
   * remove it. It is kind of a copy constructor that converts an object which
   * comprises solely persistent attributes into a full attribute. This very
   * functionality is implemented within the super type, i.e. this constructor
   * has to invoke the correponsing super type's constructor and not the super
   * type standard constructor.
   */
  Cell(const Base::PersistentCell& argument);

  /**
   * Allocates a vector for the
   */
  void initialiseStorageOnHeap();

  /**
   * Returns meta data describing the surrounding cell descriptions. The
   * routine is notably used by the automated adapters to derive adjacency
   * information on the cell level.
   */
  int getCellDescriptionsIndex() const;

  int getNumberOfADERDGCellDescriptions() const;
  int getNumberOfFiniteVolumeCellDescriptions() const;

  /**
   * Loads a ADERDGCellDescription associated
   * with this cell.
   */
  inline exahype::records::ADERDGCellDescription& getADERDGCellDescription(
      int index) {
    return ADERDGCellDescriptionHeap::getInstance().getData(
        getCellDescriptionsIndex())[index];
  }

  /**
   * Loads a ADERDGCellDescription associated
   * with this cell.
   */
  inline exahype::records::FiniteVolumesCellDescription& getFiniteVolumesCellDescription(
      int index) {
    return FiniteVolumesCellDescriptionHeap::getInstance().getData(
        getCellDescriptionsIndex())[index];
  }

  /**
   * @see setupMetaData()
   */
  void shutdownMetaData();

  /**
   * todo docu
   *
   * setupMetaData() is called if cell hasn't been properly initialised before.
   */
  void addNewCellDescription(
      const int solverNumber,
      const exahype::records::ADERDGCellDescription::Type cellType,
      const exahype::records::ADERDGCellDescription::RefinementEvent
          refinementEvent,
      const int level, const int parentIndex,
      const tarch::la::Vector<DIMENSIONS, double>& size,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre);

  void addNewCellDescription(
      const int solverNumber,
      const exahype::records::FiniteVolumesCellDescription::Type cellType,
//      const exahype::records::FiniteVolumesCellDescription::RefinementEvent
//          refinementEvent,    @todo Dominic
      const int level, const int parentIndex,
      const tarch::la::Vector<DIMENSIONS, double>& size,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre);

  /**
   * Checks if no unnecessary  memory is allocated for the cell
   * description the solverNumber is referring to.
   * If this is the case, it deallocates the unnecessarily allocated memory.
   */
  void ensureNoUnnecessaryMemoryIsAllocated(const int solverNumber);

  /**
   * Checks if all the necessary memory is allocated for the cell
   * description the solverNumber is referring to.
   * If this is not the case, it allocates the necessary
   * memory for the cell description.
   */
  void ensureNecessaryMemoryIsAllocated(const int solverNumber);

  /**
   * Determine the position of a Cell or Ancestor with respect
   * to a parent of type Ancestor.
   *
   * This method is required for the face data restriction, the
   * volume data restriction, and the FV volume data restriction.
   *
   * @todo:16/04/09:Dominic Etienne Charrier: I am not sure if this the
   * right file/class to hold this functionality.
   */
  SubcellPosition computeSubcellPositionOfCellOrAncestor(
      const exahype::records::ADERDGCellDescription& pChild) const;

  /**
   * Determine the position of a Descendant with respect
   * to a  Cell or Descendant that contains data, i.e.,
   * has at least one neighbour that is a real cell.
   *
   * This method is required for the face data prolongation, the
   * volume data prolongation, and the FV volume data prolongation.
   *
   * @todo:16/04/09:Dominic Etienne Charrier: I am not sure if this the
   * right file/class to hold this functionality.
   */
  SubcellPosition computeSubcellPositionOfDescendant(
      const exahype::records::ADERDGCellDescription& pChild) const;

  bool isInitialised() const;

  /**
   * This operation runs a dozen of assertions. It becomes nop if assertions
   * are switched off.
   */
  void validateNoNansInADERDGSolver(
    int                                  number,
    exahype::Cell&                       fineGridCell,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const std::string&                   methodTraceOfCaller
  );

  #ifdef Parallel
  void clearLoadBalancingWorkloads();
  void restrictLoadBalancingWorkloads(const Cell& childCell, bool isRemote);
  double getLocalWorkload() const;
  double getGlobalWorkload() const;
  #endif
};

#endif
