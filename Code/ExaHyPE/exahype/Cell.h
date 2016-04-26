// This file originally was created by pdt (Peano Development Toolkit) as part
// of a code based upon the Peano project by Tobias Weinzierl. For conditions
// of distribution and use of this project, please see the copyright notice at
// www.peano-framework.org. Feel free to adopt the license and authorship of
// this file and your project to your needs as long as the license is in
// agreement with the original Peano user constraints. A reference to/citation
// of  Peano and its author is highly appreciated.
#ifndef _EXAHYPE_CELL_H_
#define _EXAHYPE_CELL_H_

#include "exahype/State.h"

#include "exahype/records/Cell.h"
#include "peano/grid/Cell.h"

#include "peano/heap/DoubleHeap.h"
#include "exahype/records/ADERDGCellDescription.h"

namespace exahype {
  class Cell;

  typedef peano::heap::PlainHeap<exahype::records::ADERDGCellDescription>
  ADERDGCellDescriptionHeap;
  typedef peano::heap::PlainDoubleHeap DataHeap;
}



/**
 * @todo Dominic We should add some proper descriptions here one day.
 */
class exahype::Cell : public peano::grid::Cell<exahype::records::Cell> {
private:
  typedef class peano::grid::Cell<exahype::records::Cell> Base;

  static tarch::logging::Log _log;

public:
  typedef struct {
    int parentIndex;
    tarch::la::Vector<DIMENSIONS,int> subcellIndex;
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
   * Returns the heap index of the first ADERDGCellDescription associated
   * with this cell.
   */
  int getADERDGCellDescriptionsIndex() const;

  /**
   * Loads the ADERDGCellDescription associated
   * with this cell.
   */
  inline exahype::records::ADERDGCellDescription& getADERDGCellDescription(
      int index) {
    return ADERDGCellDescriptionHeap::getInstance().getData(
        getADERDGCellDescriptionsIndex())[index];
  }

  /**
   * Per existing cell, the configuration runs over all solvers that
   * shall be realised and links the associated cell descriptions
   * to their parents.
   */
  void addNewCellDescription(
    const int                                            solverNumber,
    const exahype::records::ADERDGCellDescription::Type  cellType,
    const exahype::records::ADERDGCellDescription::RefinementEvent refinementEvent,
    const exahype::records::ADERDGCellDescription::RefinementEvent augmentationEvent,
    const int                                            level,
    const int                                            parentIndex,
    const tarch::la::Vector<DIMENSIONS, int>&            fineGridPositionOfCell,
    const tarch::la::Vector<DIMENSIONS, double>&         size,
    const tarch::la::Vector<DIMENSIONS, double>&         cellCentre
  );

  /**
   * Per existing cell, the initialisation runs over all solvers that
   * shall be realised and allocates memory for the associated.
   * cell descriptions.
   */
  void cleanCellDescription(const int solverNumber);

  /**
   * Per existing cell, the initialisation runs over all solvers
   * and frees unused heap storage that was allocated for the cell descriptions.
   */
  void initialiseCellDescription(const int solverNumber);

//  /**
//   * Per existing cell, the initialisation has to run over all solvers that
//   * shall be realised. For a given solver/PDE a patch is to be created if
//   * the cell's level equals getMinimumTreeDepth() of if the tree is
//   * unrefined.
//   */
//  void init(const int level, const tarch::la::Vector<DIMENSIONS, double>& size,
//            const tarch::la::Vector<DIMENSIONS, double>& cellCentre);

  /**
   * Determine the position of a virtual shell with respect
   * to a  real parent cell or a virtual parent shell that contains data, i.e.,
   * has at least one neighbour that is a real cell.
   *
   * This method is required for the face data prolongation, the
   * volume data prolongation, and the FV volume data prolongation.
   *
   * @todo:16/04/09:Dominic Etienne Charrier: I am not sure if this the
   * right file/class to hold this functionality.
   */
  SubcellPosition getSubcellPositionOfVirtualShell(
      const exahype::records::ADERDGCellDescription& pChild) const;
};

#endif
