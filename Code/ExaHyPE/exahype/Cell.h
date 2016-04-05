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

// ! Begin of code for multiscalelinkedcell toolbox.
#include "peano/heap/DoubleHeap.h"
#include "exahype/records/ADERDGCellDescription.h"
// ! End of code for multiscalelinkedcell toolbox.

namespace exahype {
  class Cell;
  // ! Begin of code for multiscalelinkedcell toolbox..
  typedef peano::heap::PlainHeap<exahype::records::ADERDGCellDescription>
  ADERDGCellDescriptionHeap;
  typedef peano::heap::PlainDoubleHeap DataHeap;
  // ! End of code for multiscalelinkedcell toolbox.
}

/**
 * Blueprint for cell.
 *
 * This file has originally been created by the PDT and may be manually extended
 *to
 * the needs of your application. We do not recommend to remove anything!
 */
class exahype::Cell : public peano::grid::Cell<exahype::records::Cell> {
private:
  typedef class peano::grid::Cell<exahype::records::Cell> Base;

  static tarch::logging::Log _log;

public:
  /**
   * Type of a cell description.
   * Cell descriptions of type \p Cell hold cell and face data,
   * while the ones of type \p Shell hold only face data.
   * Both belong to the original spacetree that
   * is constructed according to solver-based refinement criteria.
   * Virtual shells hold also only face data but belong to
   * the virtual part of the augmented spacetree that
   * is created to store prolongated face data.
   */
  enum CellDescriptionType {
    Unspecified,
    RealCell,
    RealShell,
    VirtualShell
  };

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
   * with this cell and the solver with index \p solverIndex.
   */
  inline exahype::records::ADERDGCellDescription& getADERDGCellDescription(
      int solverIndex) {
    return ADERDGCellDescriptionHeap::getInstance().getData(
        getADERDGCellDescriptionsIndex())[solverIndex];
  }

  /**
   * Per existing cell, the configuration runs over all solvers that
   * shall be realised and links the associated cell descriptions
   * to their parents.
   */
  void addNewCellDescription(
      const int solverNumber,
      const exahype::Cell::CellDescriptionType cellType,
      const int level,
      const int parentIndex,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const tarch::la::Vector<DIMENSIONS, double>& size,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre);

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

  /**
   * Per existing cell, the initialisation has to run over all solvers that
   * shall be realised. For a given solver/PDE a patch is to be created if
   * the cell's level equals getMinimumTreeDepth() of if the tree is
   * unrefined.
   */
  void init(const int level, const tarch::la::Vector<DIMENSIONS, double>& size,
            const tarch::la::Vector<DIMENSIONS, double>& cellCentre);
};

#endif
