// This file originally was created by pdt (Peano Development Toolkit) as part 
// of a code based upon the Peano project by Tobias Weinzierl. For conditions 
// of distribution and use of this project, please see the copyright notice at
// www.peano-framework.org. Feel free to adopt the license and authorship of 
// this file and your project to your needs as long as the license is in 
// agreement with the original Peano user constraints. A reference to/citation 
// of  Peano and its author is highly appreciated.
#ifndef _EXAHYPE_CELL_H_ 
#define _EXAHYPE_CELL_H_

#include "exahype/records/Cell.h"
#include "peano/grid/Cell.h"

// ! Begin of code for multiscalelinkedcell toolbox.
#include "peano/heap/Heap.h"
#include "exahype/records/ADERDGCellDescription.h"
// ! End of code for multiscalelinkedcell toolbox.

namespace exahype {
  class Cell;
  // ! Begin of code for multiscalelinkedcell toolbox..
  typedef peano::heap::PlainHeap< exahype::records::ADERDGCellDescription > ADERDGCellDescriptionHeap;
  typedef peano::heap::PlainDoubleHeap DataHeap;
  // ! End of code for multiscalelinkedcell toolbox.
}


/**
 * Blueprint for cell.
 * 
 * This file has originally been created by the PDT and may be manually extended to 
 * the needs of your application. We do not recommend to remove anything!
 */
class exahype::Cell: public peano::grid::Cell< exahype::records::Cell > { 
private:
  typedef class peano::grid::Cell< exahype::records::Cell >  Base;

public:
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

  // ! Begin of code for multiscalelinkedcell toolbox
  int getADERDGCellDescriptionsIndex() const;

  inline exahype::records::ADERDGCellDescription& getADERDGCellDescription(int pdeIndex) {
    return ADERDGCellDescriptionHeap::getInstance().getData(getADERDGCellDescriptionsIndex())[pdeIndex];
  }

  void initCellWithDefaultValues();

  void initCellInComputeTree(
      const int level,
      const tarch::la::Vector<DIMENSIONS,double> size,
      const int numberOfPDEs,
      const int order,
      const int numberOfVariables);
  // ! End of code for multiscalelinkedcell toolbox
};


#endif
