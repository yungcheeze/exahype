#include "exahype/Cell.h"


exahype::Cell::Cell():
  Base() { 
  // @todo Insert your code here
}


exahype::Cell::Cell(const Base::DoNotCallStandardConstructor& value):
  Base(value) { 
  // Please do not insert anything here
}

exahype::Cell::Cell(const Base::PersistentCell& argument):
  Base(argument) {
  // @todo Insert your code here
}
