#include "ExaHyPeDataStructure/Vertex.h"
#include "peano/utils/Loop.h"
#include "peano/grid/Checkpoint.h"


exahype::Vertex::Vertex():
  Base() { 
  // @todo Insert your code here
}


exahype::Vertex::Vertex(const Base::DoNotCallStandardConstructor& value):
  Base(value) { 
  // Please do not insert anything here
}


exahype::Vertex::Vertex(const Base::PersistentVertex& argument):
  Base(argument) {
  // @todo Insert your code here
}

tarch::la::Vector<TWO_POWER_D,int>
exahype::Vertex::getCellDescriptionIndex() {
  return _vertexData.getCellDescriptionIndex();
}
