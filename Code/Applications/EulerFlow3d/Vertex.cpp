#include "EulerFlow3d/Vertex.h"
#include "peano/utils/Loop.h"
#include "peano/grid/Checkpoint.h"

#include "EulerFlow3d/multiscalelinkedcell/HangingVertexBookkeeper.h"

exahype::Vertex::Vertex():
  Base() { 
  _vertexData._persistentRecords._cellDescriptionsIndex =
      multiscalelinkedcell::HangingVertexBookkeeper::getInstance().createVertexLinkMapForNewVertex();
}


exahype::Vertex::Vertex(const Base::DoNotCallStandardConstructor& value):
  Base(value) { 
  // do nothing
}


exahype::Vertex::Vertex(const Base::PersistentVertex& argument):
  Base(argument) {
  // do nothing
}

tarch::la::Vector<TWO_POWER_D,int>&
exahype::Vertex::getCellDescriptionsIndex() {
  return _vertexData._persistentRecords._cellDescriptionsIndex;
}
