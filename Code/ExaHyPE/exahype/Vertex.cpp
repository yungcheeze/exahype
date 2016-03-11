#include "exahype/Vertex.h"
#include "peano/utils/Loop.h"
#include "peano/grid/Checkpoint.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

exahype::Vertex::Vertex() : Base() {
  _vertexData._persistentRecords._ADERDGCellDescriptionsIndex =
      multiscalelinkedcell::HangingVertexBookkeeper::getInstance()
          .createVertexLinkMapForNewVertex();
}

exahype::Vertex::Vertex(const Base::DoNotCallStandardConstructor& value)
    : Base(value) {
  // do nothing
}

exahype::Vertex::Vertex(const Base::PersistentVertex& argument)
    : Base(argument) {
  // do nothing
}

tarch::la::Vector<TWO_POWER_D, int>&
exahype::Vertex::getADERDGCellDescriptionsIndex() {
  return _vertexData._persistentRecords._ADERDGCellDescriptionsIndex;
}
