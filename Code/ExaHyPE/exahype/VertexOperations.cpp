// Do not modify any part of this file. It will be overwritten throughout the
// next pdt run.

#include "exahype/VertexOperations.h"
#include "peano/utils/Loop.h"
#include "peano/grid/Checkpoint.h"

exahype::VertexOperations::VertexOperations() {}

tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D, int>
exahype::VertexOperations::readADERDGCellDescriptionsIndex(
    const peano::grid::VertexEnumerator& enumerator,
    const Vertex* const vertices) {
  tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D, int> result;
  dfor2(x) tarch::la::slice(
      result,
      vertices[enumerator(x)]._vertexData.getADERDGCellDescriptionsIndex(),
      xScalar * TWO_POWER_D);
  enddforx return result;
}

tarch::la::Vector<TWO_POWER_D, int>
exahype::VertexOperations::readADERDGCellDescriptionsIndex(
    const Vertex& vertex) {
  return vertex._vertexData.getADERDGCellDescriptionsIndex();
}

void exahype::VertexOperations::writeADERDGCellDescriptionsIndex(
    const peano::grid::VertexEnumerator& enumerator, Vertex* const vertices,
    const tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D, int>& values) {
  dfor2(x) tarch::la::Vector<TWO_POWER_D, int> temp =
      tarch::la::slice<TWO_POWER_D>(values, xScalar * TWO_POWER_D);
  vertices[enumerator(x)]._vertexData.setADERDGCellDescriptionsIndex(temp);
  enddforx
}

void exahype::VertexOperations::writeADERDGCellDescriptionsIndex(
    Vertex& vertex, const tarch::la::Vector<TWO_POWER_D, int>& values) {
  vertex._vertexData.setADERDGCellDescriptionsIndex(values);
}
