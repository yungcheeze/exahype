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
 
#include "exahype/Vertex.h"
#include "peano/utils/Loop.h"
#include "peano/grid/Checkpoint.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

exahype::Vertex::Vertex() : Base() {
  _vertexData._persistentRecords._CellDescriptionsIndex =
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
exahype::Vertex::getCellDescriptionsIndex() {
  return _vertexData._persistentRecords._CellDescriptionsIndex;
}

bool exahype::Vertex::hasToSendMetadata(
    const exahype::State* state,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const int toRank) {
  const int srcScalar  = peano::utils::dLinearisedWithoutLookup(src,2);
  const int destScalar = peano::utils::dLinearisedWithoutLookup(dest,2);
  const tarch::la::Vector<TWO_POWER_D,int> adjacentRanks = getAdjacentRanks();

  return tarch::la::countEqualEntries(dest, src)  == 1 &&
         adjacentRanks(destScalar)   == toRank &&
         (adjacentRanks(srcScalar)   == tarch::parallel::Node::getInstance().getRank() ||
         state->isForkTriggeredForRank(adjacentRanks(srcScalar)));
}

bool exahype::Vertex::hasToReceiveMetadata(
    const exahype::State* state,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const int fromRank) {
  const int srcScalar  = peano::utils::dLinearisedWithoutLookup(src,2);
  const int destScalar = peano::utils::dLinearisedWithoutLookup(dest,2);
  const tarch::la::Vector<TWO_POWER_D,int> adjacentRanks = getAdjacentRanks();

  return tarch::la::countEqualEntries(dest, src) == 1 &&
      adjacentRanks(srcScalar)    == fromRank &&
      (adjacentRanks(destScalar)  == tarch::parallel::Node::getInstance().getRank() ||
       state->isForkingRank(adjacentRanks(destScalar)));
}
