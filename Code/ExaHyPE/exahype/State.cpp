#include "exahype/State.h"
#include "exahype/Cell.h"
#include "exahype/Vertex.h"

#include "peano/grid/Checkpoint.h"

#include "exahype/solvers/Solver.h"

#include <limits>

exahype::State::State():
Base() {
}

exahype::State::State(const Base::PersistentState& argument):
      Base(argument) {
  // do nothing
}

void exahype::State::merge(const exahype::State& anotherState) {
  // do nothing
}

void exahype::State::writeToCheckpoint( peano::grid::Checkpoint<exahype::Vertex,exahype::Cell>& checkpoint ) const {
  // do nothing
}

void exahype::State::readFromCheckpoint( const peano::grid::Checkpoint<exahype::Vertex,exahype::Cell>& checkpoint ) {
  // do nothing
}

