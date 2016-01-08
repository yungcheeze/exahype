#include "exahype/State.h"
#include "exahype/Cell.h"
#include "exahype/Vertex.h"

#include "peano/grid/Checkpoint.h"



exahype::State::State():
  Base() { 
  // @todo Insert your code here
}


exahype::State::State(const Base::PersistentState& argument):
  Base(argument) {
  // @todo Insert your code here
}


void exahype::State::writeToCheckpoint( peano::grid::Checkpoint<exahype::Vertex,exahype::Cell>& checkpoint ) const {
  // @todo Insert your code here
}

    
void exahype::State::readFromCheckpoint( const peano::grid::Checkpoint<exahype::Vertex,exahype::Cell>& checkpoint ) {
  // @todo Insert your code here
}
