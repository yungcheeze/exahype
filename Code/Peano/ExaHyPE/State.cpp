#include "ExaHyPE/State.h"
#include "ExaHyPE/Cell.h"
#include "ExaHyPE/Vertex.h"

#include "peano/grid/Checkpoint.h"



ExaHyPE::State::State():
  Base() { 
  // @todo Insert your code here
}


ExaHyPE::State::State(const Base::PersistentState& argument):
  Base(argument) {
  // @todo Insert your code here
}


void ExaHyPE::State::writeToCheckpoint( peano::grid::Checkpoint<ExaHyPE::Vertex,ExaHyPE::Cell>& checkpoint ) const {
  // @todo Insert your code here
}

    
void ExaHyPE::State::readFromCheckpoint( const peano::grid::Checkpoint<ExaHyPE::Vertex,ExaHyPE::Cell>& checkpoint ) {
  // @todo Insert your code here
}
