#include "EulerFlow3d/State.h"
#include "EulerFlow3d/Cell.h"
#include "EulerFlow3d/Vertex.h"

#include "peano/grid/Checkpoint.h"



exahype::State::State():
  Base() { 
  _stateData.setTimeStepSize(1e20);
}


exahype::State::State(const Base::PersistentState& argument):
  Base(argument) {
  // do nothing
}

void exahype::State::writeToCheckpoint( peano::grid::Checkpoint<exahype::Vertex,exahype::Cell>& checkpoint ) const {
  // do nothing
}

    
void exahype::State::readFromCheckpoint( const peano::grid::Checkpoint<exahype::Vertex,exahype::Cell>& checkpoint ) {
  // do nothing
}

void exahype::State::setTimeStepSize(const double timeStepSize) {
  _stateData.setTimeStepSize(timeStepSize);
}

double exahype::State::getTimeStepSize() const {
  return _stateData.getTimeStepSize();
}

void exahype::State::setMinimumTimeStepSizeOfBoth(exahype::State anotherState) {
  _stateData.setTimeStepSize(std::min(_stateData.getTimeStepSize(),anotherState.getTimeStepSize()));
}
