#include "EulerFlow3d/State.h"
#include "EulerFlow3d/Cell.h"
#include "EulerFlow3d/Vertex.h"

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

void exahype::State::setTimeStepSize(const double timeStepSize) {
  _stateData.setTimeStepSize(timeStepSize);
}

double exahype::State::getTimeStepSize() const {
  return _stateData.getTimeStepSize();
}

void exahype::State::setMinimumTimeStepSizeOfBoth(exahype::State anotherState) {
  _stateData.setTimeStepSize(std::min(_stateData.getTimeStepSize(),anotherState.getTimeStepSize()));
}
