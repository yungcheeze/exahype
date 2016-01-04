#include "EulerFlow/State.h"
#include "EulerFlow/Cell.h"
#include "EulerFlow/Vertex.h"

#include "peano/grid/Checkpoint.h"

#include "limits"

exahype::State::State():
  Base() { 
  _stateData.setTimeStepSize(std::numeric_limits<double>::max());
  _stateData.setOldTimeStepSize(std::numeric_limits<double>::max());
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

void exahype::State::setTimeStepSizeToMaxValue() {
  _stateData.setTimeStepSize(std::numeric_limits<double>::max());
}

void exahype::State::setTimeStepSize(const double newTimeStepSize) {
  _stateData.setTimeStepSize(newTimeStepSize);
}

double exahype::State::getTimeStepSize() const {
  return _stateData.getTimeStepSize();
}

void exahype::State::setOldTimeStepSize(const double newOldTimeStepSize) {
  _stateData.setOldTimeStepSize(newOldTimeStepSize);
}

double exahype::State::getOldTimeStepSize() const {
  return _stateData.getOldTimeStepSize();
}


void exahype::State::setMinimumTimeStepSizeOfBoth(exahype::State anotherState) {
  _stateData.setTimeStepSize(std::min(_stateData.getTimeStepSize(),anotherState.getTimeStepSize()));
}
