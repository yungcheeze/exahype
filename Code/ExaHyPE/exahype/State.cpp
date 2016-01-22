#include "exahype/State.h"
#include "exahype/Cell.h"
#include "exahype/Vertex.h"

#include "peano/grid/Checkpoint.h"

#include "limits"

exahype::State::State():
  Base() { 
  _stateData.setMaxTimeStepSize(std::numeric_limits<double>::max());
  _stateData.setOldMaxTimeStepSize(std::numeric_limits<double>::max());
}


exahype::State::State(const Base::PersistentState& argument):
  Base(argument) {
}


void exahype::State::writeToCheckpoint( peano::grid::Checkpoint<exahype::Vertex,exahype::Cell>& checkpoint ) const {
  // do nothing
}

    
void exahype::State::readFromCheckpoint( const peano::grid::Checkpoint<exahype::Vertex,exahype::Cell>& checkpoint ) {
  // do nothing
}


void exahype::State::resetAccumulatedValues() {
  _stateData.setOldMaxTimeStepSize(_stateData.getMaxTimeStepSize());

  _stateData.setMaxTimeStepSize(std::numeric_limits<double>::max());
}


void exahype::State::updateMaxTimeStepSize(const double timeStepSize) {
  _stateData.setMaxTimeStepSize( std::min(_stateData.getMaxTimeStepSize(),timeStepSize) );
}


double exahype::State::getMaxTimeStepSize() const {
  return _stateData.getOldMaxTimeStepSize();
}


void exahype::State::merge(const exahype::State& anotherState) {
  _stateData.setMaxTimeStepSize(std::min(_stateData.getMaxTimeStepSize(),anotherState._stateData.getMaxTimeStepSize()));
}


double exahype::State::getMinimalGlobalTimeStamp() const {
  return 1.0; //_stateData.getTime();
}
