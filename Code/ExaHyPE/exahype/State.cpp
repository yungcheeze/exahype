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


void exahype::State::startNewTimeStep() {
  _stateData.setOldMaxTimeStepSize(_stateData.getMaxTimeStepSize());
  // @todo Dominic Etienne Charrier
  // Commented out due to possible side effects
  //_stateData.setMinTimeStamp(std::numeric_limits<double>::max());
}


void exahype::State::resetAccumulatedValues() {
  _stateData.setMaxTimeStepSize(std::numeric_limits<double>::max());
}


void exahype::State::updateMaxTimeStepSize(const double timeStepSize) {
  _stateData.setMaxTimeStepSize( std::min(_stateData.getMaxTimeStepSize(),timeStepSize) );
}


void exahype::State::updateTimeStamp(const double timeStamp) {
  _stateData.setMinTimeStamp( std::min(_stateData.getMinTimeStamp(),timeStamp) );
}

// \todo Dominic Etienne Charrier
// I changed the definition of this file.
double exahype::State::getMaxTimeStepSize() const {
  return _stateData.getMaxTimeStepSize();
}

double exahype::State::getOldMaxTimeStepSize() const {
  return _stateData.getOldMaxTimeStepSize();
}

// \todo Dominic Etienne Charrier
// I changed the definition of this file.
void exahype::State::merge(const exahype::State& anotherState) {
  _stateData.setMaxTimeStepSize(std::min(_stateData.getMaxTimeStepSize(),anotherState._stateData.getMaxTimeStepSize()));
  _stateData.setMinTimeStamp(std::min(
      _stateData.getMinTimeStamp()+_stateData.getMaxTimeStepSize(),
      anotherState._stateData.getMinTimeStamp()+anotherState._stateData.getMaxTimeStepSize()
  ));
}


double exahype::State::getMinimalGlobalTimeStamp() const {
  return _stateData.getMinTimeStamp();
}
