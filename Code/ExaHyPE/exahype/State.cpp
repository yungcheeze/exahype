#include "exahype/State.h"
#include "exahype/Cell.h"
#include "exahype/Vertex.h"

#include "peano/grid/Checkpoint.h"

#include "exahype/solvers/Solver.h"

#include <limits>

exahype::State::State():
Base() {
  _stateData.setCurrentMinTimeStamp(std::numeric_limits<double>::max());

  _stateData.setNextMinTimeStepSize(std::numeric_limits<double>::max());
  _stateData.setCurrentMinTimeStepSize(std::numeric_limits<double>::max());
  _stateData.setPreviousMinTimeStepSize(std::numeric_limits<double>::max());
}

exahype::State::State(const Base::PersistentState& argument):
      Base(argument) {}

double exahype::State::getCurrentMinTimeStepSize() const {
  return _stateData.getCurrentMinTimeStepSize();
}

double exahype::State::getPreviousMinTimeStepSize() const {
  return _stateData.getPreviousMinTimeStepSize();
}

void exahype::State::setCurrentMinTimeStamp(double currentMinTimeStamp) {
  _stateData.setCurrentMinTimeStamp(currentMinTimeStamp);
}

double exahype::State::getCurrentMinTimeStamp() const {
  return _stateData.getCurrentMinTimeStamp();
}

void exahype::State::updateNextMinTimeStepSize(const double timeStepSize) {
  _stateData.setNextMinTimeStepSize( std::min(_stateData.getNextMinTimeStepSize(),timeStepSize) );
}

double exahype::State::getNextMinTimeStepSize() const {
  return _stateData.getNextMinTimeStepSize();
}

void exahype::State::startNewTimeStep() {
  _stateData.setPreviousMinTimeStepSize( _stateData.getCurrentMinTimeStepSize() );

  _stateData.setCurrentMinTimeStepSize ( _stateData.getNextMinTimeStepSize() );
  _stateData.setCurrentMinTimeStamp( _stateData.getCurrentMinTimeStamp()+_stateData.getNextMinTimeStepSize() );

  _stateData.setNextMinTimeStepSize( std::numeric_limits<double>::max() );

  for (
    std::vector<exahype::solvers::Solver*>::const_iterator p = solvers::RegisteredSolvers.begin();
    p != solvers::RegisteredSolvers.end();
    p++
  ) {
    (*p)->startNewTimeStep();
  }
}


void exahype::State::merge(const exahype::State& anotherState) {
  _stateData.setPreviousMinTimeStepSize( std::min(_stateData.getPreviousMinTimeStepSize(),anotherState._stateData.getPreviousMinTimeStepSize()) );
  _stateData.setCurrentMinTimeStepSize ( std::min(_stateData.getCurrentMinTimeStepSize() ,anotherState._stateData.getCurrentMinTimeStepSize ()) );
  _stateData.setCurrentMinTimeStamp( std::min(
      _stateData.getCurrentMinTimeStamp(), anotherState._stateData.getCurrentMinTimeStamp()
  ));
  _stateData.setNextMinTimeStepSize ( std::min(_stateData.getNextMinTimeStepSize() ,anotherState._stateData.getNextMinTimeStepSize ()) );

  // todo 15/02/16:Dominic Etienne Charrier
  // This does only work for consistent registries


  // @todo Fehlerhaft
/*
  int solveNumber = 0;
  for (
      SolveRegistry::iterator p = _solveRegistry.begin();
      p != _solveRegistry.end();
      p++
  ){
    _solveRegistry[ solveNumber ].merge( anotherState._solveRegistry[solveNumber] );
    solveNumber++;

    assertion(_stateData.getNextMinTimeStepSize()     <= (*p).getNextPredictorTimeStepSize());
    assertion(_stateData.getCurrentMinTimeStepSize()  <= (*p).getPredictorTimeStepSize());
    assertion(_stateData.getPreviousMinTimeStepSize() <= (*p).getCorrectorTimeStepSize());
  }
*/
}

void exahype::State::writeToCheckpoint( peano::grid::Checkpoint<exahype::Vertex,exahype::Cell>& checkpoint ) const {
  // do nothing
}

void exahype::State::readFromCheckpoint( const peano::grid::Checkpoint<exahype::Vertex,exahype::Cell>& checkpoint ) {
  // do nothing
}

