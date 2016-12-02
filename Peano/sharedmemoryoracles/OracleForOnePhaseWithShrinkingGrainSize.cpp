#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"
#include "peano/utils/Globals.h"
#include "tarch/Assertions.h"

#include <cstdlib>
#include <limits>
#include <fstream>


tarch::logging::Log  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_log( "sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize" );

const double         sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_InitialRelativeAccuracy(1e-2);
const double         sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_TimingMax( 65536.0 );


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::OracleForOnePhaseWithShrinkingGrainSize(bool learn):
  _activeMethodTrace(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling),
  _learn(learn),
  _measurements() {
}


std::string sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::toString() const {
  std::ostringstream msg;
  msg << "(biggest-problem-size=" << _biggestProblemSize
      << ",current-grain-size=" << _currentGrainSize
      << ",current-measurement=" << _currentMeasurement.toString()
      << ",previous-measured-time=" << _previousMeasuredTime
      << ",search-delta=" << _searchDelta
      << ")";
  return msg.str();
}


peano::datatraversal::autotuning::GrainSize  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelise(int problemSize, peano::datatraversal::autotuning::MethodTrace askingMethod) {
  assertion(askingMethod!=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling);
  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );

  const bool insertNewEntry                          = _measurements.count(askingMethod)==0;
  const bool resetEntryAsProblemOfSizeHasBeenUnknown = _learn && (_measurements.count(askingMethod)==0 || _measurements[askingMethod]._biggestProblemSize<problemSize );

  if ( insertNewEntry || resetEntryAsProblemOfSizeHasBeenUnknown ) {
    _measurements[askingMethod]._biggestProblemSize   = problemSize;
    _measurements[askingMethod]._currentGrainSize     = problemSize;
    _measurements[askingMethod]._currentMeasurement   = tarch::timing::Measurement( 0.0 );
    _measurements[askingMethod]._previousMeasuredTime = _TimingMax;
    _measurements[askingMethod]._searchDelta          = problemSize/2;
    logInfo(
      "parallelise()",
      "inserted new entry for " + toString(askingMethod)
      << ": " << _measurements[askingMethod].toString()
    );
  }

  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );
  assertion( _measurements.count(askingMethod)>0 );

  const bool trackTime = (_activeMethodTrace==askingMethod)
                       && (_measurements[askingMethod]._searchDelta>0);

  logDebug(
    "parallelise()",
    "will track time for " << toString(askingMethod) << "=" << trackTime <<
    " with active trace " << toString(_activeMethodTrace)
  );

  return peano::datatraversal::autotuning::GrainSize(
    _measurements[askingMethod]._currentGrainSize<problemSize ? _measurements[askingMethod]._currentGrainSize : 0,
    trackTime,
    problemSize,
    askingMethod, this
  );
}



void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::changeMeasuredMethodTrace() {
  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );
  _activeMethodTrace = peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling;

  int remainingTriesToFindSearchingTrace = 16;
  while ( _measurements.count(_activeMethodTrace)==0 ) {
    _activeMethodTrace = peano::datatraversal::autotuning::toMethodTrace( rand() % (int)(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling) );
    if (
      remainingTriesToFindSearchingTrace>0
      &&
      _measurements.count(_activeMethodTrace)==1
      &&
      _measurements[_activeMethodTrace]._searchDelta==0
    ) {
      remainingTriesToFindSearchingTrace--;
      logDebug( "changeMeasuredMethodTrace()", "skip " << toString(_activeMethodTrace) << " for measurement " << _measurements[_activeMethodTrace].toString() );
      _activeMethodTrace = peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling;
    }
  }

  assertion(_activeMethodTrace!=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling);
  logDebug( "changeMeasuredMethodTrace()", "next active method trace " << toString(_activeMethodTrace) << " with existing measurement " << _measurements[_activeMethodTrace].toString() );
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::makeAttributesLearn() {
  assertion( _measurements[_activeMethodTrace]._currentMeasurement.isAccurateValue() );
  assertion( _measurements[_activeMethodTrace]._searchDelta>0 );
  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );

  if (
    _measurements[_activeMethodTrace]._currentMeasurement.getValue() < _measurements[_activeMethodTrace]._previousMeasuredTime
  ) {
    logInfo(
      "makeAttributesLearn()",
      "found better scaling parameter choice/serial runtime for " << toString(_activeMethodTrace) << ": " << _measurements[_activeMethodTrace].toString()
    );

    // Previous search has been successful. So we might undershoot
    while ( _measurements[_activeMethodTrace]._currentGrainSize - _measurements[_activeMethodTrace]._searchDelta <= 0 ) {
      _measurements[_activeMethodTrace]._searchDelta /= 2;
    }

    _measurements[_activeMethodTrace]._currentGrainSize     -= _measurements[_activeMethodTrace]._searchDelta;
    _measurements[_activeMethodTrace]._previousMeasuredTime  = _measurements[_activeMethodTrace]._currentMeasurement.getValue();
    _measurements[_activeMethodTrace]._currentMeasurement.erase();
    _measurements[_activeMethodTrace]._currentMeasurement.increaseAccuracy(2.0);

    logInfo(
      "makeAttributesLearn()",
      "continue with " << _measurements[_activeMethodTrace].toString()
    );
  }
  else if (
    _measurements[_activeMethodTrace]._currentMeasurement.getValue() > _measurements[_activeMethodTrace]._previousMeasuredTime
  ) {
    logInfo(
      "makeAttributesLearn()",
      "parameter choice for " << toString(_activeMethodTrace) << " does not scale: " << _measurements[_activeMethodTrace].toString()
    );

    _measurements[_activeMethodTrace]._currentGrainSize     += _measurements[_activeMethodTrace]._searchDelta;
    _measurements[_activeMethodTrace]._previousMeasuredTime  = _TimingMax;
    _measurements[_activeMethodTrace]._currentMeasurement.erase();
    _measurements[_activeMethodTrace]._currentMeasurement.increaseAccuracy(2.0);

    if (_measurements[_activeMethodTrace]._searchDelta<=4) {
      _measurements[_activeMethodTrace]._searchDelta--;
    }
    else {
      _measurements[_activeMethodTrace]._searchDelta /= 2;
    }

    logInfo(
      "makeAttributesLearn()",
      "continue with " << _measurements[_activeMethodTrace].toString() << ". continue to search=" <<
      (_measurements[_activeMethodTrace]._searchDelta>0)
    );
  }

  // @todo Assertion
  if ( _measurements[_activeMethodTrace]._currentGrainSize<=0 ) {
    logInfo(
      "makeAttributesLearn()",
      "error for " << _measurements[_activeMethodTrace].toString() << ". measurement=" <<
      _measurements[_activeMethodTrace].toString()
    );
    exit(-1);
  }

  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelSectionHasTerminated(int problemSize, peano::datatraversal::autotuning::MethodTrace askingMethod, double costPerProblemElement) {
  assertion( askingMethod!=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling );

  if (_measurements[_activeMethodTrace]._currentMeasurement.getAccuracy()==0.0 ) {
    const double computeTime   = costPerProblemElement * static_cast<double>(problemSize);
    _measurements[_activeMethodTrace]._currentMeasurement.setAccuracy( computeTime * _InitialRelativeAccuracy );
    logInfo(
      "parallelSectionHasTerminated(...)",
      "fix accuracy for " << toString(askingMethod) << " to " << _measurements[_activeMethodTrace].toString()
    );
  }

  _measurements[askingMethod]._currentMeasurement.setValue( costPerProblemElement );
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::loadStatistics(const std::string& filename, int oracleNumber) {
  std::ifstream file(filename);
  std::string str = "";

  bool        tagOpen = false;
  while (std::getline(file, str)) {
    tagOpen &= str.compare("end OracleForOnePhaseWithShrinkingGrainSize")!=0;

    if (tagOpen) {
      std::string leftToken = "";
      std::string rightString = str;

      leftToken   = rightString.substr( 0, rightString.find("=") );
      rightString = rightString.substr( leftToken.size()+1 );
      peano::datatraversal::autotuning::MethodTrace  methodTrace = peano::datatraversal::autotuning::toMethodTrace(leftToken);
      logDebug( "loadStatistics(...)", "parse properties for " << toString(methodTrace) );

      leftToken   = rightString.substr( 0, rightString.find(",") );
      rightString = rightString.substr( leftToken.size()+1 );
      int   biggestProblemSize = std::stoi(leftToken);
      logDebug( "loadStatistics(...)", "got biggest problem size " << biggestProblemSize );

      leftToken   = rightString.substr( 0, rightString.find(",") );
      rightString = rightString.substr( leftToken.size()+1 );
      int   currentGrainSize = std::stoi(leftToken);
      logDebug( "loadStatistics(...)", "got current grain size " <<  currentGrainSize );

      leftToken   = rightString.substr( 0, rightString.find(",") );
      rightString = rightString.substr( leftToken.size()+1 );
      double previousMeasuredTime = std::stof(leftToken);
      logDebug( "loadStatistics(...)", "previous measured time is " <<  previousMeasuredTime );

      leftToken   = rightString.substr( 0, rightString.find(",") );
      rightString = rightString.substr( leftToken.size()+1 );
      int   searchDelta = std::stoi(leftToken);

      //leftToken   = rightString.substr( 0, rightString.find("eps=") );
      //rightString = rightString.substr( leftToken.size()+4 );
      //leftToken   = rightString.substr( 0, rightString.find(",") );
      //rightString = rightString.substr( leftToken.size()+1 );
      double accuracy = std::stof(leftToken);
      //double accuracy = 0.0;

      _measurements[methodTrace]._biggestProblemSize    = biggestProblemSize;
      _measurements[methodTrace]._currentGrainSize      = currentGrainSize;
      _measurements[methodTrace]._previousMeasuredTime  = previousMeasuredTime;
      _measurements[methodTrace]._searchDelta           = searchDelta;
      _measurements[methodTrace]._currentMeasurement.erase();
      _measurements[methodTrace]._currentMeasurement.setAccuracy( accuracy );

      logDebug( "loadStatistics(...)", "added " << toString(methodTrace) << ": "<< _measurements[methodTrace].toString() );
    }

    // Older GCC versions require an explicit cast here
    tagOpen |= str.compare( "adapter-number=" + std::to_string( (long long)oracleNumber) )==0;
  }

  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::plotStatistics(std::ostream& out, int oracleNumber) const {
  out << "begin OracleForOnePhaseWithShrinkingGrainSize" << std::endl;
  out << "adapter-number=" << oracleNumber << std::endl;

  for (auto measurement: _measurements) {
    out << peano::datatraversal::autotuning::toString(measurement.first)
        << "=" << measurement.second._biggestProblemSize
        << "," << measurement.second._currentGrainSize
        << "," << measurement.second._previousMeasuredTime
        << "," << measurement.second._searchDelta
        << "," << measurement.second._currentMeasurement.getAccuracy()
        << std::endl;
  }

  out << "end OracleForOnePhaseWithShrinkingGrainSize" << std::endl;
}


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::~OracleForOnePhaseWithShrinkingGrainSize() {
}


peano::datatraversal::autotuning::OracleForOnePhase* sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::createNewOracle() const {
  return new OracleForOnePhaseWithShrinkingGrainSize(_learn);
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::deactivateOracle() {
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::activateOracle() {
  if (_learn) {
    //
    // First check is very important. If we skip it, the second statement would
    // insert elements into map
    //
    if (
      _activeMethodTrace!=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling
      &&
      _measurements[_activeMethodTrace]._currentMeasurement.isAccurateValue()
      &&
      _measurements[_activeMethodTrace]._searchDelta>0
    ) {
      makeAttributesLearn();
    }

    if (!_measurements.empty()) {
      changeMeasuredMethodTrace();
    }
  }
}
