#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"
#include "peano/utils/Globals.h"
#include "tarch/Assertions.h"
#include "tarch/multicore/Core.h"


#include <cstdlib>
#include <limits>
#include <fstream>
#include <stdexcept>


tarch::logging::Log  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_log( "sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize" );

const double         sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_InitialRelativeAccuracy(1e-2);
const double         sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_TimingMax( 65536.0 );



sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::OracleForOnePhaseWithShrinkingGrainSize(bool learn):
  _activeMethodTrace(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling),
  _learn(learn),
  _measurements(),
  _totalRuntimeMeasurement(),
  _totalRuntimeWatch( "sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize", "OracleForOnePhaseWithShrinkingGrainSize", false ) {
}


std::string sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::toString() const {
  std::ostringstream msg;
  msg << "(biggest-problem-size=" << _biggestProblemSize
      << ",current-grain-size=" << _currentGrainSize
      << ",current-measurement=" << _currentMeasurement.toString()
      << ",previous-measured-time=" << _previousMeasuredTime
      << ",search-delta=" << _searchDelta
      << ",total-serial-time=" << _totalSerialTime;

  if (_currentMeasurement.getValue()>1e-12) {
    msg << ",estimated-speedup=" << _serialMeasurement/_currentMeasurement.getValue();
  }
  else if (_previousMeasuredTime<_TimingMax) {
    msg << ",estimated-speedup=" << _serialMeasurement/_previousMeasuredTime;
  }
  else {
    msg << ",estimated-speedup=not-available";
  }

  msg << ")";
  return msg.str();
}


peano::datatraversal::autotuning::GrainSize  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelise(int problemSize, peano::datatraversal::autotuning::MethodTrace askingMethod) {
  assertion( askingMethod != peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling );

  assertion(askingMethod!=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling);
  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );

  const bool unknownProblemSizeForScalingSetup       = _learn && _measurements.count(askingMethod)==1
                                                    && _measurements[askingMethod]._biggestProblemSize<problemSize
                                                    && _measurements[askingMethod]._currentGrainSize<_measurements[askingMethod]._biggestProblemSize;
  const bool insertNewEntry                          = _measurements.count(askingMethod)==0;
  const bool resetEntryAsProblemOfSizeHasBeenUnknown = _learn && (_measurements.count(askingMethod)==0 || _measurements[askingMethod]._biggestProblemSize<problemSize );

  if (unknownProblemSizeForScalingSetup) {
    _measurements[askingMethod]._biggestProblemSize   = problemSize;
    _measurements[askingMethod]._currentMeasurement   = tarch::timing::Measurement( 0.0 );
    _measurements[askingMethod]._previousMeasuredTime = _TimingMax;
  }
  else if ( insertNewEntry || resetEntryAsProblemOfSizeHasBeenUnknown ) {
    _measurements[askingMethod]._biggestProblemSize   = problemSize;
    _measurements[askingMethod]._currentGrainSize     = problemSize;
    _measurements[askingMethod]._currentMeasurement   = tarch::timing::Measurement( 0.0 );
    _measurements[askingMethod]._previousMeasuredTime = _TimingMax;
    _measurements[askingMethod]._searchDelta          = problemSize < tarch::multicore::Core::getInstance().getNumberOfThreads()*2 ? problemSize/2 : problemSize - problemSize / tarch::multicore::Core::getInstance().getNumberOfThreads();
    _measurements[askingMethod]._serialMeasurement    = 0.0;
    _measurements[askingMethod]._totalSerialTime      = 0.0;

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

  const auto oldActiveMethodTrace = _activeMethodTrace;
  const bool randomiseSelection   = true;

  int remainingTriesToFindSearchingTrace = 16;
  if (randomiseSelection) {
    _activeMethodTrace = peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling;

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
  }
  else {
    do {
      _activeMethodTrace = peano::datatraversal::autotuning::toMethodTrace( static_cast<int>(_activeMethodTrace) + 1 );
      if (_activeMethodTrace>=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling) {
        _activeMethodTrace = peano::datatraversal::autotuning::toMethodTrace( 0 );
      }
      remainingTriesToFindSearchingTrace--;
    }
    while (
      remainingTriesToFindSearchingTrace>0
      &&
      _measurements.count(_activeMethodTrace)==1
      &&
      _measurements[_activeMethodTrace]._searchDelta==0
    );
  }

  assertion(_measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0);
  assertion(_activeMethodTrace!=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling);

  if (
    oldActiveMethodTrace != peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling
    &&
    _measurements[_activeMethodTrace]._totalSerialTime > 0.0
    &&
    _measurements[oldActiveMethodTrace]._totalSerialTime > 0.0
    &&
    _measurements[_activeMethodTrace]._totalSerialTime < _measurements[oldActiveMethodTrace]._totalSerialTime
    &&
    _measurements[oldActiveMethodTrace]._searchDelta>0
  ) {
    logDebug( "changeMeasuredMethodTrace()", "roll back active method trace from " << toString(_activeMethodTrace) << " to " << toString(oldActiveMethodTrace) );
    _activeMethodTrace = oldActiveMethodTrace;
  }

  assertion(_measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0);
  if (
    _measurements[_activeMethodTrace]._currentMeasurement.getAccuracy() < _measurements[_activeMethodTrace]._currentMeasurement.getValue() * _measurements[_activeMethodTrace]._biggestProblemSize * _InitialRelativeAccuracy
  ) {
    _measurements[_activeMethodTrace]._currentMeasurement.increaseAccuracy(0.9);
  }


  assertion(_activeMethodTrace!=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling);
  assertion(_measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0);

  logDebug( "changeMeasuredMethodTrace()", "next active method trace " << toString(_activeMethodTrace) << " with existing measurement " << _measurements[_activeMethodTrace].toString() );
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::makeAttributesLearn() {
  assertion( _measurements[_activeMethodTrace]._currentMeasurement.isAccurateValue() );
  assertion( _measurements[_activeMethodTrace]._searchDelta>0 );
  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );

  if (_measurements[_activeMethodTrace]._currentGrainSize>=_measurements[_activeMethodTrace]._biggestProblemSize) {
    _measurements[_activeMethodTrace]._serialMeasurement  =  _measurements[_activeMethodTrace]._currentMeasurement.getValue();
  }

  if (
    _measurements[_activeMethodTrace]._currentMeasurement.getValue() < _measurements[_activeMethodTrace]._previousMeasuredTime
  ) {
    logInfo(
      "makeAttributesLearn()",
      "found better scaling parameter choice/serial runtime for " << toString(_activeMethodTrace) << ": " << _measurements[_activeMethodTrace].toString()
    );

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

    if ( _measurements[_activeMethodTrace]._biggestProblemSize<=tarch::multicore::Core::getInstance().getNumberOfThreads() ) {
      _measurements[_activeMethodTrace]._searchDelta--;
    }
    else if (
      _measurements[_activeMethodTrace]._searchDelta > tarch::multicore::Core::getInstance().getNumberOfThreads()
    ) {
      _measurements[_activeMethodTrace]._searchDelta /= tarch::multicore::Core::getInstance().getNumberOfThreads();
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

  assertion( _measurements[_activeMethodTrace]._currentGrainSize>0 );
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

  if (_measurements[_activeMethodTrace]._currentGrainSize==_measurements[_activeMethodTrace]._biggestProblemSize) {
    const double computeTime   = costPerProblemElement * static_cast<double>(problemSize);
    _measurements[_activeMethodTrace]._totalSerialTime += computeTime;
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
     std::string methodTrace;

/*     try {*/
     {
      leftToken   = rightString.substr( 0, rightString.find("=") );
      rightString = rightString.substr( leftToken.size()+1 );
      methodTrace = leftToken;
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
      logDebug( "loadStatistics(...)", "search delta is " <<  searchDelta );

      leftToken   = rightString.substr( 0, rightString.find(",") );
      rightString = rightString.substr( leftToken.size()+1 );
      double accuracy = std::stof(leftToken);
      logDebug( "loadStatistics(...)", "accuracy is " <<  accuracy << ". Ignore remainder of this line");

      leftToken   = rightString.substr( 0, rightString.find(",") );
      rightString = rightString.substr( leftToken.size()+1 );
      double totalSerialTime = std::stof(leftToken);
      logDebug( "loadStatistics(...)", "total serial time is " <<  totalSerialTime );

      leftToken   = rightString.substr( 0, rightString.find(",") );
      rightString = rightString.substr( leftToken.size()+1 );
      double serialMeasurement = std::stof(leftToken);
      logDebug( "loadStatistics(...)", "serial measurement is " <<  serialMeasurement << ". Ignore remainder of this line");

      _measurements[methodTrace]._biggestProblemSize    = biggestProblemSize;
      _measurements[methodTrace]._currentGrainSize      = currentGrainSize;
      _measurements[methodTrace]._previousMeasuredTime  = previousMeasuredTime;
      _measurements[methodTrace]._searchDelta           = searchDelta;
      _measurements[methodTrace]._currentMeasurement.erase();
      _measurements[methodTrace]._currentMeasurement.setAccuracy( accuracy );
      _measurements[methodTrace]._totalSerialTime       = totalSerialTime;
      _measurements[methodTrace]._serialMeasurement     = serialMeasurement;

      logDebug( "loadStatistics(...)", "added " << toString(methodTrace) << ": "<< _measurements[methodTrace].toString() );
     }
    /*
     catch (std::out_of_range& exception) {
      logWarning( "loadStatistics(...)", "failed to parse shared memory configuration file " << filename << " with error " << exception.what() << " for " << methodTrace << " in adapter " << oracleNumber);
     }
*/
    }

    // Older GCC versions require an explicit cast here
    tagOpen |= str.compare( "adapter-number=" + std::to_string( (long long)oracleNumber) )==0;
  }

  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::plotStatistics(std::ostream& out, int oracleNumber) const {
  out << "# " << std::endl;
  out << "# trace=biggest problem size, current grain size, previous measured time, search delta, accuracy, total serial time, serial measurement, current measurement" << std::endl;
  out << "# ----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  out << "# " << std::endl;
  out << "begin OracleForOnePhaseWithShrinkingGrainSize" << std::endl;
  out << "total-runtime=" << _totalRuntimeMeasurement.getValue() << std::endl;
  out << "adapter-number=" << oracleNumber << std::endl;

  for (auto measurement: _measurements) {
    out << peano::datatraversal::autotuning::toString(measurement.first)
        << "=" << measurement.second._biggestProblemSize
        << "," << measurement.second._currentGrainSize
        << "," << measurement.second._previousMeasuredTime
        << "," << measurement.second._searchDelta
        << "," << measurement.second._currentMeasurement.getAccuracy()
        << "," << measurement.second._totalSerialTime
        << "," << measurement.second._serialMeasurement
        << "," << measurement.second._currentMeasurement.toString();
    

    if (measurement.second._currentMeasurement.getValue()>1e-12) {
      out << ",estimated-speedup=" << measurement.second._serialMeasurement / measurement.second._currentMeasurement.getValue();
    }
    else if (measurement.second._previousMeasuredTime<_TimingMax) {
      out << ",estimated-speedup=" << measurement.second._serialMeasurement / measurement.second._previousMeasuredTime;
    }


    out << std::endl;
  }

  out << "end OracleForOnePhaseWithShrinkingGrainSize" << std::endl;
}


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::~OracleForOnePhaseWithShrinkingGrainSize() {
}


peano::datatraversal::autotuning::OracleForOnePhase* sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::createNewOracle() const {
  return new OracleForOnePhaseWithShrinkingGrainSize(_learn);
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::deactivateOracle() {
  if (_learn && holdsSearchingOracle()) {
    _totalRuntimeWatch.stopTimer();
    _totalRuntimeMeasurement.setValue( _totalRuntimeWatch.getCalendarTime() );
  }
}


bool sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::holdsSearchingOracle() const {
  bool result = false;

  for (auto measurement: _measurements) {
    result |= measurement.second._searchDelta>0;
  }

  return result;
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

    if (holdsSearchingOracle()) {
      _totalRuntimeWatch.startTimer();
    }
  }
}
