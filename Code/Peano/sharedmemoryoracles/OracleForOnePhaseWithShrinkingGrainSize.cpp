#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"
#include "peano/utils/Globals.h"
#include "tarch/Assertions.h"

#include <cstdlib>
#include <limits>
#include <fstream>


tarch::logging::Log  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_log( "sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize" );


int  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_activeMethodTrace(0);
int  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_delayBetweenTwoUpdates(0);


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::OracleForOnePhaseWithShrinkingGrainSize(
  int                                                   adapterNumber,
  const peano::datatraversal::autotuning::MethodTrace&  methodTrace
):
  _methodTrace(methodTrace),
  _adapterNumber(adapterNumber),
  _oracleIsSearching(true),
  _biggestProblemSize(0),
  _currentGrainSize(std::numeric_limits<int>::max()),
  _currentMeasurement(1.0e-0), // magic constant war vorher 1.0e-2
  _previousMeasuredTime(-1.0),
  _lastProblemSize(0.0) {
}


std::pair<int,bool> sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelise(int problemSize) {
  assertion( problemSize>0 );

  if (_biggestProblemSize<problemSize) _biggestProblemSize = problemSize;

  _lastProblemSize = problemSize;

  if (problemSize <= _currentGrainSize) {
    return std::pair<int,bool>(0,_oracleIsSearching);
  }
  else {
    return std::pair<int,bool>(_currentGrainSize,_oracleIsSearching);
  }
}



void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::changeMeasuredMethodTrace() {
  _activeMethodTrace++;
  if (_activeMethodTrace>=peano::datatraversal::autotuning::NumberOfDifferentMethodsCalling) {
    _activeMethodTrace = 0;
  }
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::makeAttributesLearn() {
  assertion( _currentMeasurement.isAccurateValue() );

  // If max problem size has increased, notably as soon as serial size is known
  if (_biggestProblemSize < _currentGrainSize) {
    assertion(_previousMeasuredTime==-1.0);
    assertion(_oracleIsSearching);
    _currentGrainSize   = _biggestProblemSize / 2 + 1;
    _oracleIsSearching  = true;

    // This is not an info as it is a standard behaviour that has to arise at least once.
    logDebug(
      "makeAttributesLearn()",
      "identified new maximum problem size and/or finished serial analysis. new grain size=" << _currentGrainSize << ", biggest-problem-size=" << _biggestProblemSize
      << ", method-trace=" << toString(_methodTrace) << ", " << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << "th adapter"
    );
  }
  else if (
    _previousMeasuredTime > _currentMeasurement.getValue() &&
    _currentGrainSize >= 2
  ) {
    _currentGrainSize /= 2;
    _oracleIsSearching  = true;

    logInfo(
      "makeAttributesLearn()",
      "algorithm part seems to scale. Continue to search for higher performance with halved grain size. new grain size=" << _currentGrainSize
      << ", oracle continues to search=" << _oracleIsSearching
     << ", biggest-problem-size=" << _biggestProblemSize
      << ", method-trace=" << toString(_methodTrace) << ", " << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << "th adapter"
    );
  }
  else if (
    _previousMeasuredTime > _currentMeasurement.getValue()
  ) {
    _oracleIsSearching  = false;

    logInfo(
      "makeAttributesLearn()",
      "algorithm does benefit from smallest grain size possible or is inherently sequential. Stop search. grain size=" << _currentGrainSize
     << ", biggest-problem-size=" << _biggestProblemSize
      << ", method-trace=" << toString(_methodTrace) << ", " << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << "th adapter"
    );
  }
  else {
    _oracleIsSearching  = false;
    _currentGrainSize  *= 2;

    logInfo(
      "makeAttributesLearn()",
      "algorithm does not benefit from small grain size currently studied. Stop search. grain size=" << _currentGrainSize
     << ", biggest-problem-size=" << _biggestProblemSize
      << ", method-trace=" << toString(_methodTrace) << ", " << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << "th adapter"
    );
  }

  _previousMeasuredTime = _currentMeasurement.getValue();
  _currentMeasurement.erase();
  _currentMeasurement.increaseAccuracy(2.0);
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelSectionHasTerminated(double elapsedCalendarTime) {
  assertion(_oracleIsSearching);
  assertion(_lastProblemSize!=0.0);

  _currentMeasurement.setValue(elapsedCalendarTime/_lastProblemSize);
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::loadStatistics(const std::string& filename) {
  std::ifstream file(filename);
  std::string str;
  std::string rightEntry = "adapter=" + std::to_string( _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 ) + ", method=" + toString(_methodTrace) + ":";
  while (std::getline(file, str)) {
    if (str.substr(0,rightEntry.size())==rightEntry) {
      str = str.substr(  str.find("oracleIsSearching=") + std::string("oracleIsSearching=").length() );
      std::string oracleIsSearching = str.substr(0,str.find(","));

      str = str.substr(  str.find("biggestProblemSize=") + std::string("biggestProblemSize=").length() );
      std::string biggestProblemSize = str.substr(0,str.find(","));

      str = str.substr(  str.find("currentGrainSize=") + std::string("currentGrainSize=").length() );
      std::string currentGrainSize = str.substr(0,str.find(","));

      str = str.substr(  str.find("previousMeasuredTime=") + std::string("previousMeasuredTime=").length() );
      std::string previousMeasuredTime = str.substr(0,str.find(","));

      str = str.substr(  str.find("lastProblemSize=") + std::string("lastProblemSize=").length() );
      std::string lastProblemSize = str.substr(0,str.find(","));

      str = str.substr(  str.find("accuracy=") + std::string("accuracy=").length() );
      std::string accuracy = str.substr(0,str.find(","));

      _oracleIsSearching    = oracleIsSearching=="1";
      _biggestProblemSize   = atoi( biggestProblemSize.c_str() );
      _currentGrainSize     = atoi( currentGrainSize.c_str() );
      _previousMeasuredTime = atof( previousMeasuredTime.c_str() );
      _lastProblemSize      = atoi( lastProblemSize.c_str() );
      _currentMeasurement   = tarch::timing::Measurement(atof(accuracy.c_str()));

      logInfo( "loadStatistics(std::string)", "found data for " << rightEntry << " " << _oracleIsSearching << "," << _biggestProblemSize << "," << _currentGrainSize << "," << _previousMeasuredTime << ":" << _lastProblemSize << "," << _currentMeasurement.toString());
    }
  }
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::plotStatistics(std::ostream& out) const {
  if (_currentMeasurement.getNumberOfMeasurements()>0 || !_oracleIsSearching) {
    out << "adapter=" << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << ", method=" << toString(_methodTrace) << ":";
    out << "oracleIsSearching="     << _oracleIsSearching
        << ",biggestProblemSize="   << _biggestProblemSize
        << ",currentGrainSize="     << _currentGrainSize
        << ",previousMeasuredTime=" << _previousMeasuredTime
        << ",lastProblemSize="      << _lastProblemSize
        << ",accuracy="             << _currentMeasurement.getAccuracy()
        << ",currentMeasurement="   << _currentMeasurement.toString();

    if (_biggestProblemSize < _currentGrainSize && _oracleIsSearching) {
      out << " (still determining serial runtime)" << std::endl;
    }
    else if (_oracleIsSearching) {
      out << " (still searching for optimal grain size)" << std::endl;
    }
    else if (_biggestProblemSize<=_currentGrainSize) {
      out << "  (does not scale, oracle is not searching anymore)" << std::endl;
    }
    else {
      out << " (scales, oracle is not searching anymore)" << std::endl;
    }
  }
}


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::~OracleForOnePhaseWithShrinkingGrainSize() {
}


peano::datatraversal::autotuning::OracleForOnePhase* sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::createNewOracle(int adapterNumber, const peano::datatraversal::autotuning::MethodTrace& methodTrace) const {
  return new OracleForOnePhaseWithShrinkingGrainSize(adapterNumber, methodTrace);
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::informAboutElapsedTimeOfLastTraversal(double elapsedTime) {
  _delayBetweenTwoUpdates++;

  if (
    _currentMeasurement.isAccurateValue()
    &&
    peano::datatraversal::autotuning::toMethodTrace( _activeMethodTrace ) == _methodTrace
    &&
    _delayBetweenTwoUpdates > peano::datatraversal::autotuning::NumberOfDifferentMethodsCalling * 32
  ) {
    makeAttributesLearn();
    changeMeasuredMethodTrace();
    _delayBetweenTwoUpdates = 0;
  }
  else if (
    _currentMeasurement.getNumberOfMeasurements()==0
  ) {
    changeMeasuredMethodTrace();
  }
}
