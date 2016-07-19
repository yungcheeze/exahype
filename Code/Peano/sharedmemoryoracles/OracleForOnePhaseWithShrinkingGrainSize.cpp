#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"
#include "peano/utils/Globals.h"
#include "tarch/Assertions.h"

#include <cstdlib>
#include <limits>
#include <fstream>


tarch::logging::Log  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_log( "sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize" );


int           sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_activeMethodTrace(0);
const double  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_InitialAccuracy(1e-4);


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::OracleForOnePhaseWithShrinkingGrainSize(
  int                                                   adapterNumber,
  const peano::datatraversal::autotuning::MethodTrace&  methodTrace
):
  _methodTrace(methodTrace),
  _adapterNumber(adapterNumber),
  _currentSearchDelta(1),
  _biggestProblemSize(0),
  _currentGrainSize(std::numeric_limits<int>::max()),
  _currentMeasurement(_InitialAccuracy),
  _previousMeasuredTime(-1.0),
  _lastProblemSize(0.0) {
}


std::pair<int,bool> sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelise(int problemSize) {
  assertion( problemSize>0 );

  _biggestProblemSize = std::max( _biggestProblemSize,problemSize );
  _lastProblemSize    = problemSize;

  if (problemSize <= _currentGrainSize) {
    return std::pair<int,bool>(0,_currentSearchDelta>0);
  }
  else {
    return std::pair<int,bool>(_currentGrainSize,_currentSearchDelta>0);
  }
}



void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::changeMeasuredMethodTrace() {
  _activeMethodTrace = rand() % peano::datatraversal::autotuning::NumberOfDifferentMethodsCalling;

  assertion1(_activeMethodTrace>=0,_activeMethodTrace);
  assertion1(_activeMethodTrace<peano::datatraversal::autotuning::NumberOfDifferentMethodsCalling,_activeMethodTrace);
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::makeAttributesLearn() {
  assertion( _currentMeasurement.isAccurateValue() );

  if (_biggestProblemSize < _currentGrainSize && _currentSearchDelta==0) {
    _currentSearchDelta   = 1;
    _biggestProblemSize   = 0;
    _currentGrainSize     = std::numeric_limits<int>::max();
    _previousMeasuredTime = -1.0;
    _lastProblemSize      = 0.0;

    _currentMeasurement   = tarch::timing::Measurement(_InitialAccuracy);

    logInfo(
      "makeAttributesLearn()",
      "found new problem size that has been unknown before. Restart analysis. grain size=" << _currentGrainSize
      << ", method-trace=" << toString(_methodTrace) << ", " << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << "th adapter"
    );
  }
  else if (_biggestProblemSize < _currentGrainSize) {
    assertion(_previousMeasuredTime==-1.0);

    _currentSearchDelta   = std::max(_biggestProblemSize/2,1);
    _currentGrainSize     = _currentSearchDelta;
    _previousMeasuredTime = _currentMeasurement.getValue();

    _currentMeasurement.erase();

    // This is not an info as it is a standard behaviour that has to arise at least once.
    logDebug(
      "makeAttributesLearn()",
      "identified new maximum problem size and/or finished serial analysis. new grain size=" << _currentGrainSize << ", biggest-problem-size=" << _biggestProblemSize
      << ", currentSearchDelta=" << _currentSearchDelta
      << ", method-trace=" << toString(_methodTrace) << ", " << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << "th adapter"
    );
  }
  else if (
    _previousMeasuredTime > _currentMeasurement.getValue() 
    &&
    _currentSearchDelta>0
  ) {
    logInfo(
      "makeAttributesLearn()",
      "found scaling parameter configuration with grain size=" << _currentGrainSize
      << ", currentSearchDelta=" << _currentSearchDelta
      << ", biggest-problem-size=" << _biggestProblemSize
      << ", method-trace=" << toString(_methodTrace) << ", " << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << "th adapter"
    );

    if (_currentGrainSize > _currentSearchDelta) {
      _currentGrainSize     -= _currentSearchDelta;
      _previousMeasuredTime  = _currentMeasurement.getValue();
      _currentMeasurement.erase();

      logInfo(
        "makeAttributesLearn()",
        "continue to study smaller grain size=" << _currentGrainSize
        << ", currentSearchDelta=" << _currentSearchDelta
        << ", biggest-problem-size=" << _biggestProblemSize
        << ", method-trace=" << toString(_methodTrace) << ", " << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << "th adapter"
      );
    } 
    else {
      _currentSearchDelta /= 2;
      _currentMeasurement.increaseAccuracy(2.0);  
      
      logInfo(
        "makeAttributesLearn()",
        "reduce delta (zero implies switch off), increase sensitivity and continue to study grain size=" << _currentGrainSize
        << ", currentSearchDelta=" << _currentSearchDelta
        << ", biggest-problem-size=" << _biggestProblemSize
        << ", method-trace=" << toString(_methodTrace) << ", " << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << "th adapter"
      );
    }
  }
  else if (
    _previousMeasuredTime < _currentMeasurement.getValue()
    &&
    _currentSearchDelta>0
  ) {
    logInfo(
      "makeAttributesLearn()",
      "studied parameter configuration does not scale. grain size=" << _currentGrainSize
      << ", currentSearchDelta=" << _currentSearchDelta
      << ", biggest-problem-size=" << _biggestProblemSize
      << ", method-trace=" << toString(_methodTrace) << ", " << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << "th adapter"
    );

    _currentGrainSize     += _currentSearchDelta;
    if (_currentGrainSize>_biggestProblemSize/2) {
      _currentGrainSize   = 0;
      _currentSearchDelta = 0;
      logInfo(
        "makeAttributesLearn()",
        "switch off multithreading for " 
        << "method-trace=" << toString(_methodTrace) << ", " << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << "th adapter"
      );
    }
    else {
      _currentSearchDelta   /=2;
      _previousMeasuredTime  = std::numeric_limits<double>::max();
 
      _currentMeasurement.erase();
      _currentMeasurement.increaseAccuracy(2.0);
 
      logInfo(
        "makeAttributesLearn()",
        "fall back to previous grain size, reduce search step size and increase sensitivity. grain size=" << _currentGrainSize
        << ", biggest-problem-size=" << _biggestProblemSize
        << ", currentSearchDelta=" << _currentSearchDelta
//        << ", accuracy=" << _currentMeasurement.toString()
        << ", method-trace=" << toString(_methodTrace) << ", " << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << "th adapter"
      );
    }
  }
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelSectionHasTerminated(double elapsedCalendarTime) {
  assertion(_currentSearchDelta>0);
  assertion(_lastProblemSize!=0.0);

  _currentMeasurement.setValue(elapsedCalendarTime/_lastProblemSize);
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::loadStatistics(const std::string& filename) {
  std::ifstream file(filename);
  std::string str;
// @todo Please remove as soon as Intel supports this properly (according to standard should work without cast)
  std::string rightEntry = "adapter=" + std::to_string( (long long)(_adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1) ) + ", method=" + toString(_methodTrace) + ":";
  while (std::getline(file, str)) {
    if (str.substr(0,rightEntry.size())==rightEntry) {
      str = str.substr(  str.find("currentSearchDelta=") + std::string("currentSearchDelta=").length() );
      std::string currentSearchDelta = str.substr(0,str.find(","));

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

      _currentSearchDelta   = atoi( currentSearchDelta.c_str() );
      _biggestProblemSize   = atoi( biggestProblemSize.c_str() );
      _currentGrainSize     = atoi( currentGrainSize.c_str() );
      _previousMeasuredTime = atof( previousMeasuredTime.c_str() );
      _lastProblemSize      = atoi( lastProblemSize.c_str() );
      _currentMeasurement   = tarch::timing::Measurement(atof(accuracy.c_str()));

      logInfo( "loadStatistics(std::string)", "found data for " << rightEntry << " " << _currentSearchDelta << "," << _biggestProblemSize << "," << _currentGrainSize << "," << _previousMeasuredTime << ":" << _lastProblemSize << "," << _currentMeasurement.toString());
    }
  }
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::plotStatistics(std::ostream& out) const {
  if (_currentMeasurement.getNumberOfMeasurements()>0) {
    out << "adapter=" << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters+1 << ", method=" << toString(_methodTrace) << ":";
    out << "currentSearchDelta="    << _currentSearchDelta
        << ",biggestProblemSize="   << _biggestProblemSize
        << ",currentGrainSize="     << _currentGrainSize
        << ",previousMeasuredTime=" << _previousMeasuredTime
        << ",lastProblemSize="      << _lastProblemSize
   //     << ",accuracy="             << _currentMeasurement._accuracy()
        << ",currentMeasurement="   << _currentMeasurement.toString();

    if (_biggestProblemSize < _currentGrainSize && _currentSearchDelta>0) {
      out << " (still determining serial runtime)" << std::endl;
    }
    else if (_currentSearchDelta>0) {
      out << " (still searching for optimal grain size)" << std::endl;
    }
    else if (_currentSearchDelta==0 && _currentGrainSize==0) {
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
  if (
    _currentMeasurement.isAccurateValue()
    &&
    peano::datatraversal::autotuning::toMethodTrace( _activeMethodTrace ) == _methodTrace
  ) {
    makeAttributesLearn();
  }

  if (
    peano::datatraversal::autotuning::toMethodTrace( _activeMethodTrace ) == _methodTrace
  ) {
    changeMeasuredMethodTrace();

    logDebug( "informAboutElapsedTimeOfLastTraversal(double)", "new active method trace is=" << peano::datatraversal::autotuning::toMethodTrace( _activeMethodTrace ) );
  }
}
