#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"
#include "peano/utils/Globals.h"
#include "tarch/Assertions.h"

#include <cstdlib>
#include <limits>
#include <fstream>


tarch::logging::Log  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_log( "sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize" );


#if defined(Asserts)
const double                       sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_InitialAccuracy(1e-2);
#else
const double                       sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_InitialAccuracy(1e-4);
#endif
const double                       sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_FinalAccuracy(1e-8);


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::OracleForOnePhaseWithShrinkingGrainSize():
  _activeMethodTrace(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling),
  _measurements() {
}


peano::datatraversal::autotuning::GrainSize  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelise(int problemSize, peano::datatraversal::autotuning::MethodTrace askingMethod) {
  if (
    _measurements.count(askingMethod)==0
    ||
    _measurements[askingMethod]._biggestProblemSize<problemSize
  ) {
    _measurements[askingMethod]._biggestProblemSize   = problemSize;
    _measurements[askingMethod]._currentGrainSize     = 0;
    _measurements[askingMethod]._currentMeasurement   = tarch::timing::Measurement(_InitialAccuracy);
    _measurements[askingMethod]._previousGrainSize    = 0;
    _measurements[askingMethod]._previousMeasuredTime = 0;
  }

  if (_measurements[askingMethod]._currentGrainSize >= problemSize) {
    return peano::datatraversal::autotuning::GrainSize(
      0,
      false,
      problemSize,
      askingMethod, nullptr
    );
  }
  else
    return peano::datatraversal::autotuning::GrainSize(
      _measurements[askingMethod]._currentGrainSize,
      _activeMethodTrace==askingMethod && !_measurements[askingMethod]._currentMeasurement.isAccurateValue(),
      problemSize,
      askingMethod, this
    );
}



void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::changeMeasuredMethodTrace() {
  _activeMethodTrace = peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling;

  while ( _measurements.count(_activeMethodTrace)==0 ) {
    _activeMethodTrace = peano::datatraversal::autotuning::toMethodTrace( rand() % (int)(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling) );
  }
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::makeAttributesLearn() {
  if (
    _measurements[_activeMethodTrace]._currentGrainSize==0
    &&
    _measurements[_activeMethodTrace]._currentMeasurement.getAccuracy()>_FinalAccuracy
  ) {
    _measurements[_activeMethodTrace]._previousGrainSize    = 0;
    _measurements[_activeMethodTrace]._previousMeasuredTime = _measurements[_activeMethodTrace]._currentMeasurement.getValue();
    _measurements[_activeMethodTrace]._currentGrainSize     = _measurements[_activeMethodTrace]._biggestProblemSize/2;
    _measurements[_activeMethodTrace]._currentMeasurement.erase();

    logInfo(
      "makeAttributesLearn()",
      "finished serial analysis, continue with"
      <<  " previous grain size=" << _measurements[_activeMethodTrace]._previousGrainSize
      << ", previous measured time=" << _measurements[_activeMethodTrace]._previousMeasuredTime
      << ", current grain size=" << _measurements[_activeMethodTrace]._currentGrainSize
      << ", current measurement=" << _measurements[_activeMethodTrace]._currentMeasurement.toString()
      << ", method-trace=" << toString(_activeMethodTrace)
    );
  }
  else if (
    _measurements[_activeMethodTrace]._currentMeasurement.getValue() < _measurements[_activeMethodTrace]._previousMeasuredTime
  ) {
    _measurements[_activeMethodTrace]._previousGrainSize     = _measurements[_activeMethodTrace]._currentGrainSize;
    _measurements[_activeMethodTrace]._previousMeasuredTime  = _measurements[_activeMethodTrace]._currentMeasurement.getValue();
    _measurements[_activeMethodTrace]._currentGrainSize     /= 2;
    _measurements[_activeMethodTrace]._currentMeasurement.erase();

    logInfo(
      "makeAttributesLearn()",
      "found better scaling parameter choice, continue with"
      <<  " previous grain size=" << _measurements[_activeMethodTrace]._previousGrainSize
      << ", previous measured time=" << _measurements[_activeMethodTrace]._previousMeasuredTime
      << ", current grain size=" << _measurements[_activeMethodTrace]._currentGrainSize
      << ", current measurement=" << _measurements[_activeMethodTrace]._currentMeasurement.toString()
      << ", method-trace=" << toString(_activeMethodTrace)
    );
  }
  else if (
    _measurements[_activeMethodTrace]._currentMeasurement.getValue() > _measurements[_activeMethodTrace]._previousMeasuredTime
    &&
    _measurements[_activeMethodTrace]._previousGrainSize==0
  ) {
    _measurements[_activeMethodTrace]._previousGrainSize     = 0;
    _measurements[_activeMethodTrace]._previousMeasuredTime  = 0;
    _measurements[_activeMethodTrace]._currentGrainSize      = 0;
    _measurements[_activeMethodTrace]._currentMeasurement.erase();
    _measurements[_activeMethodTrace]._currentMeasurement.increaseAccuracy(2.0);

    logInfo(
      "makeAttributesLearn()",
      "studied parameter configuration seems not to scale at all so switch back to serial setup, continue with"
      <<  " previous grain size=" << _measurements[_activeMethodTrace]._previousGrainSize
      << ", previous measured time=" << _measurements[_activeMethodTrace]._previousMeasuredTime
      << ", current grain size=" << _measurements[_activeMethodTrace]._currentGrainSize
      << ", current measurement=" << _measurements[_activeMethodTrace]._currentMeasurement.toString()
      << ", method-trace=" << toString(_activeMethodTrace)
    );
  }
  else if (
    _measurements[_activeMethodTrace]._currentMeasurement.getValue() > _measurements[_activeMethodTrace]._previousMeasuredTime
  ) {
    // _measurements[_activeMethodTrace]._previousGrainSize     = _measurements[_activeMethodTrace]._previousGrainSize;
    // _measurements[_activeMethodTrace]._previousMeasuredTime  = _measurements[_activeMethodTrace]._previousMeasuredTime;
    _measurements[_activeMethodTrace]._currentGrainSize      = _measurements[_activeMethodTrace]._previousGrainSize;
    _measurements[_activeMethodTrace]._currentMeasurement.erase();
    _measurements[_activeMethodTrace]._currentMeasurement.increaseAccuracy(2.0);

    logInfo(
      "makeAttributesLearn()",
      "previously studied parameter configuration seems to scale better, continue with previous parameter choice"
      <<  " previous grain size=" << _measurements[_activeMethodTrace]._previousGrainSize
      << ", previous measured time=" << _measurements[_activeMethodTrace]._previousMeasuredTime
      << ", current grain size=" << _measurements[_activeMethodTrace]._currentGrainSize
      << ", current measurement=" << _measurements[_activeMethodTrace]._currentMeasurement.toString()
      << ", method-trace=" << toString(_activeMethodTrace)
    );
  }
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelSectionHasTerminated(int problemSize, peano::datatraversal::autotuning::MethodTrace askingMethod, double costPerProblemElement) {
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

      leftToken   = rightString.substr( 0, rightString.find(",") );
      rightString = rightString.substr( leftToken.size()+1 );
      int   biggestProblemSize = std::stoi(leftToken);

      leftToken   = rightString.substr( 0, rightString.find(",") );
      rightString = rightString.substr( leftToken.size()+1 );
      int   currentGrainSize = std::stoi(leftToken);

      leftToken   = rightString.substr( 0, rightString.find(",") );
      rightString = rightString.substr( leftToken.size()+1 );
      int   previousGrainSize = std::stoi(leftToken);

      leftToken   = rightString.substr( 0, rightString.find(",") );
      rightString = rightString.substr( leftToken.size()+1 );
      double previousMeasuredTime = std::stof(leftToken);

      leftToken   = rightString.substr( 0, rightString.find("eps=") );
      rightString = rightString.substr( leftToken.size()+4 );
      leftToken   = rightString.substr( 0, rightString.find(",") );
      rightString = rightString.substr( leftToken.size()+1 );
      double accuracy = std::stof(leftToken);

      _measurements[methodTrace]._biggestProblemSize    = biggestProblemSize;
      _measurements[methodTrace]._previousGrainSize     = previousGrainSize;
      _measurements[methodTrace]._currentGrainSize      = currentGrainSize;
      _measurements[methodTrace]._previousMeasuredTime  = previousMeasuredTime;
      _measurements[methodTrace]._currentMeasurement.erase();
      _measurements[methodTrace]._currentMeasurement.setAccuracy( accuracy );
    }

    tagOpen |= str.compare( "adapter-number=" + std::to_string( oracleNumber) )==0;
  }
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::plotStatistics(std::ostream& out, int oracleNumber) const {
  out << "begin OracleForOnePhaseWithShrinkingGrainSize" << std::endl;
  out << "adapter-number=" << oracleNumber << std::endl;

  for (auto measurement: _measurements) {
    out << peano::datatraversal::autotuning::toString(measurement.first)
        << "=" << measurement.second._biggestProblemSize
        << "," << measurement.second._currentGrainSize
        << "," << measurement.second._previousGrainSize
        << "," << measurement.second._previousMeasuredTime
        << "," << measurement.second._currentMeasurement.toString()
        << std::endl;
  }

  out << "end OracleForOnePhaseWithShrinkingGrainSize" << std::endl;
}


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::~OracleForOnePhaseWithShrinkingGrainSize() {
}


peano::datatraversal::autotuning::OracleForOnePhase* sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::createNewOracle() const {
  return new OracleForOnePhaseWithShrinkingGrainSize();
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::deactivateOracle() {
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::activateOracle() {
  if (_measurements.empty()) {
    logDebug( "activateOracle()", "do not trigger any operation as measurement set is empty" );
  }
  else if (_activeMethodTrace==peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling) {
    changeMeasuredMethodTrace();
    logInfo(
      "activateOracle()", "no valid search had been triggered before. start search with "
      << " method-trace=" << toString(_activeMethodTrace)
    );
  }
  else if (_measurements[_activeMethodTrace]._currentMeasurement.isAccurateValue()) {
    makeAttributesLearn();
    changeMeasuredMethodTrace();
  }
}
