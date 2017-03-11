#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"
#include "peano/utils/Globals.h"
#include "tarch/Assertions.h"
#include "tarch/multicore/Core.h"


#include <cstdlib>
#include <limits>
#include <fstream>
#include <stdexcept>


tarch::logging::Log  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_log( "sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize" );

//const double   sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_InitialRelativeAccuracy(1e-2);
const double   sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_InitialRelativeAccuracy(1e2);
const double   sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_MaxAccuracy( 1.0 );
const double   sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_WideningFactor( 0.9 );


bool  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_hasLearnedSinceLastQuery( false );


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::OracleForOnePhaseWithShrinkingGrainSize(
  bool learn,
  bool restart,
  SelectNextStudiedMeasureTraceStrategy selectNextStudiedMeasureTraceStrategy
):
  _activeMethodTrace(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling),
  _learn(learn && tarch::multicore::Core::getInstance().getNumberOfThreads()>1),
  _restart(restart),
  _selectNextStudiedMeasureTraceStrategy(selectNextStudiedMeasureTraceStrategy),
  _measurements() {
}


bool sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::operator<(const DatabaseEntry& cmp) const {
  return ( _biggestProblemSize<cmp._biggestProblemSize );
}


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry&  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::getDatabaseEntry(
  int problemSize,
  peano::datatraversal::autotuning::MethodTrace askingMethod
) {
  if ( _measurements.count(askingMethod)==0 ) {
    _measurements.insert( std::pair<peano::datatraversal::autotuning::MethodTrace,MethodTraceData >(askingMethod,MethodTraceData()) );
    _measurements[askingMethod].push_back( DatabaseEntry(2) );
    assertion( _measurements.count(askingMethod)==1 );
    logInfo(
      "getDatabaseEntry(int)",
      "inserted trivial entry for " + peano::datatraversal::autotuning::toString(askingMethod)
      << ": " << _measurements[askingMethod].rbegin()->toString()
    );
  }

  while (_measurements[askingMethod].rbegin()->getBiggestProblemSize()<problemSize) {
    _measurements[askingMethod].push_back( DatabaseEntry(
      *_measurements[askingMethod].rbegin(),
      _measurements[askingMethod].rbegin()->getBiggestProblemSize()*2
    ));
    logInfo(
      "getDatabaseEntry(int)",
      "inserted new entry for " + peano::datatraversal::autotuning::toString(askingMethod)
      << ": " << _measurements[askingMethod].rbegin()->toString()
    );
  }

  MethodTraceData::iterator entry = _measurements[askingMethod].begin();
  while (entry->getBiggestProblemSize()<problemSize) {
    entry++;
  }

  DatabaseEntry& result = *entry;
  return result;
}


peano::datatraversal::autotuning::GrainSize  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelise(int problemSize, peano::datatraversal::autotuning::MethodTrace askingMethod) {
  assertion( askingMethod != peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling );

  assertion(askingMethod!=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling);
  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );

  auto databaseEntry = getDatabaseEntry(problemSize,askingMethod);

  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );
  assertion( _measurements.count(askingMethod)>0 );

  const bool trackTime         = (_activeMethodTrace==askingMethod) &&  (databaseEntry.isSearching());

  logDebug(
    "parallelise()",
    "will track time for " << toString(askingMethod) << "=" << trackTime <<
    " with active trace " << toString(_activeMethodTrace) <<
    " for maximum problem size " << problemSize
  );

  if (trackTime) {
    const int chosenParallelGrainSize = databaseEntry.isStudyingScalingSetup() ? databaseEntry.getCurrentGrainSize() : 0;
    return peano::datatraversal::autotuning::GrainSize(
      rand()%10==0 ? 0 : chosenParallelGrainSize,
      trackTime,
      problemSize,
      askingMethod, this
    );
  }
  else if (_activeMethodTrace==peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling) {
    return peano::datatraversal::autotuning::GrainSize(
      0,
      false,
      problemSize,
      askingMethod, this
    );
  }
  else {
    const int chosenParallelGrainSize = databaseEntry.isStudyingScalingSetup() ? databaseEntry.getCurrentGrainSize() : 0;
    return peano::datatraversal::autotuning::GrainSize(
      chosenParallelGrainSize,
      trackTime,
      problemSize,
      askingMethod, this
    );
  }
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::widenAccuracyOfCurrentlyStudiedMethodTraceAndRandomlyRestart() {
  bool widenedEntry = false;
  for (auto& p: _measurements[_activeMethodTrace]) {
    widenedEntry |= p.widenAccuracy();
  }
  if (!widenedEntry && _restart) {
    for (auto& p: _measurements[_activeMethodTrace]) {
      if (rand()%100==0) {
        p.restart();
      }
    }
  }
}


bool sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::isStableDatabaseEntrySet( const peano::datatraversal::autotuning::MethodTrace&  askingMethod) const {
  if ( _measurements.count(askingMethod)==0 ) {
    return true;
  }
  else {
    bool result             = true;

    for (auto p: _measurements.at(askingMethod)) {
      result             &= (!p.isSearching());
    }

    return result;
  }
}



void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::changeMeasuredMethodTrace() {
  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );
  assertion( !_measurements.empty() );

  // We just make it big enough to run through every variant if we used a deterministic scheme.
  int remainingTriesToFindSearchingTrace = (int)(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling);

  switch (_selectNextStudiedMeasureTraceStrategy) {
    case SelectNextStudiedMeasureTraceStrategy::Randomised:
      _activeMethodTrace             = peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling;

      while ( _measurements.count(_activeMethodTrace)==0 ) {
        _activeMethodTrace = peano::datatraversal::autotuning::toMethodTrace( rand() % (int)(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling) );
        if (
          remainingTriesToFindSearchingTrace>0
          &&
          _measurements.count(_activeMethodTrace)==1
          &&
          isStableDatabaseEntrySet(_activeMethodTrace)
        ) {
          remainingTriesToFindSearchingTrace--;
          logDebug( "changeMeasuredMethodTrace()", "skip " << toString(_activeMethodTrace) );
          _activeMethodTrace = peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling;
        }
      }
      break;

    case SelectNextStudiedMeasureTraceStrategy::Cyclic:
      do {
        _activeMethodTrace = peano::datatraversal::autotuning::toMethodTrace( static_cast<int>(_activeMethodTrace) + 1 );
        if (_activeMethodTrace>=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling) {
          _activeMethodTrace = peano::datatraversal::autotuning::toMethodTrace( 0 );
        }
        remainingTriesToFindSearchingTrace--;
      }
      while (
        _measurements.count(_activeMethodTrace)==0
        ||
        ( remainingTriesToFindSearchingTrace>0 && isStableDatabaseEntrySet(_activeMethodTrace))
      );
      break;

    case SelectNextStudiedMeasureTraceStrategy::RandomisedWithHigherPrioritiesForSmallProblemSizes:
      {
        _activeMethodTrace                               = peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling;
        static int maxProblemSizeOfLastMethodTraceChoice = 1;

        while ( _measurements.count(_activeMethodTrace)==0 ) {
          _activeMethodTrace = peano::datatraversal::autotuning::toMethodTrace( rand() % (int)(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling) );
          if (
            remainingTriesToFindSearchingTrace>0
            &&
            _measurements.count(_activeMethodTrace)==1
            &&
            (
              isStableDatabaseEntrySet(_activeMethodTrace)
              ||
              _measurements[_activeMethodTrace].back().getBiggestProblemSize() > maxProblemSizeOfLastMethodTraceChoice
            )
          ) {
            remainingTriesToFindSearchingTrace--;
            logDebug( "changeMeasuredMethodTrace()", "skip " << toString(_activeMethodTrace) << " with internal max size=" << maxProblemSizeOfLastMethodTraceChoice);
            _activeMethodTrace = peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling;
          }
        }
        maxProblemSizeOfLastMethodTraceChoice = (maxProblemSizeOfLastMethodTraceChoice + _measurements[_activeMethodTrace].back().getBiggestProblemSize()) / 2;
      }
      break;
  }

  assertion(_measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0);
  assertion(_activeMethodTrace!=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling);

  assertion(_activeMethodTrace!=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling);
  assertion(_measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0);

  logInfo(
    "changeMeasuredMethodTrace()", "next active method trace " << toString(_activeMethodTrace) << " after " <<
    ((int)(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)-remainingTriesToFindSearchingTrace) <<
    " search iterations to identify next analysed method trace"
  );
}



bool sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::isAccurate() const {
  return _currentParallelMeasurement.isAccurateValue() && _currentSerialMeasurement.isAccurateValue();
}


bool sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::isSearching() const {
  return _searchDelta>0;
}


bool sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::isScaling() const {
  return (_currentGrainSize<_biggestProblemSize && _searchDelta==0)
      || (_searchDelta>0 && _currentParallelMeasurement.getValue() < _currentSerialMeasurement.getValue() );
}


bool sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::isStudyingScalingSetup() const {
  return _currentGrainSize<_biggestProblemSize;
}


int sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::getCurrentGrainSize() const {
  return _currentGrainSize;
}


int  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::getBiggestProblemSize() const {
  return _biggestProblemSize;
}



void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::learn() {
  assertion( _currentSerialMeasurement.isAccurateValue() );
  assertion( _currentParallelMeasurement.isAccurateValue() );
  assertion( _searchDelta>0 );

  const double newSpeedup = _currentSerialMeasurement.getValue() / _currentParallelMeasurement.getValue();

  if ( newSpeedup > _previousSpeedup ) {
    logInfo( "learn()", "found better scaling parameter choice/serial runtime for: " << toString() );

    while ( _currentGrainSize - _searchDelta <= 0 && _searchDelta>0 ) {
      _searchDelta /= 2;
    }

    _currentGrainSize     -= _searchDelta;
    _currentParallelMeasurement.erase();
    _currentParallelMeasurement.increaseAccuracy(2.0);
    _currentSerialMeasurement.increaseAccuracy(2.0);
    _previousSpeedup       = newSpeedup;
  }
  else {
    logInfo( "learn()", "parameter choice for " << toString() << " does not scale" );

    _currentGrainSize     += _searchDelta;
    _previousSpeedup       = 0;
    _currentParallelMeasurement.erase();
    _currentParallelMeasurement.increaseAccuracy(2.0);
    _currentSerialMeasurement.increaseAccuracy(2.0);
    _previousSpeedup       = 1.0;

    if ( _biggestProblemSize<=tarch::multicore::Core::getInstance().getNumberOfThreads() ) {
      _searchDelta--;
    }
    else if (
      _searchDelta > tarch::multicore::Core::getInstance().getNumberOfThreads()
    ) {
      _searchDelta /= std::max(2,tarch::multicore::Core::getInstance().getNumberOfThreads());
    }
    else {
      _searchDelta /= 2;
    }
  }

  if (_currentGrainSize>_biggestProblemSize/2) {
    _searchDelta = 0;
    logInfo( "learn()", "stop search for " << toString() << " as it searches in non-scaling regime");
  }
  else {
    logInfo( "learn()", "continue with " << toString() );
  }

  assertion1(_currentGrainSize>0,  toString() );
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelSectionHasTerminated(int problemSize, int grainSize, peano::datatraversal::autotuning::MethodTrace askingMethod, double costPerProblemElement) {
  assertion( askingMethod!=peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling );

  auto& databaseEntry = getDatabaseEntry(problemSize,askingMethod);
  if (!databaseEntry.isAccuracyInitialised() ) {
    const double computeTime   = costPerProblemElement * static_cast<double>(problemSize);
    databaseEntry.initAccuracy(computeTime);
    logDebug(
      "parallelSectionHasTerminated(...)",
      "fix accuracy for " << toString(askingMethod) << " to " << databaseEntry.toString()
    );
  }

  databaseEntry.setValue( grainSize, costPerProblemElement );
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::setValue(int grainSize, double value) {
  if (grainSize==0) {
    _currentSerialMeasurement.setValue( value );
  }
  else {
    _currentParallelMeasurement.setValue( value );
  }
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::initAccuracy(double serialComputeTime) {
  _currentSerialMeasurement.setAccuracy(   serialComputeTime * _InitialRelativeAccuracy );
  _currentParallelMeasurement.setAccuracy( serialComputeTime * _InitialRelativeAccuracy );
}


bool sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::widenAccuracy() {
  bool widenThisEntry = false;

  if (_currentSerialMeasurement.getAccuracy()!=0.0 && _currentSerialMeasurement.getAccuracy()<_MaxAccuracy ) {
    _currentSerialMeasurement.increaseAccuracy( _WideningFactor );
    widenThisEntry = true;
  }

  if (_currentParallelMeasurement.getAccuracy()!=0.0 && _currentParallelMeasurement.getAccuracy()<_MaxAccuracy ) {
    _currentParallelMeasurement.increaseAccuracy( _WideningFactor );
    widenThisEntry = true;
  }

  return widenThisEntry;
}


std::string sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::toString() const {
  std::ostringstream msg;
  msg <<        _biggestProblemSize
      << "," << _currentGrainSize
      << "," << _searchDelta
      << "," << _currentSerialMeasurement.getAccuracy()
      << "," << _currentParallelMeasurement.getAccuracy()
      << "," << _previousSpeedup
      << "," << _currentSerialMeasurement.toString()
      << "," << _currentParallelMeasurement.toString()
      << ")";
  return msg.str();
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::plotStatistics(std::ostream& out, int oracleNumber) const {
  out << "# " << std::endl;
  out << "# trace=biggest problem size, current grain size, previous measured time, search delta, accuracy, current measurement" << std::endl;
  out << "# -------------------------------------------------------------------------------------------------------------------" << std::endl;
  out << "# " << std::endl;
  out << "begin OracleForOnePhaseWithShrinkingGrainSize" << std::endl;
  out << "initial-relative-accuracy=" << _InitialRelativeAccuracy << std::endl;
  out << "adapter-number=" << oracleNumber << std::endl;

  for (auto measurement: _measurements)
  for (auto p: measurement.second) {
    out << peano::datatraversal::autotuning::toString(measurement.first)
        << "=" << p.toString() << std::endl;
  }

  out << "end OracleForOnePhaseWithShrinkingGrainSize" << std::endl;

  logDebug( "plotStatistics(std::ostream,int)", "piped statistics for oracle no " << oracleNumber );
}


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::DatabaseEntry(int problemSize) {
  assertion(problemSize>0);

  if (
    problemSize < tarch::multicore::Core::getInstance().getNumberOfThreads()*2
    ||
    tarch::multicore::Core::getInstance().getNumberOfThreads() <= 2
  ) {
    _searchDelta = problemSize/2;
  }
  else {
    _searchDelta = problemSize - problemSize / tarch::multicore::Core::getInstance().getNumberOfThreads();
  }

  _biggestProblemSize         = problemSize;
  _currentGrainSize           = problemSize - _searchDelta;
  _currentSerialMeasurement   = tarch::timing::Measurement( 0.0 );
  _currentParallelMeasurement = tarch::timing::Measurement( 0.0 );
  _previousSpeedup            = 1.0;
}


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::DatabaseEntry( const DatabaseEntry& prototype, int newProblemSize ):
  DatabaseEntry(newProblemSize) {
  assertion(prototype._biggestProblemSize<newProblemSize);
  if (prototype._currentGrainSize < prototype._biggestProblemSize && prototype.isAccuracyInitialised() ) {
    // These are all integers and it should be possible to  determine all
    // settings through integer arithmetics if we do not determine a rescaling
    // factor (integer) once and then add it but first multiply and then
    // divide. However, if we do so, we quickly run into integer overflows:
    //     _currentGrainSize = prototype._currentGrainSize * newProblemSize / prototype._biggestProblemSize;
    //     _searchDelta      = prototype._searchDelta      * newProblemSize / prototype._biggestProblemSize;
    //
    // So we have to work with first computing the ratio though it might introduce truncation errors:
    int scaleUp = newProblemSize / prototype._biggestProblemSize;
    if (scaleUp==0) {
      scaleUp = 1;
    }

    _currentGrainSize = prototype._currentGrainSize * scaleUp;
    _searchDelta      = prototype._searchDelta      * scaleUp;
  }

  if (_searchDelta==0) {
    restart();
  }

  assertion1(_currentGrainSize<=_biggestProblemSize, toString() );
  assertion1(_currentGrainSize>0,                    toString() );
  assertion1(_biggestProblemSize>0,                  toString() );
  assertion1(_searchDelta>=0,                        toString() );
}


bool sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::isAccuracyInitialised() const {
  return _currentParallelMeasurement.getAccuracy()>0.0;
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::restart() {
  assertion1( _currentGrainSize>0, toString() );

  if (
    _biggestProblemSize>tarch::multicore::Core::getInstance().getNumberOfThreads()
    &&
    _currentGrainSize<_biggestProblemSize/4
  ) {
    _searchDelta        = _currentGrainSize;
    _currentGrainSize   = _currentGrainSize * 2;
  }
  else {
    _searchDelta        = _biggestProblemSize/2;
    _currentGrainSize   = _biggestProblemSize;
  }

  _currentSerialMeasurement   = tarch::timing::Measurement( 0.0 );
  _currentParallelMeasurement = tarch::timing::Measurement( 0.0 );
  _previousSpeedup            = 1.0;

  if (_searchDelta<=0) {
    _searchDelta = _currentGrainSize/2;
  }

  logInfo( "restart(...)", "restarted " << toString() );
}


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::DatabaseEntry::DatabaseEntry(std::string  inputLine) {
  std::string leftToken    = "";
  std::string rightString  = inputLine;
  std::string methodTrace;

  leftToken   = rightString.substr( 0, rightString.find(",") );
  rightString = rightString.substr( leftToken.size()+1 );
  _biggestProblemSize = std::stoi(leftToken);
  logDebug( "DatabaseEntry(std::string)", "got biggest problem size " << _biggestProblemSize );

  leftToken   = rightString.substr( 0, rightString.find(",") );
  rightString = rightString.substr( leftToken.size()+1 );
  _currentGrainSize = std::stoi(leftToken);
  logDebug( "DatabaseEntry(std::string)", "got current grain size " <<  _currentGrainSize );

  leftToken   = rightString.substr( 0, rightString.find(",") );
  rightString = rightString.substr( leftToken.size()+1 );
  _searchDelta = std::stoi(leftToken);
  logDebug( "DatabaseEntry(std::string)", "search delta is " <<  _searchDelta );

  leftToken   = rightString.substr( 0, rightString.find(",") );
  rightString = rightString.substr( leftToken.size()+1 );
  double serialAccuracy = std::stof(leftToken);
  logDebug( "DatabaseEntry(std::string)", "serial accuracy is " <<  serialAccuracy << ". Ignore remainder of this line");

  leftToken   = rightString.substr( 0, rightString.find(",") );
  rightString = rightString.substr( leftToken.size()+1 );
  double parallelAccuracy = std::stof(leftToken);
  logDebug( "DatabaseEntry(std::string)", "parallel accuracy is " <<  parallelAccuracy << ". Ignore remainder of this line");

  leftToken   = rightString.substr( 0, rightString.find(",") );
  rightString = rightString.substr( leftToken.size()+1 );
  _previousSpeedup = std::stof(leftToken);
  logDebug( "DatabaseEntry(std::string)", "previous speedup is " <<  _previousSpeedup << ". Ignore remainder of this line");

  _currentSerialMeasurement.erase();
  _currentSerialMeasurement.setAccuracy( serialAccuracy );

  _currentParallelMeasurement.erase();
  _currentParallelMeasurement.setAccuracy( parallelAccuracy );

  bool isValid = _biggestProblemSize>0
              && _currentGrainSize>0
              && _currentGrainSize<=_biggestProblemSize
              && _searchDelta>=0
	            && serialAccuracy>=0.0
              && parallelAccuracy>=0.0
              && _previousSpeedup>0.0
	            ;

  if (!isValid) {
    if (_biggestProblemSize<=0) {
      _biggestProblemSize = 65536;
      logError( "DatabaseEntry(std::string)", "unable to parse file entries w.r.t. biggest problem size" );
    }
    logError( "DatabaseEntry(std::string)", "input file seems to have been corrupted. Restart entry " << toString() << ". Corrupted line: " << inputLine );
    _currentGrainSize           = _biggestProblemSize;
    _searchDelta                = 0;
    _currentSerialMeasurement   = tarch::timing::Measurement( 0.0 );
    _currentParallelMeasurement = tarch::timing::Measurement( 0.0 );
    _previousSpeedup            = 1.0;
    restart();
  }
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::loadStatistics(const std::string& filename, int oracleNumber) {
  std::ifstream file(filename);
  std::string str = "";

  bool        tagOpen = false;
  while (std::getline(file, str)) {
    tagOpen &= str.compare("end OracleForOnePhaseWithShrinkingGrainSize")!=0;

    if (tagOpen) {
      #ifdef __EXCEPTIONS
      try {
      #endif

      std::string leftToken = "";
      std::string rightString = str;

      leftToken   = rightString.substr( 0, rightString.find("=") );
      rightString = rightString.substr( leftToken.size()+1 );
      peano::datatraversal::autotuning::MethodTrace  methodTrace = peano::datatraversal::autotuning::toMethodTrace(leftToken);
      logDebug( "loadStatistics(...)", "parse properties for " << toString(methodTrace) );

      DatabaseEntry newEntry(rightString);

      _measurements[methodTrace].push_back( newEntry );

      logDebug( "loadStatistics(...)", "added " << newEntry.toString() << " for " << toString(methodTrace) );
      #ifdef __EXCEPTIONS
      }
      catch (std::out_of_range& exception) {
        logError(
          "loadStatistics(...)",
          "failed to parse shared memory configuration file " << filename << " with error in " << exception.what() << " in adapter " << oracleNumber
        );
        logError(
          "loadStatistics(...)",
          "flawed string: " << str
        );
      }
      #endif
    }

    // Older GCC versions require an explicit cast here
    tagOpen |= str.compare( "adapter-number=" + std::to_string( (long long)oracleNumber) )==0;
  }

  assertion( _measurements.count(peano::datatraversal::autotuning::MethodTrace::NumberOfDifferentMethodsCalling)==0 );
}


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::~OracleForOnePhaseWithShrinkingGrainSize() {
}


peano::datatraversal::autotuning::OracleForOnePhase* sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::createNewOracle() const {
  return new OracleForOnePhaseWithShrinkingGrainSize(_learn,_restart,_selectNextStudiedMeasureTraceStrategy);
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
    ) {
      bool          oneMeasurementDidTerminateSearch = false;
      DatabaseEntry learningEntry(1);
      for (auto& p: _measurements[_activeMethodTrace]) {
        if ( p.isAccurate() && p.isSearching() ) {
          logInfo( "activateOracle()", "found entry that should learn for trace " << toString(_activeMethodTrace) );
          p.learn();
          if ( !p.isSearching() ) {
            logInfo( "activateOracle()", "entry fixed, so propagate its data: " << p.toString() );
            oneMeasurementDidTerminateSearch = true;
            learningEntry                    = p;
          }
        }
        else if (
          oneMeasurementDidTerminateSearch
          &&
          ( p.isSearching() || !p.isScaling() )
        ) {
          // propagate data
          p = DatabaseEntry( learningEntry, p.getBiggestProblemSize() );
          logInfo( "activateOracle()", "have propagated solution from " << learningEntry.toString() << " into " << p.toString() );
        }
      }
    }

    if (!_measurements.empty() ) {
      changeMeasuredMethodTrace();
      widenAccuracyOfCurrentlyStudiedMethodTraceAndRandomlyRestart();
    }
  }
}



bool sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::hasLearnedSinceLastQuery() {
  bool result = _hasLearnedSinceLastQuery;
  _hasLearnedSinceLastQuery = false;
  return result;
}
