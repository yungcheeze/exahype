#include "sharedmemoryoracles/OracleForOnePhaseWithGrainSizeSampling.h"
#include "tarch/Assertions.h"
#include "peano/utils/Globals.h"

#include <cstdlib>
#include <math.h>
#include <fstream>


tarch::logging::Log  sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling::_log( "sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling" );


sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling::OracleForOnePhaseWithGrainSizeSampling(
  int                 numberOfSamples,
  bool                useThreadPipelining,
  bool                logarithmicDistribution,
  int                 adapterNumber,
  const peano::datatraversal::autotuning::MethodTrace&  methodTrace
):
  _numberOfSamples(numberOfSamples),
  _logarithmicDistribution(logarithmicDistribution),
  _executionTimes(),
  _currentMeasurement(0),
  _adapterNumber(adapterNumber),
  _methodTrace(methodTrace),
  _useThreadPipelining(useThreadPipelining) {
  assertion( _numberOfSamples>0 );
}


sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling::~OracleForOnePhaseWithGrainSizeSampling() {
}


std::pair<int,bool> sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling::parallelise(int problemSize) {
  logTraceInWith1Argument( "parallelise(...)", problemSize );

  assertion( _currentMeasurement==0 );

  int result;

  if (
    !_useThreadPipelining &&
    (
      _methodTrace==peano::datatraversal::autotuning::PipelineAscendTask ||
      _methodTrace==peano::datatraversal::autotuning::PipelineDescendTask
    )
  ) {
    result = 0;
  }
  else if(_logarithmicDistribution){
    result = static_cast<int>(pow(pow(0.5, 0.25), rand() % _numberOfSamples) * problemSize);
  } else {
    if (problemSize<_numberOfSamples) {
      result = rand() % problemSize;
    }
    else {
      result = (problemSize/_numberOfSamples) * (rand() % _numberOfSamples);
    }
  }

  _executionTimes[problemSize][result].setAccuracy( std::numeric_limits<double>::max() );
  _currentMeasurement = &(_executionTimes[problemSize][result]);

  logTraceOutWith1Argument( "parallelise(...)", result );

  return std::pair<int,bool>(result,true);
}


void sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling::parallelSectionHasTerminated(double elapsedCalendarTime) {
  assertion( _currentMeasurement!=0 );

  _currentMeasurement->setValue(elapsedCalendarTime);
  _currentMeasurement = 0;
}



void sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling::loadStatistics(const std::string& filename) {
  assertionMsg( false, "not yet implemented" );
}


void sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling::plotStatistics(std::ostream& out) const {
  if (_executionTimes.empty()) return;

  out << "adapter=" << _adapterNumber-peano::datatraversal::autotuning::NumberOfPredefinedAdapters
      << ", phase="  << _methodTrace
      << ", name="   << peano::datatraversal::autotuning::toString(_methodTrace)
      << ", no-of-problem-sizes-studied=" << _executionTimes.size();

  for (ExecutionTimeSamplingDatabase::const_iterator p=_executionTimes.begin(); p!=_executionTimes.end(); p++) {
    out << ", problem-size=" << p->first
        << ", no-of-grain-sizes-studied=" << p->second.size() << ":";
    for (
      ExecutionTimeDatabase::const_iterator measurements = p->second.begin();
      measurements != p->second.end();
      measurements++
    ) {
      out << "(grain-size=" << measurements->first << ", runtime=" << measurements->second.toString() << ")";
    }
  }
}


peano::datatraversal::autotuning::OracleForOnePhase* sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling::createNewOracle(int adapterNumber, const peano::datatraversal::autotuning::MethodTrace& methodTrace) const {
  return new OracleForOnePhaseWithGrainSizeSampling(_numberOfSamples, _useThreadPipelining, _logarithmicDistribution, adapterNumber, methodTrace);
}


void sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling::informAboutElapsedTimeOfLastTraversal(double elapsedTime) {
}
