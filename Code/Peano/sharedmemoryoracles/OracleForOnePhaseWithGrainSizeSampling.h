// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _SHARED_MEMORY_ORAClES_ORACLE_FOR_ONE_PHASE_WITH_GRAIN_SIZE_SAMPLING_H_
#define _SHARED_MEMORY_ORAClES_ORACLE_FOR_ONE_PHASE_WITH_GRAIN_SIZE_SAMPLING_H_


#include "tarch/logging/Log.h"
#include "peano/datatraversal/autotuning/OracleForOnePhase.h"
#include "tarch/timing/Measurement.h"


#include <map>
#include <vector>


namespace sharedmemoryoracles {
  class OracleForOnePhaseWithGrainSizeSampling;
}


/**
 * Grain Size Sampling
 *
 * This oracle samples the grain size stochastically and allows us to get a
 * first impression about the distribution of well-suited grain size choices.
 * In practice it probably should not be used - it will not give you a fast
 * code.
 *
 * <h1>Postprocessing</h1>
 *
 * Though the sampling writes text output (on demand), it is not that easy to
 * understand. With the oracle comes a Pylab script that creates a couple of
 * plots that visualise the sampled data. It accepts the sampling output file
 * as argument and gives a new file with the additional extension .html.
 *
 * adapters=6,phases=10  fehlt noch
 *
 * There is an example sampling file in this directory called
 * sampling.example-properties. It can be used to check what the output should
 * look like.
 *
 * @author Tobias Weinzierl
 */
class sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling: public peano::datatraversal::autotuning::OracleForOnePhase {
  private:
    static tarch::logging::Log  _log;

    /**
     * Map selected grain sizes to measurements.
     */
    typedef std::map< int, tarch::timing::Measurement >             ExecutionTimeDatabase;

    /**
     * Holds for each input problem size the corresponding set of samplings.
     */
    typedef std::map< int, ExecutionTimeDatabase >   ExecutionTimeSamplingDatabase;

    int                            _numberOfSamples;
    bool                           _logarithmicDistribution;
    ExecutionTimeSamplingDatabase  _executionTimes;
    tarch::timing::Measurement*    _currentMeasurement;

    /**
     * Only for statistics.
     */
    const int            _adapterNumber;

    /**
     * Only for statistics.
     */
    const peano::datatraversal::autotuning::MethodTrace    _methodTrace;

    const bool           _useThreadPipelining;
  public:
    /**
     * @param numberOfSamples The number of samples that the Oracle will use.
     * @param logarithmicDistribution Decides whether the samples will be distributed in an logarithmic (true) or
     * an equidistant (false) way.
     */
    OracleForOnePhaseWithGrainSizeSampling(int numberOfSamples, bool useThreadPipelining, bool logarithmicDistribution, int adapterNumber=-1, const peano::datatraversal::autotuning::MethodTrace& methodTrace = peano::datatraversal::autotuning::NumberOfDifferentMethodsCalling);

    virtual ~OracleForOnePhaseWithGrainSizeSampling();

    /**
     * Determine grain size to sample
     *
     * If the method trace is ascend or descend, the function either returns 1
     * or 0 as these operations just distinguish between two variants.
     * Otherwise, we either use a logarithmic or an equidistant sample
     * distribution. However, if we use equidistant and the problem size is
     * smaller than the number of samples, we need a different computing rule
     * (division of integers otherwise gives a spacing of 0 between different
     * sample values).
     */
    virtual std::pair<int,bool> parallelise(int problemSize);

    /**
     */
    virtual void parallelSectionHasTerminated(double elapsedCalendarTime);

    virtual void plotStatistics(std::ostream& out) const;
    virtual void loadStatistics(const std::string& filename);

    virtual void informAboutElapsedTimeOfLastTraversal(double elapsedTime);

    /**
     * For this oracle type, the adapter number is completely irrelevant.
     */
    virtual peano::datatraversal::autotuning::OracleForOnePhase* createNewOracle(int adapterNumber, const peano::datatraversal::autotuning::MethodTrace& methodTrace) const;
};


#endif
