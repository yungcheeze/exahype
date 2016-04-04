// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _SHARED_MEMORY_ORAClES_ORACLE_FOR_ONE_PHASE_WITH_SHRINKING_GRAIN_SIZE_H_
#define _SHARED_MEMORY_ORAClES_ORACLE_FOR_ONE_PHASE_WITH_SHRINKING_GRAIN_SIZE_H_


#include "tarch/logging/Log.h"
#include "peano/datatraversal/autotuning/OracleForOnePhase.h"
#include "tarch/timing/Measurement.h"


#include <map>


namespace sharedmemoryoracles {
  class OracleForOnePhaseWithShrinkingGrainSize;
}


/**
 * Oracle With Shrinking Grain Size
 *
 * This is a very simple oracle that runs through three different states:
 *
 * - First, it tries to find out a reasonable time per unknown in serial mode.
 * - Second, it sets the optimal grain size to half of the biggest grain size
 *   found so far, and halves this grain size iteratively as long as the
 *   runtime improves.
 * - If the grain size reduction leads to an increase in runtime, the grain size
 *   is doubled and this is the optimal grain size returned on each request.
 *
 * With each grain size modification, the oracle also makes the measurements
 * more precise, i.e. it waits longer for valid data.
 *
 * The fundamental assumption for the runtime here is that the speedup over
 * grain size is something like a parabula, i.e. has an unique maximum somewhere
 * between 0 and the total grain size. We assume that a scaling code part
 * performs better if we set the grain size to half of the problem size and then
 * we try to approximate the best grain size iteratively from the right, i.e. our
 * best grain size decreases monotonically.
 *
 * This oracle relies on a workload per cell that is more or less constant. The
 * oracle does not work for particle methods, e.g., where the workload does not
 * directly correlate to the number of cells.
 *
 * @author Tobias Weinzierl
 */
class sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize: public peano::datatraversal::autotuning::OracleForOnePhase {
  private:
    static tarch::logging::Log                           _log;

    /**
     * We never do optimise all traces. We only do it with one trace at a time.
     */
    static int                                           _activeMethodTrace;
    static int                                           _delayBetweenTwoUpdates;

    const peano::datatraversal::autotuning::MethodTrace  _methodTrace;
    const int                                            _adapterNumber;

    bool                                                 _oracleIsSearching;
    int                                                  _biggestProblemSize;
    int                                                  _currentGrainSize;
    tarch::timing::Measurement                           _currentMeasurement;
    double                                               _previousMeasuredTime;
    double                                               _lastProblemSize;

    void makeAttributesLearn();

    /**
     * We do only measure one method trace at a time.
     */
    void changeMeasuredMethodTrace();
  public:
    /**
     * Oracle Constructor
     */
    OracleForOnePhaseWithShrinkingGrainSize(int adapterNumber=-1, const peano::datatraversal::autotuning::MethodTrace& methodTrace = peano::datatraversal::autotuning::NumberOfDifferentMethodsCalling);

    virtual ~OracleForOnePhaseWithShrinkingGrainSize();

    virtual std::pair<int,bool> parallelise(int problemSize);
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
