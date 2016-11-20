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
 * <h2> Usage </h2>
 *
 * The usuage per se is trivial: simple pass an instance of this class to the
 * oracle:
 *
 * <pre>
peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
//    new peano::datatraversal::autotuning::OracleForOnePhaseDummy(true)
  new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize()
);
   </pre>
 *
 * This solution however has one major flaw: The oracle starts to search for
 * good grain sizes everytime even if there have good grain sizes been
 * identified already. Notably, it does not store insight persistently
 * in-between two program runs.
 *
 * To allow the oracle to store its data persistently and to reload properties
 * from older runs, you have to use the oracle's load and plotStatistics
 * commands. Just after you have invoked setOracle(), add a statement alike
 * <pre>
  peano::datatraversal::autotuning::Oracle::getInstance().loadStatistics( "sharedmemory.grain-sizes" );
 </pre>
 *
 * In return, you have to store your data before you quit the program.
 * plotStatistics() allows you to pipe all information directly into a file.
 *
 * <pre>
  peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics( "sharedmemory.grain-sizes" );
 </pre>
 *
 *
 * If you apply your code to multiple problems, it might be reasonable to run
 * the autotuning for different characteristic problem sets. Often, a good
 * grain size depends on the input parameters, too.
 *
 * <h2> Analysing the output </h2>
 *
 * The output file (sharedmemory.gain-sizes in the above example) is a plain
 * text file. It is however non-trivial to read. Therefore, I provide a Python
 * script in this directory to analyse the files. Just invoke
 *
 *
 * <h2> MPI+X </h2>
 *
 *
 * @author Tobias Weinzierl
 */
class sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize: public peano::datatraversal::autotuning::OracleForOnePhase {
  private:
    static tarch::logging::Log                           _log;
    static const double                                  _InitialAccuracy;
    static const double                                  _FinalAccuracy;

    /**
     * We never do optimise all traces. We only do it with one trace at a time.
     */
    peano::datatraversal::autotuning::MethodTrace        _activeMethodTrace;

    struct DatabaseEntry {
      int                          _biggestProblemSize;
      int                          _currentGrainSize;
      tarch::timing::Measurement   _currentMeasurement;
      int                          _previousGrainSize;
      double                       _previousMeasuredTime;
    };


    /**
     * Memorise for each pair of adapter number and trace number whether this
     * configuration actually has been used. We do sample only those
     * combinations that actually are used.
     */
    std::map<peano::datatraversal::autotuning::MethodTrace,DatabaseEntry>               _measurements;

    void makeAttributesLearn();

    void changeMeasuredMethodTrace();
  public:
    /**
     * Oracle Constructor
     */
    OracleForOnePhaseWithShrinkingGrainSize();

    virtual ~OracleForOnePhaseWithShrinkingGrainSize();

    peano::datatraversal::autotuning::GrainSize parallelise(int problemSize, peano::datatraversal::autotuning::MethodTrace askingMethod) override;
    void parallelSectionHasTerminated(int problemSize, peano::datatraversal::autotuning::MethodTrace askingMethod, double costPerProblemElement) override;
    void plotStatistics(std::ostream& out, int oracleNumber) const override;
    void loadStatistics(const std::string& filename, int oracleNumber) override;

    /**
     * Nop
     */
    void deactivateOracle() override;
    void activateOracle() override;
    peano::datatraversal::autotuning::OracleForOnePhase* createNewOracle() const override;
};


#endif
