// This file originally was created by pdt (Peano Development Toolkit) as part
// of a code based upon the Peano project by Tobias Weinzierl. For conditions
// of distribution and use of this project, please see the copyright notice at
// www.peano-framework.org. Feel free to adopt the license and authorship of
// this file and your project to your needs as long as the license is in
// agreement with the original Peano user constraints. A reference to/citation
// of  Peano and its author is highly appreciated.
#ifndef _EXAHYPE_RUNNERS_RUNNER_H_
#define _EXAHYPE_RUNNERS_RUNNER_H_

#include "exahype/Parser.h"
#include "tarch/logging/Log.h"

#include "exahype/State.h"

namespace exahype {
namespace runners {
class Runner;
}
namespace repositories {
class Repository;
}
}

/**
 * Runner
 *
 */
class exahype::runners::Runner {
 private:
  static tarch::logging::Log _log;

  const exahype::Parser& _parser;

  /**
   * Setup the oracles for the shared memory parallelisation. Different
   * oracles can be employed:
   *
   * - If no autotuning is used and no valid properties file is provided and
   *   the code is compiled with -DPerformanceAnalysis, we use the grain size
   *   sampling
   * - If no autotuning is used and no valid properties file is provided, we
   *   use the default oracle coming along with the Peano kernel
   * - If autotuning is enabled and no valid properties file is provided, we
   *
   *
   * <h2>Invocation sequence</h2>
   *
   * It is important that we init the shared memory environment after we have
   * created the repository. See Orace::loadStatistics().
   */
  void initSharedMemoryConfiguration();

  /**
   * The shared memory environment has to be set up before we create the
   * repository.
   */
  void initDistributedMemoryConfiguration();
  void shutdownSharedMemoryConfiguration();
  void shutdownDistributedMemoryConfiguration();

  int runAsMaster(exahype::repositories::Repository& repository);

#ifdef Parallel
  int runAsWorker(exahype::repositories::Repository& repository);

  /**
   * If the master node calls runGlobalStep() on the repository, all MPI
   * ranks besides rank 0 invoke this operation no matter whether they are
   * idle or not. Hence, please call this operation manually within
   * runAsMaster() if you require the master to execute the same function
   * as well.
   */
  void runGlobalStep();
#endif

  void initSolvers();
  void startNewTimeStep(int n);

  /**
   * Do one time step where all phases are actually fused into one traversal
   *
   * @param plot      Do plot the way along
   */
  void runOneTimeStampWithFusedAlgorithmicSteps(
      exahype::repositories::Repository& repository, bool plot);

  bool setAccurateTimeStepSizesIfStabilityConditionWasHarmed();

  void recomputePredictorIfNecessary(
      exahype::repositories::Repository& repository);

  /**
   * Do one time step but actually use a couple of iterations to do so.
   *
   *
   * @param plot      Do plot in the after the corrector has been applied
   */
  void runOneTimeStampWithFourSeparateAlgorithmicSteps(
      exahype::repositories::Repository& repository, bool plot);

 public:
  Runner(const Parser& parser);
  virtual ~Runner();

  /**
   * Run
   */
  int run();
};

#endif
