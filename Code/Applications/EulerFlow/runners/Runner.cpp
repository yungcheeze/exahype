#include "EulerFlow/runners/Runner.h"

#include "EulerFlow/repositories/Repository.h"
#include "EulerFlow/repositories/RepositoryFactory.h"

#include "peano/utils/UserInterface.h"

#include "tarch/Assertions.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/autotuning/OracleForOnePhaseDummy.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"

#include "tarch/multicore/MulticoreDefinitions.h"

// @todo Remove this include as soon as you've created your real-world geometry
#include "peano/geometry/Hexahedron.h"

// ! Begin of code for DG method
#include "EulerFlow/Constants.h"
// ! End of code for DG method

tarch::logging::Log exahype::runners::Runner::_log("exahype::runners::Runner");

exahype::runners::Runner::Runner() {
  // @todo Insert your code here
}

exahype::runners::Runner::~Runner() {
  // @todo Insert your code here
}

int exahype::runners::Runner::run() {
  // @todo Insert your geometry generation here and adopt the repository
  //       generation to your needs. There is a dummy implementation to allow
  //       for a quick start, but this is really very dummy (it generates
  //       solely a sphere computational domain and basically does nothing with
  //       it).

  // Start of dummy implementation
  peano::geometry::Hexahedron geometry(
      tarch::la::Vector<DIMENSIONS, double>(1.0),
      tarch::la::Vector<DIMENSIONS, double>(0.0));
  exahype::repositories::Repository* repository =
      exahype::repositories::RepositoryFactory::getInstance()
          .createWithSTDStackImplementation(
              geometry,
              tarch::la::Vector<DIMENSIONS, double>(1.0),  // domainSize,
              tarch::la::Vector<DIMENSIONS, double>(
                  0.0)  // computationalDomainOffset
              );
  // End of dummy implementation

  int result = 0;
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    result = runAsMaster(*repository);
  }
#ifdef Parallel
  else {
    result = runAsWorker(*repository);
  }
#endif

  delete repository;

  return result;
}

int exahype::runners::Runner::runAsMaster(
    exahype::repositories::Repository& repository) {
  peano::utils::UserInterface userInterface;
  userInterface.writeHeader();

#ifdef SharedMemoryParallelisation
  // We assume that the workload per cell is that big that we can set the
  // enterCell
  // grain size to 1 as well as the minimum grain size. All the other values
  // remain
  // the default values.

  peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
      new peano::datatraversal::autotuning::OracleForOnePhaseDummy(true));
#endif
  // The space-tree is initialised with 1 coarse grid cell on level 1 and 3^d
  // fine grid cells on level 2.
  for (int coarseGridLevel = 1;
       coarseGridLevel < EXAHYPE_INITIAL_GLOBAL_REFINEMENT_LEVEL;
       coarseGridLevel++) {
    repository.switchToInitialGrid();
    do {
      repository.iterate();
    } while (!repository.getState().isGridBalanced());
  }
  repository.switchToGridExport();  // export the grid
  repository.iterate();

  repository
      .switchToPatchInitialisation();  // initialise the cell descriptions;
  repository.iterate();

  /*
   * Apply the initial conditions.
   * Then, compute the the initial current time step size.
   */
  repository.switchToInitialConditionAndExportAndGlobalTimeStepComputation();
  repository.iterate();

  /*
   * Compute current first predictor based on current time step size.
   * Set current time step size as old time step size of next iteration.
   * Compute the current time step size of the next iteration.
   */
  repository.switchToPredictorAndGlobalTimeStepComputation();
  repository.iterate();
  logInfo("runAsMaster(...)", "[ExaHyPE] global time step="
                                  << 0 << ", global time step size (seconds)="
                                  << repository.getState().getTimeStepSize());
  logInfo("runAsMaster(...)",
          "[ExaHyPE] old global time step size (seconds)="
              << repository.getState().getOldTimeStepSize());

  for (int n = 1; n < 400; n++) {
    /*
     * Exchange the fluctuations.
     */
    repository.switchToFaceDataExchange();
    repository.iterate();

    /*
     * The two adapters that are embedded in the if clause below perform the
     *following steps:
     *
     * 1. Perform the corrector step using the old update and the old global
     *time step size.
     *    This is a leaf-cell-local operation. Thus we immediately obtain the
     *leaf-cell-local current solution.
     * 2. Perform the predictor step using the leaf-cell-local current solution
     *and the current global time step size.
     * 3. Compute the leaf-cell-local time step sizes
     * 4. (Optionally) Export the leaf-cell-local current solution.
     * 5. After the traversal, set the global current time step size as the new
     *old global time step size.
     *    Find the minimum leaf-cell-local time step size and set it as the new
     *current
     *    global time step size.
     */
    if (n % EXAHYPE_PLOTTING_STRIDE == 0) {
      repository
          .switchToCorrectorAndPredictorAndGlobalTimeStepComputationAndExport();
    } else {
      repository.switchToCorrectorAndPredictorAndGlobalTimeStepComputation();
    }
    repository.iterate();

    logInfo("runAsMaster(...)", "[ExaHyPE] global time step="
                                    << n + 1
                                    << ", global time step size (seconds)="
                                    << repository.getState().getTimeStepSize());
    logInfo("runAsMaster(...)",
            "[ExaHyPE] old global time step size (seconds)="
                << repository.getState().getOldTimeStepSize());
  }

  repository.logIterationStatistics();
  repository.terminate();
  // ! End of code for DG method

  return 0;
}
