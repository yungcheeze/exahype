#include "exahype/runners/Runner.h"
#include "exahype/repositories/Repository.h"
#include "exahype/repositories/RepositoryFactory.h"


#include "tarch/Assertions.h"
#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "tarch/multicore/Core.h"
#include "tarch/multicore/MulticoreDefinitions.h"

#include "peano/utils/UserInterface.h"
#include "peano/geometry/Hexahedron.h"


#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/autotuning/OracleForOnePhaseDummy.h"

#include "sharedmemoryoracles/OracleForOnePhaseWithGrainSizeSampling.h"

#include "exahype/plotters/Plotter.h"


tarch::logging::Log  exahype::runners::Runner::_log( "exahype::runners::Runner" );



exahype::runners::Runner::Runner(const Parser& parser):
  _parser(parser) {
}


exahype::runners::Runner::~Runner() {
}


void exahype::runners::Runner::initSharedMemoryConfiguration() {
#ifdef SharedMemoryParallelisation
  const int numberOfThreads = _parser.getNumberOfThreads();
  tarch::multicore::Core::getInstance().configure(numberOfThreads);

  std::ifstream f(_parser.getMulticorePropertiesFile().c_str());
  bool multicorePropertiesFileDoesExist = f.good();
  f.close();

  switch (_parser.getMulticoreOracleType()) {
    case Parser::Dummy:
      logInfo( "initSharedMemoryConfiguration()", "use dummy shared memory oracle" );
      peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
          new peano::datatraversal::autotuning::OracleForOnePhaseDummy(true)
      );
      break;
    case Parser::Autotuning:
      logInfo( "initSharedMemoryConfiguration()", "use autotuning shared memory oracle" );
      break;
    case Parser::GrainSizeSampling:
      logInfo( "initSharedMemoryConfiguration()", "use shared memory oracle sampling" );
      peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
          new sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling(
              10,
              true, // useThreadPipelining,
              true  // logarithmicDistribution
          )
      );
      break;
  }
#endif
}


void exahype::runners::Runner::shutdownSharedMemoryConfiguration() {
#ifdef SharedMemoryParallelisation
  switch (_parser.getMulticoreOracleType()) {
    case Parser::Dummy:
      break;
    case Parser::Autotuning:
      break;
    case Parser::GrainSizeSampling:
      logInfo( "shutdownSharedMemoryConfiguration()", "wrote statistics into file " << _parser.getMulticorePropertiesFile() );
      peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics(_parser.getMulticorePropertiesFile());
      break;
  }
#endif
}



int exahype::runners::Runner::run() {
  initSharedMemoryConfiguration();

  logInfo("run(...)", "create computational domain at " << _parser.getOffset() << " of width/size " << _parser.getSize() );

  peano::geometry::Hexahedron geometry(
    _parser.getSize(),
    tarch::la::Vector<DIMENSIONS,double>( _parser.getOffset() )
  );

  exahype::repositories::Repository* repository = 
    exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
      geometry,
      tarch::la::Vector<DIMENSIONS,double>( _parser.getSize() ),
      tarch::la::Vector<DIMENSIONS,double>( _parser.getOffset() )
    );

  int result = 0;
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    result = runAsMaster( *repository );
  }
#ifdef Parallel
  else {
    result = runAsWorker( *repository );
  }
#endif

  shutdownSharedMemoryConfiguration();

  delete repository;

  return result;
}


int exahype::runners::Runner::runAsMaster(exahype::repositories::Repository& repository) {
  peano::utils::UserInterface userInterface;
  userInterface.writeHeader();

  // The space-tree is initialised with 1 coarse grid cell on level 1 and 3^d fine grid cells on level 2.
  repository.switchToInitialGrid();
  do {
    repository.iterate();
  } while (!repository.getState().isGridBalanced());

  repository.switchToPatchInitialisation();       // initialise the cell descriptions;
  repository.iterate();

  /*
   * Apply the initial conditions.
   * Then, compute the the initial current time step size.
   */
  repository.getState().updateTimeStamp(0.0);
  repository.switchToInitialConditionAndGlobalTimeStepComputation();
  repository.iterate();

  logInfo(
    "runAsMaster(...)",
    "\t t_min     ="  << repository.getState().getMinimalGlobalTimeStamp() <<
    "\t dt_max    =" << repository.getState().getMaxTimeStepSize() <<
    "\t dt_max_old=" << repository.getState().getOldMaxTimeStepSize()
  );

  /*
   * Compute current first predictor based on current time step size.
   * Set current time step size as old time step size of next iteration.
   * Compute the current time step size of the next iteration.
   */
  // todo Dominic Etienne Charrier
  // I changed code here
  repository.switchToPredictorAndGlobalTimeStepComputation();
  repository.iterate();

  logInfo(
    "runAsMaster(...)",
    "\t t_min     ="  << repository.getState().getMinimalGlobalTimeStamp() <<
    "\t dt_max    =" << repository.getState().getMaxTimeStepSize() <<
    "\t dt_max_old=" << repository.getState().getOldMaxTimeStepSize()
  );

  const double simulationEndTime = _parser.getSimulationEndTime();
  int n=1;

  while ( repository.getState().getMinimalGlobalTimeStamp()<simulationEndTime ) {
    if (exahype::plotters::isAPlotterActive(repository.getState().getMinimalGlobalTimeStamp())) {
      repository.switchToPlot();
      repository.iterate();
      exahype::plotters::finishedPlotting();
      logInfo( "runAsMaster(...)", "all snapshots written" );
    }

    // todo Dominic Etienne Charrier
    // Swap of time steps has to be performed by adapter
    // since it already computes the new time step size.
    if ( _parser.fuseAlgorithmicSteps() ) {
      runOneTimeStampWithFusedAlgorithmicSteps(repository);
    }
    else {
      runOneTimeStampWithFourSeparateAlgorithmicSteps(repository);
    }

    logInfo(
      "runAsMaster(...)",
      "step " << n <<
      "\t t_min     ="  << repository.getState().getMinimalGlobalTimeStamp() <<
      "\t dt_max    =" << repository.getState().getMaxTimeStepSize() <<
      "\t dt_max_old=" << repository.getState().getOldMaxTimeStepSize()
    );

    n++;
    #if defined(Debug) || defined(Asserts)
    logInfo( "runAsMaster(...)", "state=" << repository.getState().toString() );
    #endif
  }

  repository.logIterationStatistics();
  repository.terminate();
  // ! End of code for DG method

  return 0;
}


void exahype::runners::Runner::runOneTimeStampWithFusedAlgorithmicSteps(exahype::repositories::Repository& repository) {
  /*
   * Exchange the fluctuations.
   */
  repository.switchToFaceDataExchange();
  repository.iterate();

  /*
   * The two adapters that are embedded in the if clause below perform the following steps:
   *
   * 1. Perform the corrector step using the old update and the old global time step size.
   *    This is a cell-local operation. Thus we immediately obtain the cell-local current solution.
   * 2. Perform the predictor step using the cell-local current solution and the current global time step size.
   * 3. Compute the cell-local time step sizes
   * 4. (Optionally) Export the cell-local current solution.
   * 5. After the traversal, set the global current time step size as the new old global time step size.
   *    Find the minimum cell-local time step size and set it as the new current
   *    global time step size.
   */
  repository.switchToCorrectorAndPredictorAndGlobalTimeStepComputation();
  repository.iterate();
}


void exahype::runners::Runner::runOneTimeStampWithFourSeparateAlgorithmicSteps(exahype::repositories::Repository& repository) {
  assertionMsg( false, "not implemented yet. Dominic, can you please do it?" );

  // todo Dominic Etienne Charrier
  // Only one time step (group) is used in this case.
}
