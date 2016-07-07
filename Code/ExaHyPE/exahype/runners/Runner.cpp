/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#include "exahype/runners/Runner.h"
#include "exahype/repositories/Repository.h"
#include "exahype/repositories/RepositoryFactory.h"

#include "tarch/Assertions.h"

#include "tarch/logging/CommandLineLogger.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"

#include "tarch/multicore/Core.h"
#include "tarch/multicore/MulticoreDefinitions.h"

#include "peano/parallel/JoinDataBufferPool.h"
#include "peano/parallel/JoinDataBufferPool.h"
#include "peano/parallel/SendReceiveBufferPool.h"
#include "peano/parallel/loadbalancing/Oracle.h"
#include "peano/parallel/loadbalancing/OracleForOnePhaseWithGreedyPartitioning.h"
#include "tarch/parallel/FCFSNodePoolStrategy.h"

#include "peano/geometry/Hexahedron.h"
#include "peano/utils/UserInterface.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/autotuning/OracleForOnePhaseDummy.h"

#include "sharedmemoryoracles/OracleForOnePhaseWithGrainSizeSampling.h"
#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"

#include "exahype/plotters/Plotter.h"

#include "exahype/solvers/Solver.h"

tarch::logging::Log exahype::runners::Runner::_log("exahype::runners::Runner");

exahype::runners::Runner::Runner(const Parser& parser) : _parser(parser) {}

exahype::runners::Runner::~Runner() {}

void exahype::runners::Runner::initDistributedMemoryConfiguration() {
  // @todo evtl. fehlen hier die Includes
  /*
    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      tarch::parallel::NodePool::getInstance().setStrategy(
        new mpibalancing::FairNodePoolStrategy(6)
      );
    }
    #else
  */
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    tarch::parallel::NodePool::getInstance().setStrategy(
        new tarch::parallel::FCFSNodePoolStrategy());
  }

  // @todo evtl. fehlen hier die Includes
  /*
    peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
      new mpibalancing::StaticBalancing(true)
    );
    #else
  */
  peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
      new peano::parallel::loadbalancing::
          OracleForOnePhaseWithGreedyPartitioning(true));

  tarch::parallel::NodePool::getInstance().restart();
  tarch::parallel::NodePool::getInstance().waitForAllNodesToBecomeIdle();

#if defined(Debug) || defined(Asserts)
  tarch::parallel::Node::getInstance().setDeadlockTimeOut(120 * 4);
  tarch::parallel::Node::getInstance().setTimeOutWarning(60 * 4);
#else
  tarch::parallel::Node::getInstance().setDeadlockTimeOut(120);
  tarch::parallel::Node::getInstance().setTimeOutWarning(60);
#endif

  const int bufferSize = 64;
  peano::parallel::SendReceiveBufferPool::getInstance().setBufferSize(
      bufferSize);
  peano::parallel::JoinDataBufferPool::getInstance().setBufferSize(bufferSize);
}

void exahype::runners::Runner::shutdownDistributedMemoryConfiguration() {
  tarch::parallel::NodePool::getInstance().terminate();
  exahype::repositories::RepositoryFactory::getInstance()
      .shutdownAllParallelDatatypes();
}

void exahype::runners::Runner::initSharedMemoryConfiguration() {
#ifdef SharedMemoryParallelisation
  const int numberOfThreads = _parser.getNumberOfThreads();
  tarch::multicore::Core::getInstance().configure(numberOfThreads);

  switch (_parser.getMulticoreOracleType()) {
    case Parser::Dummy:
      logInfo("initSharedMemoryConfiguration()",
              "use dummy shared memory oracle");
      peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
          new peano::datatraversal::autotuning::OracleForOnePhaseDummy(
              true, false,
              0)  // @todo Vasco bitte mal auf 1 setzen und nochmal durchjagen
          );
      break;
    case Parser::Autotuning:
      logInfo("initSharedMemoryConfiguration()",
              "use autotuning shared memory oracle");
      peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
          new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize());
      break;
    case Parser::GrainSizeSampling:
      logInfo("initSharedMemoryConfiguration()",
              "use shared memory oracle sampling");
      peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
          new sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling(
              32,
              false,  // useThreadPipelining,
              true    // logarithmicDistribution
              ));
      break;
  }

  std::ifstream f(_parser.getMulticorePropertiesFile().c_str());
  if (f.good()) {
    peano::datatraversal::autotuning::Oracle::getInstance().loadStatistics(
        _parser.getMulticorePropertiesFile());
  }
  f.close();
#endif
}

void exahype::runners::Runner::shutdownSharedMemoryConfiguration() {
#ifdef SharedMemoryParallelisation
  switch (_parser.getMulticoreOracleType()) {
    case Parser::Dummy:
      break;
    case Parser::Autotuning:
    case Parser::GrainSizeSampling:
      logInfo("shutdownSharedMemoryConfiguration()",
              "wrote statistics into file "
                  << _parser.getMulticorePropertiesFile());
      peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics(
          _parser.getMulticorePropertiesFile());
      break;
  }
#endif
}


exahype::repositories::Repository* exahype::runners::Runner::createRepository() const {
  // Geometry is static as we need it survive the whole simulation time.
  static peano::geometry::Hexahedron geometry(
    _parser.getDomainSize(),
    tarch::la::Vector<DIMENSIONS, double>(_parser.getOffset()));

  logInfo(
    "run(...)",
    "create computational domain at " << _parser.getOffset() <<
    " of width/size " << _parser.getDomainSize() <<
    ". bounding box has size " << _parser.getBoundingBoxSize() );

  return exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
    geometry,
    tarch::la::Vector<DIMENSIONS, double>(_parser.getBoundingBoxSize()),
    tarch::la::Vector<DIMENSIONS, double>(_parser.getOffset()));
}


int exahype::runners::Runner::run() {

  exahype::repositories::Repository* repository = createRepository();

  initDistributedMemoryConfiguration();
  initSharedMemoryConfiguration();

  int result = 0;
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    result = runAsMaster(*repository);
  }
  #ifdef Parallel
  else {
    result = runAsWorker(*repository);
  }
  #endif

  shutdownSharedMemoryConfiguration();
  shutdownDistributedMemoryConfiguration();

  delete repository;

  return result;
}

int exahype::runners::Runner::runAsMaster(
    exahype::repositories::Repository& repository) {
  peano::utils::UserInterface::writeHeader();

  /*
   * Build up the initial space tree.
   */
  repository.switchToAugmentedAMRGrid();
  int gridSetupIterations = 0;
  do {
    repository.iterate();
    gridSetupIterations++;
  } while (!repository.getState().isGridBalanced());
  repository.iterate();  // We need one extra iteration.
  gridSetupIterations++;
  repository.iterate();  // We need one extra iteration.
  gridSetupIterations++;

  logInfo("runAsMaster()",
          "grid setup iterations=" << gridSetupIterations << ", max-level="
                                   << repository.getState().getMaxLevel());
#ifdef Parallel
  logInfo("runAsMaster()",
          "number of working ranks=" << tarch::parallel::NodePool::getInstance()
                                            .getNumberOfWorkingNodes());
  logInfo(
      "runAsMaster()",
      "number of idle ranks="
          << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes());
#endif
  repository.switchToSolutionUpdateAndGlobalTimeStepComputation();
  repository.iterate();
#if defined(Debug) || defined(Asserts)
  startNewTimeStep(-1,true);
#else
  startNewTimeStep(-1,false);
#endif
  /*
   * Set the time stamps of the solvers to the initial value again.
   *
   * !!! Rationale
   *
   * The time step size computation
   * sets the predictor time stamp to the value
   * of the predictor time stamp plus the admissible
   * time step size on each patch for each solver.
   */
  initSolverTimeStamps();

  /*
   * Compute current first predictor based on current time step size.
   * Set current time step size as old time step size of next iteration.
   * Compute the current time step size of the next iteration.
   */
  bool plot = exahype::plotters::isAPlotterActive(
      solvers::Solver::getMinSolverTimeStamp());
  if (plot) {
    repository.switchToPredictorAndPlotAndGlobalTimeStepComputation();
  }
  else {
    repository.switchToPredictorAndGlobalTimeStepComputation();
  }
  repository.iterate();
  startNewTimeStep(0,true);

  // #if DIMENSIONS==2
  //  repository.switchToPlotAugmentedAMRGrid();
  //  repository.iterate();
  // #endif

  const double simulationEndTime = _parser.getSimulationEndTime();
  int n = 1;

  while ((solvers::Solver::getMinSolverTimeStamp() < simulationEndTime) &&
         tarch::la::greater(solvers::Solver::getMinSolverTimeStepSize(), 0.0)) {
    bool plot = exahype::plotters::isAPlotterActive(
        solvers::Solver::getMinSolverTimeStamp());

    if (_parser.fuseAlgorithmicSteps()) {
      runOneTimeStampWithFusedAlgorithmicSteps(repository, plot);
      recomputePredictorIfNecessary(repository);
      startNewTimeStep(n,true);
    } else {
      runOneTimeStampWithFourSeparateAlgorithmicSteps(repository, plot);
      startNewTimeStep(n,true);
    }

    n++;
    logDebug("runAsMaster(...)", "state=" << repository.getState().toString());
  }

  repository.logIterationStatistics(true);
  repository.terminate();

  return 0;
}

void exahype::runners::Runner::initSolverTimeStamps() {
  // todo 16/02/26:Dominic Etienne Charrier: The initial time stamp
  // should be set by the user in his solver sub class.
  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    // todo:16/03/04:Dominic Charrier
    p->setMinPredictorTimeStamp(
        0.0);  // introduce reset method that sets both to t=zero
    p->setMinCorrectorTimeStamp(
        0.0);  // introduce reset method that sets both to t=zero
  }
}

void exahype::runners::Runner::startNewTimeStep(int n,bool printInfo) {
  double currentMinTimeStamp = std::numeric_limits<double>::max();
  double currentMinTimeStepSize = std::numeric_limits<double>::max();
  double nextMinTimeStepSize = std::numeric_limits<double>::max();

  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    p->startNewTimeStep();

    currentMinTimeStamp =
        std::min(currentMinTimeStamp, p->getMinCorrectorTimeStamp());
    currentMinTimeStepSize =
        std::min(currentMinTimeStepSize, p->getMinCorrectorTimeStepSize());
    nextMinTimeStepSize =
        std::min(nextMinTimeStepSize, p->getMinPredictorTimeStepSize());
  }

  if (printInfo) {
    logInfo("startNewTimeStep(...)",
            "step " << n << "\t t_min          =" << currentMinTimeStamp);

    logInfo("startNewTimeStep(...)",
            "\t\t dt_min         =" << currentMinTimeStepSize);

    logInfo("startNewTimeStep(...)",
            "\t\t next dt_min    =" << nextMinTimeStepSize);
  }
#if defined(Debug) || defined(Asserts)
  tarch::logging::CommandLineLogger::getInstance()
      .closeOutputStreamAndReopenNewOne();
#endif
}

void exahype::runners::Runner::runOneTimeStampWithFusedAlgorithmicSteps(
    exahype::repositories::Repository& repository, bool plot) {
  /*
   * The adapter below performs the following steps:
   *
   * 1. Exchange the fluctuations using the predictor computed in the previous
   *sweep
   *    and the corrector time stemp size.
   * 2. Perform the corrector step using the corrector update and the corrector
   *time step size.
   *    This is a cell-local operation. Thus we immediately obtain the
   *cell-local current solution.
   * 3. Perform the predictor step using the cell-local current solution and the
   *predictor time step size.
   * 4. Compute the cell-local time step sizes
   */

  if (plot) {
    repository.switchToADERDGTimeStepAndPlot();
  } else {
    repository.switchToADERDGTimeStep();
  }

  repository.iterate();
}

// @todo 16/02/29:Dominic Etienne Charrier
// @Tobias: This should move into solver class, or not?
// The function does only make sense for optimistic time stepping
bool exahype::runners::Runner::
    setAccurateTimeStepSizesIfStabilityConditionWasHarmed() {
  bool cflConditionWasViolated = false;

  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    bool solverTimeStepSizeIsInstable = (p->getMinPredictorTimeStepSize() >
                                         p->getMinNextPredictorTimeStepSize());

    // todo 16/02/26:Dominic Etienne Charrier: The initial time stamp
    // introduce reset method that sets both to t=0.99*... make alpha=0.99
    // solver variable
    if (solverTimeStepSizeIsInstable) {
      p->updateMinNextPredictorTimeStepSize(
          0.99 * p->getMinNextPredictorTimeStepSize());  // set next
                                                         // predictor time
                                                         // step size
      p->setMinPredictorTimeStepSize(
          0.99 * p->getMinPredictorTimeStepSize());  // set next corrector
                                                     // time step size
    } else {
      p->updateMinNextPredictorTimeStepSize(
          .5 * (p->getMinPredictorTimeStepSize() +
                p->getMinNextPredictorTimeStepSize()));
    }

    cflConditionWasViolated =
        cflConditionWasViolated | solverTimeStepSizeIsInstable;
  }

  return cflConditionWasViolated;  // | tooDiffusive;
}

void exahype::runners::Runner::recomputePredictorIfNecessary(
    exahype::repositories::Repository& repository) {
  // Must be evaluated before we start a new time step
  bool stabilityConditionWasHarmed =
      setAccurateTimeStepSizesIfStabilityConditionWasHarmed();
  // Note that it is important to switch the time step sizes, i.e,
  // start a new time step, before we recompute the predictor.

  if (stabilityConditionWasHarmed) {
    logInfo("startNewTimeStep(...)",
            "\t\t Space-time predictor must be recomputed.");

    repository.switchToPredictorRerun();
    repository.iterate();
  }
}

void exahype::runners::Runner::runOneTimeStampWithFourSeparateAlgorithmicSteps(
    exahype::repositories::Repository& repository, bool plot) {
  // Only one time step (predictor vs. corrector) is used in this case.
  repository.switchToFaceDataExchange();  // Riemann -> face2face
  repository.iterate();
  repository.switchToCorrector();  // Face to cell
  repository.iterate();

  if (plot) {
    repository.switchToPlotAndGlobalTimeStepComputation();  // Inside cell
  } else {
    repository.switchToGlobalTimeStepComputation();  // Inside cell
  }
  repository.iterate();

  repository.switchToPredictor();  // Cell onto faces
  repository.iterate();
}
