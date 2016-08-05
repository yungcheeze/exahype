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
#include "tarch/parallel/FCFSNodePoolStrategy.h"

#include "tarch/multicore/Core.h"
#include "tarch/multicore/MulticoreDefinitions.h"

#include "peano/parallel/JoinDataBufferPool.h"
#include "peano/parallel/JoinDataBufferPool.h"
#include "peano/parallel/SendReceiveBufferPool.h"
#include "peano/parallel/loadbalancing/Oracle.h"
#include "peano/parallel/loadbalancing/OracleForOnePhaseWithGreedyPartitioning.h"

#include "peano/geometry/Hexahedron.h"

#include "peano/utils/UserInterface.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/autotuning/OracleForOnePhaseDummy.h"

#include "sharedmemoryoracles/OracleForOnePhaseWithGrainSizeSampling.h"
#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"

#include "mpibalancing/GreedyBalancing.h"
#include "mpibalancing/StaticBalancing.h"

#include "exahype/plotters/Plotter.h"

#include "exahype/solvers/ADERDGSolver.h"

tarch::logging::Log exahype::runners::Runner::_log("exahype::runners::Runner");

exahype::runners::Runner::Runner(const Parser& parser) : _parser(parser) {}

exahype::runners::Runner::~Runner() {}

void exahype::runners::Runner::initDistributedMemoryConfiguration() {
  #ifdef Parallel
  std::string configuration = _parser.getMPIConfiguration();
  if (_parser.getMPILoadBalancingType()==Parser::MPILoadBalancingType::Static) {
    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      if (configuration.find( "FCFS" )!=std::string::npos ) {
        tarch::parallel::NodePool::getInstance().setStrategy(
          new tarch::parallel::FCFSNodePoolStrategy()
        );
        logInfo("initDistributedMemoryConfiguration()", "load balancing relies on FCFS answering strategy");
      }
      // @todo evtl. fehlen hier die Includes
      /*
        if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
          tarch::parallel::NodePool::getInstance().setStrategy(
            new mpibalancing::FairNodePoolStrategy(6)
          );
        }
        #else
      */
      else {
        logError("initDistributedMemoryConfiguration()", "no valid load balancing answering strategy specified: use FCFS");
        tarch::parallel::NodePool::getInstance().setStrategy(
          new tarch::parallel::FCFSNodePoolStrategy()
        );
      }
    }

    if ( configuration.find( "greedy" )!=std::string::npos ) {
      logInfo("initDistributedMemoryConfiguration()", "use greedy load balancing without joins (mpibalancing/GreedyBalancing)");
      peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
        new mpibalancing::GreedyBalancing(1,3)
      );
    }
    else if ( configuration.find( "hotspot" )!=std::string::npos ) {
      logInfo("initDistributedMemoryConfiguration()", "use global hotspot elimination without joins (mpibalancing/StaticBalancing)");
      peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
        new mpibalancing::StaticBalancing(false)
      );
    }
    else {
      logError("initDistributedMemoryConfiguration()", "no valid load balancing configured. Use greedy load balancing without joins");
      peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
        new peano::parallel::loadbalancing::OracleForOnePhaseWithGreedyPartitioning(false)
      );
    }
  } // end of static load balancing


  tarch::parallel::NodePool::getInstance().restart();
  tarch::parallel::NodePool::getInstance().waitForAllNodesToBecomeIdle();

  tarch::parallel::Node::getInstance().setDeadlockTimeOut(_parser.getMPITimeOut());
  tarch::parallel::Node::getInstance().setTimeOutWarning(_parser.getMPITimeOut()/2);
  logInfo("initDistributedMemoryConfiguration()", "use MPI time out of " << _parser.getMPITimeOut() << " (warn after half the timeout span)");

  const int bufferSize = _parser.getMPIBufferSize();
  peano::parallel::SendReceiveBufferPool::getInstance().setBufferSize(bufferSize);
  peano::parallel::JoinDataBufferPool::getInstance().setBufferSize(bufferSize);
  logInfo("initDistributedMemoryConfiguration()", "use MPI buffer size of " << bufferSize);
  #endif
}


void exahype::runners::Runner::shutdownDistributedMemoryConfiguration() {
  #ifdef Parallel
  tarch::parallel::NodePool::getInstance().terminate();
  exahype::repositories::RepositoryFactory::getInstance().shutdownAllParallelDatatypes();
  #endif
}

void exahype::runners::Runner::initSharedMemoryConfiguration() {
#ifdef SharedMemoryParallelisation
  const int numberOfThreads = _parser.getNumberOfThreads();
  tarch::multicore::Core::getInstance().configure(numberOfThreads);

  switch (_parser.getMulticoreOracleType()) {
    case Parser::MulticoreOracleType::Dummy:
      logInfo("initSharedMemoryConfiguration()",
              "use dummy shared memory oracle");
      peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
          new peano::datatraversal::autotuning::OracleForOnePhaseDummy(
              true, false,
              0)  // @todo Vasco bitte mal auf 1 setzen und nochmal durchjagen
          );
      break;
    case Parser::MulticoreOracleType::Autotuning:
      logInfo("initSharedMemoryConfiguration()",
              "use autotuning shared memory oracle");
      peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
          new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize());
      break;
    case Parser::MulticoreOracleType::GrainSizeSampling:
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
    case Parser::MulticoreOracleType::Dummy:
      break;
    case Parser::MulticoreOracleType::Autotuning:
    case Parser::MulticoreOracleType::GrainSizeSampling:
      #ifdef Parallel
      if (tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getNumberOfNodes()-1) {
        logInfo("shutdownSharedMemoryConfiguration()",
          "wrote statistics into file " << _parser.getMulticorePropertiesFile()
	  << ". Dump from all other ranks subpressed to avoid file races"
	);
        peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics(
          _parser.getMulticorePropertiesFile());
      }
      #else
      logInfo("shutdownSharedMemoryConfiguration()",
              "wrote statistics into file "
                  << _parser.getMulticorePropertiesFile());
      peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics(
          _parser.getMulticorePropertiesFile());
      #endif
      break;
  }
#endif
}


exahype::repositories::Repository* exahype::runners::Runner::createRepository() const {
  // Geometry is static as we need it survive the whole simulation time.
  static peano::geometry::Hexahedron geometry(
    _parser.getDomainSize(),
    _parser.getOffset());

  logDebug(
    "run(...)",
    "create computational domain at " << _parser.getOffset() <<
    " of width/size " << _parser.getDomainSize() <<
    ". bounding box has size " << _parser.getBoundingBoxSize() );

  #ifdef Parallel
  return exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
    geometry,
    _parser.getBoundingBoxSize()*9.0/7.0,
    _parser.getOffset()-1.0/7.0*_parser.getBoundingBoxSize()
  );
  #else
  return exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
    geometry,
    _parser.getBoundingBoxSize(),
    _parser.getOffset()
  );
  #endif
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

void exahype::runners::Runner::createGrid(exahype::repositories::Repository& repository) {
  #ifdef Parallel
  const bool UseStationaryCriterion = tarch::parallel::Node::getInstance().getNumberOfNodes()==1;
  #else
  const bool UseStationaryCriterion = true;
  #endif

  int gridSetupIterations = 0;
  repository.switchToAugmentedAMRGrid();

  do {
    repository.iterate();
    gridSetupIterations++;
    repository.iterate();
    gridSetupIterations++;
    logInfo("runAsMaster()",
      "grid setup iteration #" << gridSetupIterations <<
      ", max-level=" << repository.getState().getMaxLevel() <<
      ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
  }
  while (
   ( UseStationaryCriterion && !repository.getState().isGridStationary())
   ||
   (!UseStationaryCriterion && !repository.getState().isGridBalanced())
  );

  logInfo("createGrid(Repository)", "finished grid setup after " << gridSetupIterations << " iterations" );

  if (tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()>0) {
    logWarning( "createGrid(Repository)", "there are still " << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes() << " ranks idle" )
  }
}


int exahype::runners::Runner::runAsMaster(exahype::repositories::Repository& repository) {
  peano::utils::UserInterface::writeHeader();

  initSolverTimeStamps();
  createGrid(repository);

  #if defined(Dim2) && defined(Asserts)
  repository.switchToPlotAugmentedAMRGrid();
  repository.iterate();
  #endif

  /*
   * Set ADER-DG corrector time stamp and finite volumes time stamp.
   * Compute ADER-DG corrector time step size implicitly and finite volumes time step size.
   * (Implicitly means here that we set the predictor time step size but after the next newTimeStep(...)
   * the corrector time step size is set as this predictor time step size.)
   */
  initSolverTimeStamps();
  repository.switchToSolutionAdjustmentAndGlobalTimeStepComputation();
  repository.iterate();
  #if defined(Debug) || defined(Asserts)
  startNewTimeStep(-1,true,true);
  #else
  startNewTimeStep(-1,true,false);
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
      solvers::Solver::getMinSolverTimeStampOfAllSolvers());
  if (plot) {
    repository.switchToPredictorAndPlotAndGlobalTimeStepComputation();
  }
  else {
    repository.switchToPredictorAndGlobalTimeStepComputation();
  }
  repository.iterate();

  /*
   * Set ADER-DG predictor time stamp and do not touch the finite volumes time stamp.
   * Compute ADER-DG predictor time step size and do not touch finite volumes time step size.
   */
  startNewTimeStep(0,false,true);

  const double simulationEndTime = _parser.getSimulationEndTime();
  int n = 1;

  logDebug("runAsMaster(...)","min solver time stamp: "     << solvers::Solver::getMinSolverTimeStampOfAllSolvers()); // change to log debug
  logDebug("runAsMaster(...)","min solver time step size: " << solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers());

  while ((solvers::Solver::getMinSolverTimeStampOfAllSolvers() < simulationEndTime) &&
         tarch::la::greater(solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers(), 0.0)) {
    bool plot = exahype::plotters::isAPlotterActive(
        solvers::Solver::getMinSolverTimeStampOfAllSolvers());

    if (_parser.getFuseAlgorithmicSteps()) {
      runOneTimeStampWithFusedAlgorithmicSteps(repository, plot);
      recomputePredictorIfNecessary(repository,_parser.getFuseAlgorithmicStepsFactor());
      startNewTimeStep(n,true,true);
    } else {
      runOneTimeStampWithFourSeparateAlgorithmicSteps(repository, plot);
      startNewTimeStep(n,true,true);
    }

    n++;
    logDebug("runAsMaster(...)", "state=" << repository.getState().toString());
  }
  if ( tarch::la::equals(solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers(), 0.0)) {
    logWarning("runAsMaster(...)","Minimum solver time step size is zero (up to machine precision).");
  }

  repository.logIterationStatistics(true);
  repository.terminate();

  return 0;
}

void exahype::runners::Runner::initSolverTimeStamps() {
  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    p->initInitialTimeStamp(0.0);
  }
}

void exahype::runners::Runner::initFiniteVolumesSolverTimeStamps() {
  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    if (p->getType()==exahype::solvers::Solver::Type::FiniteVolumes) {
      p->initInitialTimeStamp(0.0);
    }
  }
}

void exahype::runners::Runner::startNewTimeStep(int n,bool startNewFiniteVolumesTimeStep,bool printInfo) {
  double currentMinTimeStamp = std::numeric_limits<double>::max();
  double currentMinTimeStepSize = std::numeric_limits<double>::max();
  double nextMinTimeStepSize = std::numeric_limits<double>::max();

  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    if (startNewFiniteVolumesTimeStep || p->getType()==exahype::solvers::Solver::Type::ADER_DG) {
      p->startNewTimeStep();
    }
    currentMinTimeStamp =
        std::min(currentMinTimeStamp, p->getMinTimeStamp());
    currentMinTimeStepSize =
        std::min(currentMinTimeStepSize, p->getMinTimeStepSize());
    nextMinTimeStepSize =
        std::min(nextMinTimeStepSize, p->getNextMinTimeStepSize());
  }

  if (printInfo) {
    logInfo("startNewTimeStep(...)",
            "step " << n << "\t t_min          =" << currentMinTimeStamp);

    logInfo("startNewTimeStep(...)",
            "\t\t dt_min         =" << currentMinTimeStepSize);

    logDebug("startNewTimeStep(...)",
            "\t\t next dt_min    =" << nextMinTimeStepSize); // only interesting for ADER-DG. Prints MAX_DOUBLE for finite volumes.
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
bool exahype::runners::Runner::setStableTimeStepSizesIfStabilityConditionWasHarmed(double factor) {
  bool cflConditionWasViolated = false;

  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    if (p->getType()==exahype::solvers::Solver::Type::ADER_DG) {
        bool solverTimeStepSizeIsInstable =
          static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinPredictorTimeStepSize()
          >
          static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinNextPredictorTimeStepSize();

        if (solverTimeStepSizeIsInstable) {
          static_cast<exahype::solvers::ADERDGSolver*>(p)->updateMinNextPredictorTimeStepSize(
              factor * static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinNextPredictorTimeStepSize());
          static_cast<exahype::solvers::ADERDGSolver*>(p)->setMinPredictorTimeStepSize(
              factor * static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinPredictorTimeStepSize());
        } else {
          static_cast<exahype::solvers::ADERDGSolver*>(p)->updateMinNextPredictorTimeStepSize(
              .5 * (static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinPredictorTimeStepSize() +
                  static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinNextPredictorTimeStepSize()));
        }

        cflConditionWasViolated =
            cflConditionWasViolated | solverTimeStepSizeIsInstable;
    }
  }

  return cflConditionWasViolated;  // | tooDiffusive;
}

void exahype::runners::Runner::recomputePredictorIfNecessary(
    exahype::repositories::Repository& repository,double factor) {
  // Must be evaluated before we start a new time step
  bool stabilityConditionWasHarmed = setStableTimeStepSizesIfStabilityConditionWasHarmed(factor);
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
  repository.switchToRiemannSolver();  // Riemann -> face2face
  repository.iterate();

  if (plot) {
    repository.switchToCorrectorAndPlot();  // Face to cell + Inside cell
  } else {
    repository.switchToCorrector();  // Face to cell + Inside cell
  }
  repository.iterate();

  repository.switchToPredictor();  // Cell onto faces
  repository.iterate();
}
