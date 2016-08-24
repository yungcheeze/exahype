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

#include "../../../Peano/mpibalancing/HotspotBalancing.h"
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

#ifdef Parallel
#include "mpibalancing/GreedyBalancing.h"
#include "mpibalancing/FairNodePoolStrategy.h"
#endif
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
      else if (configuration.find( "fair" )!=std::string::npos ) {
        int ranksPerNode = static_cast<int>(exahype::Parser::getValueFromPropertyString(configuration,"ranks_per_node"));
        if (ranksPerNode<=0) {
          logError( "initDistributedMemoryConfiguration()", "please inform fair balancing how many ranks per node you use through value \"ranks_per_node:XXX\". Read value " << ranksPerNode << " is invalid" );
          ranksPerNode = 1;
        }
        if ( ranksPerNode>=tarch::parallel::Node::getInstance().getNumberOfNodes() ) {
          logWarning( "initDistributedMemoryConfiguration()", "value \"ranks_per_node:XXX\" exceeds total rank count. Reset to 1" );
          ranksPerNode = 1;
        }
        tarch::parallel::NodePool::getInstance().setStrategy(
            new mpibalancing::FairNodePoolStrategy(ranksPerNode)
        );
        logInfo("initDistributedMemoryConfiguration()", "load balancing relies on fair answering strategy with " << ranksPerNode << " rank(s) per node") ;
      }
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
          new mpibalancing::GreedyBalancing(getCoarsestGridLevelOfAllSolvers(),getCoarsestGridLevelOfAllSolvers())
      );
    }
    else if ( configuration.find( "hotspot" )!=std::string::npos ) {
      logInfo("initDistributedMemoryConfiguration()", "use global hotspot elimination without joins (mpibalancing/StaticBalancing)");
      peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
          new mpibalancing::HotspotBalancing(false,getCoarsestGridLevelOfAllSolvers())
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


int exahype::runners::Runner::getCoarsestGridLevelOfAllSolvers() const {
  double boundingBox = _parser.getBoundingBoxSize()(0);
  double hMax        = exahype::solvers::Solver::getCoarsestMeshSizeOfAllSolvers();

  int    result      = 1;
  double currenthMax = std::numeric_limits<double>::max();
  while (currenthMax>hMax) {
    currenthMax = boundingBox / threePowI(result);
    result++;
  }

  logDebug( "getCoarsestGridLevelOfAllSolvers()", "regular grid depth of " << result << " creates grid with h_max=" << currenthMax );
  return std::max(3,result);
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
      ". bounding box has size " << _parser.getBoundingBoxSize() <<
      ". grid regular up to level " << getCoarsestGridLevelOfAllSolvers() << " (level 1 is coarsest available cell in tree)" );
#ifdef Parallel
  const double boundingBoxScaling = static_cast<double>(getCoarsestGridLevelOfAllSolvers()) / (static_cast<double>(getCoarsestGridLevelOfAllSolvers())-2);
  assertion4(boundingBoxScaling>=1.0, boundingBoxScaling, getCoarsestGridLevelOfAllSolvers(), _parser.getDomainSize(), _parser.getBoundingBoxSize() );
  const double boundingBoxShift   = (1.0-boundingBoxScaling)/2.0;
  assertion5(boundingBoxShift<=0.0, boundingBoxScaling, getCoarsestGridLevelOfAllSolvers(), _parser.getDomainSize(), _parser.getBoundingBoxSize(), boundingBoxScaling );

  logInfo(
      "run(...)",
      "increase domain artificially by " << boundingBoxScaling << " and shift bounding box by " << boundingBoxShift << " to simplify load balancing along boundary");
  return exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
      geometry,
      _parser.getBoundingBoxSize()*boundingBoxScaling,
      _parser.getOffset()+boundingBoxShift*_parser.getBoundingBoxSize()
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

  int gridSetupIterationsToRun = 3;
  while (gridSetupIterationsToRun>0) {
    repository.iterate();
    gridSetupIterations++;

    if ( UseStationaryCriterion && repository.getState().isGridStationary() ) {
      gridSetupIterationsToRun--;
    }
    else if ( !repository.getState().isGridBalanced() && tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()>0 ) {
      gridSetupIterationsToRun=3;  // we need at least 3 sweeps to recover from ongoing balancing
    }
    else if ( !repository.getState().isGridBalanced()  ) {
      gridSetupIterationsToRun=1;  // one additional step to get adjacency right
    }
    else {
      gridSetupIterationsToRun--;
    }

#if defined(TrackGridStatistics) && defined(Asserts)
    logInfo("runAsMaster()",
        "grid setup iteration #" << gridSetupIterations <<
        ", max-level=" << repository.getState().getMaxLevel() <<
        ", state=" << repository.getState().toString() <<
        ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
#elif defined(Asserts)
    logInfo("runAsMaster()",
        "grid setup iteration #" << gridSetupIterations <<
        ", state=" << repository.getState().toString() <<
        ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
#elif defined(TrackGridStatistics)
    logInfo("runAsMaster()",
        "grid setup iteration #" << gridSetupIterations <<
        ", max-level=" << repository.getState().getMaxLevel() <<
        ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
#else
    logInfo("runAsMaster()",
        "grid setup iteration #" << gridSetupIterations <<
        ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
#endif
  }

  logInfo("createGrid(Repository)", "finished grid setup after " << gridSetupIterations << " iterations" );

  if (tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()>0) {
    logWarning( "createGrid(Repository)", "there are still " << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes() << " ranks idle" )
  }

#ifdef Parallel
  // Might be too restrictive for later runs. Remove but keep warning from above
  assertion( tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()==0 );
#endif
}


int exahype::runners::Runner::runAsMaster(exahype::repositories::Repository& repository) {
  peano::utils::UserInterface::writeHeader();

  initSolverTimeStamps();
  createGrid(repository);
  /*
   * Set ADER-DG corrector time stamp and finite volumes time stamp.
   * Compute ADER-DG corrector time step size implicitly and finite volumes time step size.
   * (Implicitly means here that we set the predictor time step size but after the next newTimeStep(...)
   * the corrector time step size is set as this predictor time step size.)
   *
   * Note that it is important that we run SolutionAdjustmentAnd
   * GlobalTimeStepComputation directly after the grid setup
   * since we receive here the metadata
   * that was sent in the last iteration of the grid setup.
   */
  initSolverTimeStamps();
  repository.switchToSolutionAdjustmentAndGlobalTimeStepComputation();
  repository.iterate();

#if defined(Dim2) && defined(Asserts)
  repository.switchToPlotAugmentedAMRGrid();
  repository.iterate();
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
   * Reset the time stamps of the finite volumes solvers.
   *
   * !!! Rationale
   * Unlike for the rearranged ADER-DG scheme, we only
   * need one warm up iteration for the finite volumes solvers.
   * But since we have to perform two time step computations
   * to warm up the ADER-DG schemes,
   * the finite volumes solvers think they are already
   * advanced by one time step after the computation
   * of the second time step size.
   */
  initFiniteVolumesSolverTimeStamps();

  /*
   * Finally print the initial time step info.
   */
  printTimeStepInfo(-1);

  const double simulationEndTime = _parser.getSimulationEndTime();

  logDebug("runAsMaster(...)","min solver time stamp: "     << solvers::Solver::getMinSolverTimeStampOfAllSolvers());
  logDebug("runAsMaster(...)","min solver time step size: " << solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers());

  while ((solvers::Solver::getMinSolverTimeStampOfAllSolvers() < simulationEndTime) &&
      tarch::la::greater(solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers(), 0.0)) {
    bool plot = exahype::plotters::isAPlotterActive(
        solvers::Solver::getMinSolverTimeStampOfAllSolvers());

    if (_parser.getFuseAlgorithmicSteps()) {
      int numberOfStepsToRun = 1;
      if (plot) {
        numberOfStepsToRun = 0;
      }
      else if (solvers::Solver::allSolversUseTimeSteppingScheme(solvers::Solver::TimeStepping::GlobalFixed)) {
        /**
         * This computation is optimistic. If we were pessimistic, we had to
         * use the max solver time step size. However, this is not necessary
         * here, as we half the time steps anyway.
         */
    	if (solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers()>0.0) {
    	  const double timeIntervalTillNextPlot = std::min(exahype::plotters::getTimeOfNextPlot(),simulationEndTime) - solvers::Solver::getMaxSolverTimeStampOfAllSolvers();
          numberOfStepsToRun = std::floor( timeIntervalTillNextPlot / solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers() );
          logDebug("runAsMaster(...)", "number of possible time steps=" << numberOfStepsToRun << " with time till next plot=" << timeIntervalTillNextPlot );
    	}
        numberOfStepsToRun = numberOfStepsToRun<2 ? 1 : numberOfStepsToRun/2;
      }

      runOneTimeStampWithFusedAlgorithmicSteps(repository, numberOfStepsToRun);
      recomputePredictorIfNecessary(repository,_parser.getFuseAlgorithmicStepsFactor());
      printTimeStepInfo(numberOfStepsToRun);
    } else {
      runOneTimeStampWithThreeSeparateAlgorithmicSteps(repository, plot);
      printTimeStepInfo(1);
    }

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

void exahype::runners::Runner::printTimeStepInfo(int numberOfStepsRanSinceLastCall) {
  double currentMinTimeStamp    = std::numeric_limits<double>::max();
  double currentMinTimeStepSize = std::numeric_limits<double>::max();
  double nextMinTimeStepSize    = std::numeric_limits<double>::max();

  static int n = 0;
  if (numberOfStepsRanSinceLastCall==0) {
    n++;
  }
  else if (numberOfStepsRanSinceLastCall>0) {
    n+=numberOfStepsRanSinceLastCall;
  }

  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    currentMinTimeStamp =
        std::min(currentMinTimeStamp, p->getMinTimeStamp());
    currentMinTimeStepSize =
        std::min(currentMinTimeStepSize, p->getMinTimeStepSize());
    nextMinTimeStepSize =
        std::min(nextMinTimeStepSize, p->getNextMinTimeStepSize());
  }

  logInfo("startNewTimeStep(...)",
      "step " << n << "\t t_min          =" << currentMinTimeStamp);

  logInfo("startNewTimeStep(...)",
      "\t\t dt_min         =" << currentMinTimeStepSize);

  logDebug("startNewTimeStep(...)",
      "\t\t memoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB"); 

  logDebug("startNewTimeStep(...)",
      "\t\t next dt_min    =" << nextMinTimeStepSize); // Only interesting for ADER-DG. Prints MAX_DOUBLE for finite volumes.
#if defined(Debug) || defined(Asserts)
  tarch::logging::CommandLineLogger::getInstance().closeOutputStreamAndReopenNewOne();
#endif
}

void exahype::runners::Runner::runOneTimeStampWithFusedAlgorithmicSteps(
    exahype::repositories::Repository& repository, int numberOfStepsToRun) {
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

  if (numberOfStepsToRun==0) {
    repository.switchToADERDGTimeStepAndPlot();
    repository.iterate();
  } else {
    repository.switchToADERDGTimeStep();
    repository.iterate(numberOfStepsToRun);
  }
}

bool exahype::runners::Runner::setStableTimeStepSizesIfStabilityConditionWasHarmed(double factor) {
  assertion1(tarch::la::smallerEquals(factor,1.) && tarch::la::greater(factor,0.),"Factor must be smaller or equal to 1.");
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

      cflConditionWasViolated = cflConditionWasViolated | solverTimeStepSizeIsInstable;
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

void exahype::runners::Runner::runOneTimeStampWithThreeSeparateAlgorithmicSteps(
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
