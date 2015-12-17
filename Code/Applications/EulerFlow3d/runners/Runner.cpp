#include "EulerFlow3d/runners/Runner.h"


#include "EulerFlow3d/repositories/Repository.h"
#include "EulerFlow3d/repositories/RepositoryFactory.h"

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
#include "EulerFlow3d/Constants.h"
// ! End of code for DG method

tarch::logging::Log  exahype::runners::Runner::_log( "exahype::runners::Runner" );

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
      tarch::la::Vector<DIMENSIONS,double>(1.0),
      tarch::la::Vector<DIMENSIONS,double>(0.0)
  );
  exahype::repositories::Repository* repository = 
      exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
          geometry,
          tarch::la::Vector<DIMENSIONS,double>(1.0),   // domainSize,
          tarch::la::Vector<DIMENSIONS,double>(0.0)    // computationalDomainOffset
      );
  // End of dummy implementation

  int result = 0;
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    result = runAsMaster( *repository );
  }
#ifdef Parallel
  else {
    result = runAsWorker( *repository );
  }
#endif

  delete repository;

  return result;
}

int exahype::runners::Runner::runAsMaster(exahype::repositories::Repository& repository) {
  peano::utils::UserInterface userInterface;
  userInterface.writeHeader();

#ifdef SharedMemoryParallelisation
  // We assume that the workload per cell is that big that we can set the enterCell
  // grain size to 1 as well as the minimum grain size. All the other values remain
  // the default values.

  peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
      new peano::datatraversal::autotuning::OracleForOnePhaseDummy(true));
#endif
  // The space-tree is initialised with 1 coarse grid cell on level 1 and 3^d fine grid cells on level 2.
  for (int coarseGridLevel=1; coarseGridLevel<EXAHYPE_INITIAL_GLOBAL_REFINEMENT_LEVEL; coarseGridLevel++) {
    repository.switchToInitialGrid();
    do {
      repository.iterate();
    } while (!repository.getState().isGridBalanced());
  }
  repository.switchToGridExport();                // export the grid
  repository.iterate();

  repository.switchToPatchInitialisation();       // initialize the cell descriptions;
  repository.iterate();

  repository.switchToInitialConditionAndExport(); // initialize the fields of the cell descriptions, i.e., the initial values.
  repository.iterate();

  repository.getState().setTimeStepSize(1e20);
  repository.switchToGlobalTimeStepComputation();
  repository.iterate();
  logInfo("runAsMaster(...)", "global time step (seconds)=" <<
          repository.getState().getTimeStepSize() );

  for (int n=0; n<100; n++) {
    // predictor
    repository.switchToPredictor();
    repository.iterate();

    // face data exchange
    repository.switchToFaceDataExchange();
    repository.iterate();

    // corrector
    if (n%EXAHYPE_PLOTTING_STRIDE==0) {
      repository.switchToCorrectorAndExport();
    } else {
      repository.switchToCorrector();
    }
    repository.iterate();

    // global reduction
    repository.getState().setTimeStepSize(1e20);
    repository.switchToGlobalTimeStepComputation();
    repository.iterate();
    logInfo("runAsMaster(...)", "[ExaHyPE] global time step="<< n <<", global time step size (seconds)=" <<
            repository.getState().getTimeStepSize() );
  }

  repository.logIterationStatistics();
  repository.terminate();
  // ! End of code for DG method

  return 0;
}
