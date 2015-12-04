#include "EulerFlow3d/runners/Runner.h"


#include "EulerFlow3d/repositories/Repository.h"
#include "EulerFlow3d/repositories/RepositoryFactory.h"

#include "peano/utils/UserInterface.h"

#include "tarch/Assertions.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"

// @todo Remove this include as soon as you've created your real-world geometry
#include "peano/geometry/Hexahedron.h" 

// ! Begin of code for DG method
#include "EulerFlow3d/Constants.h"
// ! End of code for DG method

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

  // ! Begin of code for DG method
  const double CFL             = 0.1;                                                     // CFL number for a 0-th order DG method
  const double minimumMeshSize = 1./std::pow(3.,EXAHYPE_INITIAL_GLOBAL_REFINEMENT_LEVEL);  // tripartioning; this is the smallest fine grid mesh size after initial refinement
  const double maximumVelocity = 1.;                                                       // flux/CFL condition specific value (~max. eigenvalue). (MPI_REDUCTION necessary in order to obtain global time scale)

  const double timeStepSize    = CFL * minimumMeshSize/maximumVelocity;
  const int    maximumTimeStep = std::ceil(EXAHYPE_SIMULATION_TIME/timeStepSize);
  const int    plottingStride  = std::floor(maximumTimeStep/EXAHYPE_NUMBER_OF_PLOTS);

  repository.getState().setTimeStepSize(1e-3);

  // The space-tree is initialised with 1 coarse grid cell on level 1 and 3^d fine grid cells on level 2.
  for (int coarseGridLevel=1; coarseGridLevel<EXAHYPE_INITIAL_GLOBAL_REFINEMENT_LEVEL; coarseGridLevel++) {
    repository.switchToInitialGrid();
    do {
      repository.iterate();
    } while (!repository.getState().isGridStationary());
  }

  repository.switchToGridExport();          // export the grid
  repository.iterate();

  repository.switchToPatchInit();           // initialize the cell descriptions;
  repository.iterate();
  repository.switchToInitialCondition();    // initialize the fields of the cell descriptions, i.e., the initial values.
  repository.iterate();

  for (int n=0; n<240; n++) {
    repository.switchToPredictor();
    repository.iterate();
    repository.switchToCorrector();
    repository.iterate();
    repository.switchToSolutionExport();
    repository.iterate();
//    if (plottingStride>0 && n%plottingStride==0) {
//
//    } else {
//      repository.switchToTimeStep();
//    }
    repository.iterate();
  }

  repository.logIterationStatistics();
  repository.terminate();
  // ! End of code for DG method

  return 0;
}
