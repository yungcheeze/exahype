#include "AdvectionDG0/runners/Runner.h"


#include "AdvectionDG0/repositories/Repository.h"
#include "AdvectionDG0/repositories/RepositoryFactory.h"

#include "peano/utils/UserInterface.h"

#include "tarch/Assertions.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"

// @todo Remove this include as soon as you've created your real-world geometry
#include "peano/geometry/Hexahedron.h" 

#include "AdvectionDG0/Constants.h"

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
  
  // ! Start of code for DG methdo
  const double CFL             = 1.;                                    // CFL number for a 0-th order DG method
  const double minimumMeshSize = 1./pow(3.,INITIAL_GLOBAL_REFINEMENT);  // tripartioning this is the smallest mesh size after initial refinement
  const double maximumVelocity = 1.;                                    // flux/CFL condition specific value (~max. eigenvalue). (MPI_REDUCTION necessary for global time scale)

  const double timeStepSize    = CFL * maximumVelocity/minimumMeshSize; // maximum flux value

  repository.getState().setTimeStepSize(timeStepSize);                  // application-specific value

  repository.switchToCreateGrid();

  for (int i=0; i<0)

  repository.logIterationStatistics();
  repository.terminate();
  // ! End of code for DG method

  return 0;
}
