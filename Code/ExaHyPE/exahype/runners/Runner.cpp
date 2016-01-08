#include "exahype/runners/Runner.h"


#include "exahype/repositories/Repository.h"
#include "exahype/repositories/RepositoryFactory.h"

#include "peano/utils/UserInterface.h"

#include "tarch/Assertions.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"


// @todo Remove this include as soon as you've created your real-world geometry
#include "peano/geometry/Hexahedron.h" 


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

  // @todo Insert your code here
  
  // Start of dummy implementation
  
  repository.switchToInitialGrid(); repository.iterate();
  repository.switchToGridExport(); repository.iterate();
  repository.switchToPatchInitialisation(); repository.iterate();
  repository.switchToPatchInitialisationAndExport(); repository.iterate();
  repository.switchToFaceDataExchange(); repository.iterate();
  repository.switchToInitialConditionAndGlobalTimeStepComputation(); repository.iterate();
  repository.switchToInitialConditionAndExportAndGlobalTimeStepComputation(); repository.iterate();
  repository.switchToPredictorAndGlobalTimeStepComputation(); repository.iterate();
  repository.switchToCorrectorAndPredictorAndGlobalTimeStepComputation(); repository.iterate();
  repository.switchToCorrectorAndPredictorAndGlobalTimeStepComputationAndExport(); repository.iterate();

 
 
  repository.logIterationStatistics();
  repository.terminate();
  // End of dummy implementation

  return 0;
}
