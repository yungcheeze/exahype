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


tarch::logging::Log  exahype::runners::Runner::_log( "exahype::runners::Runner" );



exahype::runners::Runner::Runner(const Parser& parser):
  _parser(parser) {
}


exahype::runners::Runner::~Runner() {
}


void exahype::runners::Runner::setupComputationalDomain() {
  peano::geometry::Hexahedron geometry(
    tarch::la::Vector<DIMENSIONS,double>( _parser.getSize() ),
    _parser.getOffset()
   );
}


void exahype::runners::Runner::initSharedMemoryConfiguration() {
  #ifdef SharedMemoryParallelisation
  const int numberOfThreads = parser.getNumberOfThreads();
  tarch::multicore::Core::getInstance().configure(numberOfThreads);



  ifstream f(name.c_str());
  bool multicorePropertiesFileDoesExist = f.good();
      f.close();


  if (_parser.useMulticoreAutotuning()) {
//      if (_parser.getMulticorePropertiesFile())
      // if no autotuning and profile, dann


  }
  else {
    #ifdef PerformanceAnalysis
      grain size sampling oracle
    #else

      dummy oracle
    #endif
  }
  #endif
}


void exahype::runners::Runner::shutdownSharedMemoryConfiguration() {
  #ifdef SharedMemoryParallelisation
  if (_parser.useMulticoreAutotuning()) {

      // zurueckschreiben der Daten in File
//      if (_parser.getMulticorePropertiesFile())
      // if no autotuning and profile, dann


  }
//#ifdef PerformanceAnalysis
  grain size sampling oracle
#endif
}



int exahype::runners::Runner::run() {
  setupComputationalDomain();

  initSharedMemoryConfiguration();

  exahype::repositories::Repository* repository = 
    exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
      geometry,
      tarch::la::Vector<DIMENSIONS,double>( _parser.getSize() ),
      _parser.getOffset()
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
