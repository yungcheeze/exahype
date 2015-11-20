#include "ExaHyPE/runners/Runner.h"


#ifdef Parallel
#include "peano/utils/Globals.h"
#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/parallel/NodePool.h"
#include "peano/parallel/messages/ForkMessage.h"
#include "ExaHyPE/repositories/Repository.h"


int ExaHyPE::runners::Runner::runAsWorker(ExaHyPE::repositories::Repository& repository) {
  int newMasterNode = tarch::parallel::NodePool::getInstance().waitForJob(); 
  while ( newMasterNode != tarch::parallel::NodePool::JobRequestMessageAnswerValues::Terminate ) {
    if ( newMasterNode >= tarch::parallel::NodePool::JobRequestMessageAnswerValues::NewMaster ) {
      peano::parallel::messages::ForkMessage forkMessage;
      forkMessage.receive(tarch::parallel::NodePool::getInstance().getMasterRank(),tarch::parallel::NodePool::getInstance().getTagForForkMessages(), true, ReceiveIterationControlMessagesBlocking);

      repository.restart(
        forkMessage.getH(),
        forkMessage.getDomainOffset(),
        forkMessage.getLevel(),
        forkMessage.getPositionOfFineGridCellRelativeToCoarseGridCell()
      );
  
      bool continueToIterate = true;
      while (continueToIterate) {
        switch (repository.continueToIterate()) {
          case ExaHyPE::repositories::Repository::Continue:
            repository.iterate();  
            break;
          case ExaHyPE::repositories::Repository::Terminate:
            continueToIterate = false;  
            break;
          case ExaHyPE::repositories::Repository::RunGlobalStep:
            runGlobalStep();  
            break;
        }
      }
    
      // insert your postprocessing here
      // -------------------------------  

      // -------------------------------

      repository.terminate();
    }
    else if ( newMasterNode == tarch::parallel::NodePool::JobRequestMessageAnswerValues::RunAllNodes ) {
      runGlobalStep();  
    }
    newMasterNode = tarch::parallel::NodePool::getInstance().waitForJob(); 
  }
  return 0;
}


void ExaHyPE::runners::Runner::runGlobalStep() {
    // insert yourcode here
    // -------------------------------

    // -------------------------------
}
#endif
