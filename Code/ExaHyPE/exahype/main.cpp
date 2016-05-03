#include "tarch/logging/Log.h"
#include "tarch/tests/TestCaseRegistry.h"
#include "tarch/logging/CommandLineLogger.h"
#include "tarch/logging/LogFilterFileReader.h"
#include "tarch/parallel/Node.h"

#include "peano/peano.h"

#include "exahype/Parser.h"
#include "exahype/runners/Runner.h"

#include "kernels/KernelCalls.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

tarch::logging::Log _log("");

int main(int argc, char** argv) {
  peano::fillLookupTables();

  //
  //   Setup environment
  // =====================
  //
  int parallelSetup = peano::initParallelEnvironment(&argc, &argv);
  if (parallelSetup != 0) {
#ifdef Parallel
    // Please do not use the logging if MPI doesn't work properly.
    std::cerr << "mpi initialisation wasn't successful. Application shut down"
              << std::endl;
#else
    _log.error("main()",
               "mpi initialisation wasn't successful. Application shut down");
#endif
    return parallelSetup;
  }

  int sharedMemorySetup = peano::initSharedMemoryEnvironment();
  if (sharedMemorySetup != 0) {
    logError("main()",
             "shared memory initialisation wasn't successful. Application shut "
             "down");
    return sharedMemorySetup;
  }

  //
  //   Parse config file
  // =====================
  //
  if (argc != 2) {
    logError("main()", "Usage: ./ExaHyPE config-file");
    return -1;
  }

  exahype::Parser parser;
  parser.readFile(argv[1]);

  if (!parser.isValid()) {
    logError("main()", "invalid config file. Quit");
    return -2;
  }

  // @todo 04/02/16:Dominic Etienne Charrier
  // initGauss.. and initDG.. must be called before
  // we can run any kernel tests.
  //
  //   Init registries and lookup tables
  // =====================================
  //
  kernels::initSolvers(parser);
  kernels::initGaussLegendreNodesAndWeights();
  kernels::initDGMatrices();

  //
  //   Configure the logging
  // =========================
  //
  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
#ifdef Parallel
  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      true,   // logMachineName
      true,   // logMessageType
      true,   // logTrace
      "exahype.log-file");
#else
  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      false,  // logMachineName
      true,   // logMessageType
      true,   // logTrace
      "exahype.log-file");
#endif

  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  if (!tarch::logging::LogFilterFileReader::parsePlainTextFile(
          "exahype.log-filter")) {
    tarch::logging::CommandLineLogger::getInstance().clearFilterList();
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry("info", false));
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry(
            "info", -1, "peano::grid", true));
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry("debug", true));
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry("debug", -1,
                                                             "exahype", false));
  }

//
//   Run tests
// =============
//
#if defined(Debug) || defined(Asserts)
  tarch::tests::TestCaseRegistry::getInstance().getTestCaseCollection().run();
  int testExitCode = tarch::tests::TestCaseRegistry::getInstance()
                         .getTestCaseCollection()
                         .getNumberOfErrors();

  if (testExitCode != 0) {
    logError("main()", "unit tests failed. Quit");
    return -2;
  }
#endif

  exahype::runners::Runner runner(parser);
  int programExitCode = runner.run();

  if (programExitCode == 0) {
#ifdef Parallel
    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      logInfo("main()", "Peano terminates successfully");
    }
#else
    logInfo("main()", "Peano terminates successfully");
#endif
  } else {
    logInfo("main()", "quit with error code " << programExitCode);
  }

  peano::shutdownParallelEnvironment();
  peano::shutdownSharedMemoryEnvironment();
  peano::releaseCachedData();

  kernels::freeDGMatrices();

  return programExitCode;
}
