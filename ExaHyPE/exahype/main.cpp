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

#include <vector>
#include <string>
#include <cstdlib> // getenv

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
  if (argc < 2) {
    logError("main()", "Usage: ./ExaHyPE config-file [additional args passed to Solver...]");
    return -1;
  }

  exahype::Parser parser;
  parser.readFile(argv[1]);

  if (!parser.isValid()) {
    logError("main()", "invalid config file. Quit");
    return -2;
  }
  
  // Collect all command line arguments for the Solvers
  std::vector<std::string> cmdlineargs(argv + 1, argv + argc);

  //
  //   Init solver registries
  // =====================================
  //
  kernels::initSolvers(parser, cmdlineargs);

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
  #elif defined(Asserts) || defined(Debug)
  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      false,  // logMachineName
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
      false,   // logTrace
      "");
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
if(! std::getenv("EXAHYPE_SKIP_TESTS")) { // cf issue #74
  tarch::tests::TestCaseRegistry::getInstance().getTestCaseCollection().run();
  int testExitCode = tarch::tests::TestCaseRegistry::getInstance()
                         .getTestCaseCollection()
                         .getNumberOfErrors();

  if (testExitCode != 0) {
    logError("main()", "unit tests failed. Quit.");
    return -2;
  }
} else {
  logInfo("main()", "Skipping tests as EXAHYPE_SKIP_TESTS is set."
     "We do so because tests are broken in the moment and nobody repairs them.");
} // end if getenv(EXAHYPE_SKIP_TESTS)
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

  kernels::finalise();

  return programExitCode;
}
