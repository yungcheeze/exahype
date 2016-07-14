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
#include "exahype/Parser.h"

#include "tarch/Assertions.h"

#include <fstream>

#include <stdio.h>
#include <string.h>

#include "tarch/la/ScalarOperations.h"

tarch::logging::Log exahype::Parser::_log("exahype::Parser");


double exahype::Parser::getValueFromPropertyString( const std::string& parameterString, const std::string& key ) {
  std::size_t startIndex      = parameterString.find(key);
  startIndex      = parameterString.find(":",startIndex);
  std::size_t endIndexBracket = parameterString.find("}",startIndex+1);
  std::size_t endIndexComma   = parameterString.find(",",startIndex+1);

  std::size_t endIndex        = endIndexBracket<endIndexComma ? endIndexBracket : endIndexComma;

  std::string substring       = parameterString.substr(startIndex+1,endIndex-startIndex-1);

  double result;
  std::istringstream ss(substring);
  ss >> result;

  if (ss) {
    return result;
  }
  else {
    return std::numeric_limits<double>::quiet_NaN();
  }
}


exahype::Parser::Parser() {
  _identifier2Type.insert (
      std::pair<std::string,exahype::solvers::Solver::Type> ("ADER-DG", exahype::solvers::Solver::Type::ADER_DG) );

  _identifier2TimeStepping.insert (
      std::pair<std::string,exahype::solvers::Solver::TimeStepping> ("global", exahype::solvers::Solver::TimeStepping::Global) );
}

void exahype::Parser::readFile(const std::string& filename) {
  const int MAX_CHARS_PER_LINE = 512;
  const char* const DELIMITER = " =\t";

  _tokenStream.clear();

  std::ifstream inputFile;
  inputFile.open(filename.c_str());
  if (!inputFile.good()) {
    logError("readFile(String)", "can not open file " << filename);
    _tokenStream.clear();
    return;
  }

  bool currentlyReadsComment = false;

  while (!inputFile.eof()) {
    char lineBuffer[MAX_CHARS_PER_LINE];
    inputFile.getline(lineBuffer, MAX_CHARS_PER_LINE);

    char* token;

    // parse the line
    token = strtok(lineBuffer, DELIMITER);
    int currentTokenInLine = 0;
    while (token) {
      std::string newToken = token;
      if (currentlyReadsComment && newToken.find("*/") != std::string::npos) {
        currentlyReadsComment = false;
      } else if (!currentlyReadsComment &&
          newToken.find("/*") != std::string::npos) {
        currentlyReadsComment = true;
      } else if (!currentlyReadsComment) {
        logDebug("readFile(String)", "got token " << newToken);
        _tokenStream.push_back(newToken);
      }
      token = strtok(0, DELIMITER);  // subsequent tokens
      currentTokenInLine++;
    }
  }
}


bool exahype::Parser::isValid() const { return !_tokenStream.empty(); }


std::string exahype::Parser::getTokenAfter(std::string token,
                                           int additionalTokensToSkip) const {
  assertion(isValid());
  int currentToken = 0;
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
      _tokenStream[currentToken] != token) {
    currentToken++;
  }
  currentToken += (additionalTokensToSkip + 1);
  if (currentToken < static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  } else
    return "notoken";
}


std::string exahype::Parser::getTokenAfter(std::string token0,
                                           std::string token1,
                                           int additionalTokensToSkip) const {
  assertion(isValid());
  int currentToken = 0;
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
      _tokenStream[currentToken] != token0) {
    currentToken++;
  }
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
      _tokenStream[currentToken] != token1) {
    currentToken++;
  }
  currentToken += (additionalTokensToSkip + 1);
  if (currentToken < static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  } else
    return "notoken";
}

std::string exahype::Parser::getTokenAfter(std::string token0, int occurance0,
                                           int additionalTokensToSkip) const {
  assertion(isValid());
  assertion(occurance0 > 0);
  int currentToken = 0;
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
      (_tokenStream[currentToken] != token0 || occurance0 > 1)) {
    if (_tokenStream[currentToken] == token0) occurance0--;
    currentToken++;
  }
  currentToken += (additionalTokensToSkip + 1);
  if (currentToken < static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  } else
    return "notoken";
}

std::string exahype::Parser::getTokenAfter(std::string token0, int occurance0,
                                           std::string token1, int occurance1,
                                           int additionalTokensToSkip) const {
  assertion(isValid());
  assertion(occurance0 > 0);
  assertion(occurance1 > 0);
  int currentToken = 0;
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
      (_tokenStream[currentToken] != token0 || occurance0 > 1)) {
    if (_tokenStream[currentToken] == token0) occurance0--;
    currentToken++;
  }
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
      (_tokenStream[currentToken] != token1 || occurance1 > 1)) {
    if (_tokenStream[currentToken] == token1) occurance1--;
    currentToken++;
  }
  currentToken += (additionalTokensToSkip + 1);
  if (currentToken < static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  } else
    return "notoken";
}

int exahype::Parser::getNumberOfThreads() const {
  assertion(isValid());
  std::string token = getTokenAfter("shared-memory", "cores");
  logDebug("getNumberOfThreads()", "found token " << token);
  int result = atoi(token.c_str());
  if (result == 0) {
    logError("getNumberOfThreads()",
             "Invalid number of cores set: "
             << token << ". Use one core, i.e. switch off multithreading");
    result = 1;
  }
  return result;
}

tarch::la::Vector<DIMENSIONS, double> exahype::Parser::getDomainSize() const {
  assertion(isValid());
  std::string token;
  tarch::la::Vector<DIMENSIONS, double> result;
  token = getTokenAfter("computational-domain", "width", 0);
  result(0) = atof(token.c_str());
  token = getTokenAfter("computational-domain", "width", 1);
  result(1) = atof(token.c_str());
#if DIMENSIONS == 3
  token = getTokenAfter("computational-domain", "width", 2);
  result(2) = atof(token.c_str());
#endif
  return result;
}


tarch::la::Vector<DIMENSIONS, double> exahype::Parser::getBoundingBoxSize() const {
  tarch::la::Vector<DIMENSIONS, double> domainSize = getDomainSize();
  double longestH = tarch::la::max( domainSize );
  return tarch::la::Vector<DIMENSIONS, double>(longestH);
}


tarch::la::Vector<DIMENSIONS, double> exahype::Parser::getOffset() const {
  assertion(isValid());
  std::string token;
  tarch::la::Vector<DIMENSIONS, double> result;
  token = getTokenAfter("computational-domain", "offset", 0);
  result(0) = atof(token.c_str());
  token = getTokenAfter("computational-domain", "offset", 1);
  result(1) = atof(token.c_str());
#if DIMENSIONS == 3
  token = getTokenAfter("computational-domain", "offset", 2);
  result(2) = atof(token.c_str());
#endif
  logDebug("getSize()", "found offset " << result);
  return result;
}

std::string exahype::Parser::getMulticorePropertiesFile() const {
  std::string result = getTokenAfter("shared-memory", "properties-file");
  logDebug("getMulticorePropertiesFile()", "found token " << result);
  return result;
}


exahype::Parser::MPILoadBalancingType exahype::Parser::getMPILoadBalancingType() const {
  std::string token = getTokenAfter("distributed-memory", "identifier");
  exahype::Parser::MPILoadBalancingType result = MPILoadBalancingType::Static;
  if (token.compare("static_load_balancing") == 0) {
    result = MPILoadBalancingType::Static;
  }
  else {
    logError("getMPILoadBalancingType()", "Invalid distributed memory identifier " << token );
  }
  return result;
}


std::string exahype::Parser::getMPIConfiguration() const {
  return getTokenAfter("distributed-memory", "configure");
}


int exahype::Parser::getMPIBufferSize() const {
  std::string token = getTokenAfter("distributed-memory", "buffer-size");
  int result = atoi(token.c_str());
  if (result<=0) {
    logError("getMPIBufferSize()", "Invalid MPI buffer size " << token << ". reset to 64");
    result = 64;
  }
  return result;
}


int exahype::Parser::getMPITimeOut() const {
  std::string token = getTokenAfter("distributed-memory", "timeout");
  int result = atoi(token.c_str());
  if (result<=0) {
    logError("getMPIBufferSize()", "Invalid MPI timeout value " << token << ". reset to 0, i.e. no timeout");
    result = 0;
  }
  return result;
}


exahype::Parser::MulticoreOracleType exahype::Parser::getMulticoreOracleType()
const {
  std::string token = getTokenAfter("shared-memory", "identifier");
  exahype::Parser::MulticoreOracleType result = MulticoreOracleType::Dummy;
  if (token.compare("dummy") == 0) {
    result = MulticoreOracleType::Dummy;
  } else if (token.compare("autotuning") == 0) {
    result = MulticoreOracleType::Autotuning;
  } else if (token.compare("sampling") == 0) {
    result = MulticoreOracleType::GrainSizeSampling;
  } else {
    logError("getMulticoreOracleType()",
             "Invalid shared memory identifier "
             << token << ". Use dummy, autotuning, sampling. Set to dummy");
    result = MulticoreOracleType::Dummy;
  }
  return result;
}

double exahype::Parser::getSimulationEndTime() const {
  assertion(isValid());
  std::string token = getTokenAfter("computational-domain", "end-time");
  logDebug("getSimulationEndTime()", "found token " << token);
  double result = atof(token.c_str());
  if (result <= 0) {
    logError("getSimulationEndTime()",
             "Invalid simulation end-time: " << token << ". Use 1.0");
    result = 1.0;
  }
  return result;
}

bool exahype::Parser::getFuseAlgorithmicSteps() const {
  assertion(isValid());
  std::string token = getTokenAfter("optimisation", "fuse-algorithmic-steps");
  logDebug("getFuseAlgorithmicSteps()", "found fuse-algorithmic-steps" << token);
  return token.compare("on") == 0;
}

double exahype::Parser::getFuseAlgorithmicStepsFactor() const {
  assertion(isValid());
  std::string token = getTokenAfter("optimisation", "fuse-algorithmic-steps-factor");
  double result = atof(token.c_str());
  logDebug("getFuseAlgorithmicStepsFactor()", "found fuse-algorithmic-steps-factor " << token);

  if (result < 0.0 || result > 1.0) {
    logError("getFuseAlgorithmicStepsFactor()","'fuse-algorithmic-steps-factor': Value must be greater than zero and smaller than one.");
    std::cerr.flush();
    assert(false);
    exit(ASSERTION_EXIT_CODE);
  }

  return result;
}

exahype::solvers::Solver::Type exahype::Parser::getType(int solverNumber) const {
  assertion(isValid());
  std::string token;
  exahype::solvers::Solver::Type result;
  token = getTokenAfter("solver", solverNumber * 2 + 1, 0);
  if (_identifier2Type.find(token)!=_identifier2Type.end()) {
    result = _identifier2Type.at(token);
    logDebug("getType()", "found type " << result);
    return result;
  } else {
    logError("getType()","'" << getIdentifier(solverNumber) << "': 'type': Value '" << token << "' is invalid. See the ExaHyPE documentation for valid values.");
    std::cerr.flush();
    assert(false);
    exit(ASSERTION_EXIT_CODE);
  }
}

std::string exahype::Parser::getIdentifier(int solverNumber) const {
  assertion(isValid());
  std::string token;
  token = getTokenAfter("solver", solverNumber * 2 + 1, 1);
  logDebug("getIdentifier()", "found identifier " << token);
  return token;
}

int exahype::Parser::getVariables(int solverNumber) const {
  assertion(isValid());
  std::string token;
  int result;
  token = getTokenAfter("solver", solverNumber * 2 + 1, "variables", 1);
  result = atoi(token.c_str());

  if (result < 1) {
    logError("getVariables()","'" << getIdentifier(solverNumber) << "': 'variables': Value must be greater than zero.");
    std::cerr.flush();
    assert(false);
    exit(ASSERTION_EXIT_CODE);
  }

  logDebug("getVariables()", "found variables " << result);
  return result;
}

int exahype::Parser::getParameters(int solverNumber) const {
  assertion(isValid());
  std::string token;
  int result;
  token = getTokenAfter("solver", solverNumber * 2 + 1, "parameters", 1);
  result = atoi(token.c_str());

  if (result < 0) {
    logError("getParameters()","'" << getIdentifier(solverNumber) << "': 'parameters': Value must not be negative.");
    std::cerr.flush();
    assert(false);
    exit(ASSERTION_EXIT_CODE);
  }

  logDebug("getParameters()", "found parameters " << result);
  return result;
}

int exahype::Parser::getOrder(int solverNumber) const {
  assertion(isValid());
  std::string token;
  int result;
  token = getTokenAfter("solver", solverNumber * 2 + 1, "order", 1);
  result = atof(token.c_str());

  if (result < 0) {
    logError("getOrder()","'" << getIdentifier(solverNumber) << "': 'order': Value must not be negative.");
    std::cerr.flush();
    assert(false);
    exit(ASSERTION_EXIT_CODE);
  }

  logDebug("getOrder()", "found order " << result);
  return result;
}

double exahype::Parser::getMaximumMeshSize(int solverNumber) const {
  assertion(isValid());
  std::string token;
  double result;
  token = getTokenAfter("solver", solverNumber * 2 + 1, "maximum-mesh-size", 1, 0);
  result = atof(token.c_str());
  if (tarch::la::smallerEquals(result,0.0)) {
    logError("getMaximumMeshSize()","'" << getIdentifier(solverNumber) << "': 'maximum-mesh-size': Value must be greater than zero.");
    std::cerr.flush();
    assert(false);
    exit(ASSERTION_EXIT_CODE);
  }

  logDebug("getMaximumMeshSize()", "found maximum mesh size " << result);
  return result;
}

exahype::solvers::Solver::TimeStepping exahype::Parser::getTimeStepping(int solverNumber) const {
  assertion(isValid());
  std::string token;
  exahype::solvers::Solver::TimeStepping result;
  token = getTokenAfter("solver", solverNumber * 2 + 1, "time-stepping", 1);
  if (_identifier2TimeStepping.find(token)!=_identifier2TimeStepping.end()) {
    result = _identifier2TimeStepping.at(token);
    logDebug("getTimeStepping()", "found TimeStepping " << result);
    return result;
  } else {
    logError("getTimeStepping()","'" << getIdentifier(solverNumber) << "': 'time-stepping': Value '" << token << "' is invalid. See the ExaHyPE documentation for valid values.");
    std::cerr.flush();
    assert(false);
    exit(ASSERTION_EXIT_CODE);
  }
  return exahype::solvers::Solver::TimeStepping::Global;
}

double exahype::Parser::getFirstSnapshotTimeForPlotter(
    int solverNumber, int plotterNumber) const {
  assertion(isValid());
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber * 2 + 1, "plot",
                                    plotterNumber * 2 + 1, 2);
  logDebug("getFirstSnapshotTimeForPlotter()", "found token " << token);
  assertion3(token.compare("notoken") != 0, token, solverNumber, plotterNumber);
  return atof(token.c_str());
}

double exahype::Parser::getRepeatTimeForPlotter(int solverNumber,
                                                int plotterNumber) const {
  assertion(isValid());
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber * 2 + 1, "plot",
                                    plotterNumber * 2 + 1, 4);
  logDebug("getRepeatTimeForPlotter()", "found token " << token);
  assertion3(token.compare("notoken") != 0, token, solverNumber, plotterNumber);
  return atof(token.c_str());
}

std::string exahype::Parser::getIdentifierForPlotter(int solverNumber,
                                                     int plotterNumber) const {
  assertion(isValid());
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber * 2 + 1, "plot",
                                    plotterNumber * 2 + 1);
  logDebug("getIdentifierForPlotter()", "found token " << token);
  assertion3(token.compare("notoken") != 0, token, solverNumber, plotterNumber);
  return token;
}

std::string exahype::Parser::getFilenameForPlotter(int solverNumber,
                                                   int plotterNumber) const {
  assertion(isValid());
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber * 2 + 1, "plot",
                                    plotterNumber * 2 + 1, 6);
  logDebug("getFilenameForPlotter()", "found token " << token);
  assertion3(token.compare("notoken") != 0, token, solverNumber, plotterNumber);
  return token;
}

std::string exahype::Parser::getSelectorForPlotter(int solverNumber, int plotterNumber) const {
  assertion(isValid());
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber * 2 + 1, "plot",
                                    plotterNumber * 2 + 1, 8);
  logDebug("getSelectorForPlotter()", "found token " << token);
  assertion3(token.compare("notoken") != 0, token, solverNumber, plotterNumber);
  return (token != "notoken") ? token : "{}";
}

std::string exahype::Parser::getProfilerIdentifier() const {
  assertion(isValid());
  std::string token = getTokenAfter("profiling", "profiler");
  logDebug("getProfilerIdentifier()", "found token" << token);
  return (token != "notoken") ? token : "NoOpProfiler";
}

std::string exahype::Parser::getMetricsIdentifierList() const {
  assertion(isValid());
  std::string token = getTokenAfter("profiling", "metrics");
  logDebug("getMetricsIdentifierList()", "found token " << token);
  return (token != "notoken") ? token : "{}";
}

void exahype::Parser::logSolverDetails(int solverNumber) const {
  logInfo("logSolverDetails()", "Solver "               << getTokenAfter("solver", solverNumber * 2 + 1, 0) << " " << getIdentifier(solverNumber) << ":");
  logInfo("logSolverDetails()", "variables:\t\t"        << getVariables(solverNumber));
  logInfo("logSolverDetails()", "parameters:\t"         << getParameters(solverNumber));
  logInfo("logSolverDetails()", "order:\t\t"            << getOrder(solverNumber));
  logInfo("logSolverDetails()", "maximum-mesh-size:\t"  << getMaximumMeshSize(solverNumber));
  logInfo("logSolverDetails()", "time-stepping:\t"      << getTokenAfter("solver", solverNumber * 2 + 1, "time-stepping", 1));
}


void exahype::Parser::checkSolverConsistency(int solverNumber) const {
  assertion1(solverNumber < static_cast<int>(exahype::solvers::RegisteredSolvers.size()),solverNumber);
  exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[solverNumber];

  bool recompile       = false;
  bool runToolkitAgain = false;
  if (solver->getType() != getType(solverNumber)) {
    logError("checkIfSolverIsConsistent","'" << getIdentifier(solverNumber) <<
             "': Solver type in specification file" <<
             "differs from implementation solver type.");
    recompile = true;
  }

  if (solver->getIdentifier().compare(getIdentifier(solverNumber))) {
    logError("checkSolverConsistency","'" << getIdentifier(solverNumber) <<
             "': Identifier in specification file " <<
             "('" << getIdentifier(solverNumber) << "') differs from identifier used in implementation ('" << solver->getIdentifier() << "').");
    recompile = true;
  }

  if (solver->getNumberOfVariables() != getVariables(solverNumber)) {
    logError("checkSolverConsistency","'" << getIdentifier(solverNumber) <<
             "': Value for 'variables' in specification file" <<
             "('" << getVariables(solverNumber) << "') differs from number of variables used in implementation file ('" << solver->getNumberOfVariables() << "').");
    recompile = true;
  }

  if (solver->getNumberOfParameters() != getParameters(solverNumber)) {
    logError("checkSolverConsistency","'" << getIdentifier(solverNumber) <<
             "': Value for field 'parameters' in specification file" <<
             "('" << getParameters(solverNumber) << "') differs from  number of parameters used in implementation file ('" << solver->getNumberOfParameters() << "').");
    recompile = true;
  }

  if (solver->getNodesPerCoordinateAxis() != getOrder(solverNumber)+1) {
    logError("checkSolverConsistency","'" << getIdentifier(solverNumber) <<
             "': Value for field 'order' in specification file " <<
             "('" << getOrder(solverNumber) << "') differs from value used in implementation file ('" << solver->getNodesPerCoordinateAxis()-1 << "'). ");
    runToolkitAgain = true;
  }

  if (runToolkitAgain) {
    logError("checkSolverConsistency","Please (1) run the Toolkit again, and (2) recompile!");
    std::cerr.flush();
    exit(1);
  }

  if (recompile) {
    logError("checkSolverConsistency","Please (1) adjust the specification file (*.exahype) or the file '" << solver->getIdentifier() << ".cpp' accordingly, and (2) recompile!");
    std::cerr.flush();
    exit(1);
  }
}
