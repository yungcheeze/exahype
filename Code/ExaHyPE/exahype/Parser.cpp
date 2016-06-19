#include "exahype/Parser.h"
#include "tarch/Assertions.h"

#include <fstream>

#include <stdio.h>
#include <string.h>

tarch::logging::Log exahype::Parser::_log("exahype::Parser");

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

double exahype::Parser::getSize() const {
  assertion(isValid());
  std::string token = getTokenAfter("computational-domain", "width");
  logDebug("getSize()", "found token " << token);
  double result = atof(token.c_str());
  if (result <= 0) {
    logError("getSize()", "Invalid width of computational domain: "
                              << token << ". Use unit cube");
    result = 1.0;
  }
  return result;
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

exahype::Parser::MulticoreOracleType exahype::Parser::getMulticoreOracleType()
    const {
  std::string token = getTokenAfter("shared-memory", "identifier");
  exahype::Parser::MulticoreOracleType result = Dummy;
  if (token.compare("dummy") == 0) {
    result = Dummy;
  } else if (token.compare("autotuning") == 0) {
    result = Autotuning;
  } else if (token.compare("sampling") == 0) {
    result = GrainSizeSampling;
  } else {
    logError("getMulticoreOracleType()",
             "Invalid shared memory identifier "
                 << token << ". Use dummy, autotuning, sampling. Set to dummy");
    result = Dummy;
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

bool exahype::Parser::fuseAlgorithmicSteps() const {
  assertion(isValid());
  std::string token = getTokenAfter("optimisation", "fuse-algorithmic-steps");
  logDebug("fuseAlgorithmicSteps()", "found token " << token);
  return token.compare("on") == 0;
}

double exahype::Parser::getFirstSnapshotTimeForPlotter(
    int solverNumber, int plotterNumber) const {
  assertion(isValid());
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber * 2 + 1, "plot",
                                    plotterNumber * 2 + 1, 2);
  logDebug("fuseAlgorithmicSteps()", "found token " << token);
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
  logDebug("fuseAlgorithmicSteps()", "found token " << token);
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
  logDebug("fuseAlgorithmicSteps()", "found token " << token);
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
  logDebug("fuseAlgorithmicSteps()", "found token " << token);
  assertion3(token.compare("notoken") != 0, token, solverNumber, plotterNumber);
  return token;
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
