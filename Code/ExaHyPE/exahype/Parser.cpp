#include "exahype/Parser.h"
#include "tarch/Assertions.h"

#include <fstream>

#include <stdio.h>
#include <string.h>



tarch::logging::Log exahype::Parser::_log( "exahype::Parser" );


void exahype::Parser::readFile( const std::string& filename ) {
  const int MAX_CHARS_PER_LINE = 512;
  const int MAX_TOKENS_PER_LINE = 20;
  const char* const DELIMITER = " =";

  _tokenStream.clear();

  std::ifstream inputFile;
  inputFile.open( filename.c_str() );
  if (!inputFile.good()) {
    logError( "readFile(String)", "can not open file " << filename );
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
      if (currentlyReadsComment && newToken.find("*/")!=std::string::npos) {
        currentlyReadsComment = false;
      }
      else if (!currentlyReadsComment && newToken.find("/*")!=std::string::npos) {
        currentlyReadsComment = true;
      }
      else if (!currentlyReadsComment) {
        logDebug( "readFile(String)", "got token " << newToken );
        _tokenStream.push_back( newToken );
      }
      token = strtok(0, DELIMITER); // subsequent tokens
      currentTokenInLine++;
    }
  }

}


bool exahype::Parser::isValid() const {
  return  !_tokenStream.empty();
}


std::string exahype::Parser::getTokenAfter( std::string token, int additionalTokensToSkip ) const {
  assertion( isValid() );
  int currentToken = 0;
  while (_tokenStream[currentToken]!=token && currentToken<static_cast<int>(_tokenStream.size())) {
    currentToken++;
  }
  currentToken += (additionalTokensToSkip+1);
  if ( currentToken<static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  }
  else return "notoken";
}


std::string exahype::Parser::getTokenAfter( std::string token0, std::string token1, int additionalTokensToSkip ) const {
  assertion( isValid() );
  int currentToken = 0;
  while (_tokenStream[currentToken]!=token0 && currentToken<static_cast<int>(_tokenStream.size())) {
    currentToken++;
  }
  while (_tokenStream[currentToken]!=token1 && currentToken<static_cast<int>(_tokenStream.size())) {
    currentToken++;
  }
  currentToken += (additionalTokensToSkip+1);
  if ( currentToken<static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  }
  else return "notoken";
}


int exahype::Parser::getNumberOfThreads() {
  assertion( isValid() );
  std::string token = getTokenAfter("shared-memory","cores");
  logDebug( "getNumberOfThreads()", "found token " << token );
  int result = atoi( token.c_str() );
  if (result==0) {
    logError( "getNumberOfThreads()", "Invalid number of cores set: " << token << ". Use one core, i.e. switch off multithreading" );
    result = 1;
  }
  return result;
}


double exahype::Parser::getSize() const {
  assertion( isValid() );
  std::string token = getTokenAfter("computational-domain","size");
  // @todo change into Debug
  logInfo( "getSize()", "found token " << token );
  double result = atof( token.c_str() );
  if (result<=0) {
    logError( "getSize()", "Invalid size of computational domain: " << token << ". Use unit cube" );
    result = 1.0;
  }
  return result;
}


tarch::la::Vector<DIMENSIONS,double> exahype::Parser::getOffset() const {
  assertion( isValid() );
  std::string token;
  tarch::la::Vector<DIMENSIONS,double> result;
  token     = getTokenAfter("computational-domain","offset", 0);
  result(0) = atof( token.c_str() );
  token     = getTokenAfter("computational-domain","offset", 1);
  result(1) = atof( token.c_str() );
  token     = getTokenAfter("computational-domain","offset", 2);
  result(2) = atof( token.c_str() );
  // @todo change into Debug
  logInfo( "getSize()", "found offset " << result );
  return result;
}



std::string exahype::Parser::getMulticorePropertiesFile() const {
  std::string result =  getTokenAfter("shared-memory","properties-file");
  // @todo Change into Debug
  logInfo( "getMulticorePropertiesFile()", "found token " << result );
  return result;
}


bool exahype::Parser::useMulticoreAutotuning() const {
  std::string token =  getTokenAfter("shared-memory","properties-file");
  // @todo Change into Debug
  logInfo( "useMulticoreAutotuning()", "found token " << token );
  if (token.compare("off")==0) {
    return false;
  }
  else if (token.compare("on")==0) {
    return true;
  }
  else {
    logWarning( "useMulticoreAutotuning()", "Invalid autotuning option: " << token << ". Switch off autotuning.");
    return false;
  }
}
