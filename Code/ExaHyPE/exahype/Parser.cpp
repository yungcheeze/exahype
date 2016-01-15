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


std::string exahype::Parser::getTokenAfter( std::string token ) const {
  assertion( isValid() );
  int currentToken = 0;
  while (_tokenStream[currentToken]!=token && currentToken<static_cast<int>(_tokenStream.size())) {
    currentToken++;
  }
  currentToken++;
  if ( currentToken<static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  }
  else return "notoken";
}


std::string exahype::Parser::getTokenAfter( std::string token0, std::string token1 ) const {
  assertion( isValid() );
  int currentToken = 0;
  while (_tokenStream[currentToken]!=token0 && currentToken<static_cast<int>(_tokenStream.size())) {
    currentToken++;
  }
  while (_tokenStream[currentToken]!=token1 && currentToken<static_cast<int>(_tokenStream.size())) {
    currentToken++;
  }
  currentToken++;
  if ( currentToken<static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  }
  else return "notoken";
}


int exahype::Parser::getNumberOfThreads() {
  assertion( isValid() );
  std::string token = getTokenAfter("shared-memory","cores");
  logInfo( "getNumberOfThreads()", "found token " << token );
  int result = atoi( token.c_str() );
  if (result==0) {
    logInfo( "getNumberOfThreads()", "Invalid number of cores set: " << token << ". Use one core, i.e. switch off multithreading" );
    result = 1;
  }
  return result;
}
