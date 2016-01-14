#include "exahype/Parser.h"

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
        logInfo( "readFile(String)", "got token " << newToken );
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



/*
string name,age,salary,hoursWorked,randomText;

while(getline(readFile,line))   {
    stringstream iss(line);
    getline(iss, name, ':');
    getline(iss, age, '-');
    getline(iss, salary, ',');
    getline(iss, hoursWorked, '[');
    getline(iss, randomText, ']');
}
readFile.close();


a




using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <cstring>
*/
/*

const int MAX_CHARS_PER_LINE = 512;
const int MAX_TOKENS_PER_LINE = 20;
const char* const DELIMITER = " ";

int main()
{
  // create a file-reading object
  ifstream fin;
  fin.open("data.txt"); // open a file
  if (!fin.good())
    return 1; // exit if file not found

  // read each line of the file
  while (!fin.eof())
  {
    // read an entire line into memory
    char buf[MAX_CHARS_PER_LINE];
    fin.getline(buf, MAX_CHARS_PER_LINE);

    // parse the line into blank-delimited tokens
    int n = 0; // a for-loop index

    // array to store memory addresses of the tokens in buf
    const char* token;

    // parse the line
    token = strtok(buf, DELIMITER); // first token
    if (token) // zero if line is blank
    {
      for (n = 1; n < MAX_TOKENS_PER_LINE; n++)
      {
        token[n] = strtok(0, DELIMITER); // subsequent tokens
        if (!token[n]) break; // no more tokens
      }
    }

    // process (print) the tokens
    for (int i = 0; i < n; i++) // n = #of tokens
      cout << "Token[" << i << "] = " << token[i] << endl;
    cout << endl;
  }
}


*/
