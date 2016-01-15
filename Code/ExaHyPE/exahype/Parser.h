#ifndef EXAHYPE_PARSER
#define EXAHYPE_PARSER

namespace exahype {
  class Parser;
}


#include <vector>
#include "tarch/logging/Log.h"


/**
 * ExaHyPE command line parser
 *
 * @author Tobias Weinzierl
 */
class exahype::Parser {
  private:
    static tarch::logging::Log _log;

    std::vector< std::string > _tokenStream;

    /**
     * @return "notoken" if not found.
     */
    std::string getTokenAfter( std::string token ) const;
    std::string getTokenAfter( std::string token0, std::string token1 ) const;
  public:
    void readFile( const std::string& filename );

    bool isValid() const;

    /**
     * @return How many threads is the code supposed to use?
     */
    int getNumberOfThreads();
};

#endif

