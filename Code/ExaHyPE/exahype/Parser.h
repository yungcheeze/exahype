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
  public:
    void readFile( const std::string& filename );

    bool isValid() const;
};

#endif

