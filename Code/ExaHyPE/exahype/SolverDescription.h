#ifndef _EXAHYPE_SOLVER_DESCRIPTION_H_
#define _EXAHYPE_SOLVER_DESCRIPTION_H_


#include <string>
#include <vector>


namespace exahype {
  class SolverDescription;
}



class exahype::SolverDescription {
  private:
    /**
     * Each solver has an identifier/name. It is used for debug purposes only.
     */
    const std::string _identifier;

    /**
     * Each solver has a kernel number that says which kernel is to be
     * executed. Typically this is an ascending index starting from 0.
     */
    const int         _kernelNumber;
  public:
    static std::vector<SolverDescription> ExistingSolvers;

    SolverDescription(const std::string& identifier, int kernelNumber);

//    kernelIdentifier
};

#endif

