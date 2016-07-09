// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================


#include <memory>

#include "exahype/profilers/Profiler.h"
#include "exahype/solvers/Solver.h"


namespace SRHD{
  class SRHDSolver;
}


class SRHD::SRHDSolver: public exahype::solvers::Solver {
  public:
    SRHDSolver(const std::string& identifier, exahype::solvers::Solver::Type type, int kernelNumber, int numberOfVariables, int numberOfParameters, int nodesPerCoordinateAxis, tarch::la::Vector<DIMENSIONS,double> maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler); 
    int getMinimumTreeDepth() const override;
    void spaceTimePredictor(double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) override; 
    void solutionUpdate(double* luh, const double* const lduh, const double dt) override;
    void volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) override;
    double stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx ) override;
    void solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) override;
    bool hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t) override;
    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) override;
    void faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) override;
    void faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) override;
    void volumeUnknownsProlongation(double* luhFine, const double* luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) override;
    void volumeUnknownsRestriction(double* luhCoarse, const double* luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) override;
  private:
    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);
    static void flux(const double* const Q, double** F);
    static void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);
};


