#include "GRMHDSolver.h"

#include "GRMHDSolver_Variables.h"

#include "Fortran/PDE.h"
#include "Fortran/InitialData.h"

/****** This is the
 FFFF III N   N III TTTTTT EEEE     V     V  OOO  L    U   U M   M EEEE 
 F     I  NN  N  I    TT   E        V     V O   O L    U   U MM MM E    
 FFF   I  N N N  I    TT   EEE       V   V  O   O L    U   U M M M EEE  
 F     I  N  NN  I    TT   E          V V   O   O L    U   U M   M E    
 F    III N   N III   TT   EEEE        V     OOO  LLLL  UUU  M   M EEEE 
 Version of GRMHD. Don't get confused
 ******/

const double excision_radius = 1.0;

void GRMHD::GRMHDSolver::init(std::vector<std::string>& cmdlineargs) {
  static tarch::logging::Log _log("GRMHDSolver");
  logInfo("init()", "FV Godunov + GRMHD starting.");
}

bool GRMHD::GRMHDSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  using namespace tarch::la;
  // Do excision only in 3D.
  bool insideExcisionBall = norm2(center) < excision_radius;
  //bool insideExcisionBall = false;
  return tarch::la::equals(t,0.0) || insideExcisionBall;

}

void GRMHD::GRMHDSolver::adjustSolution(const double* const x,const double w,const double t,const double dt, double* Q) {
  initialdata_(x, &t, Q);
}

exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void GRMHD::GRMHDSolver::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void GRMHD::GRMHDSolver::flux(const double* const Q, double** F) {
  pdeflux_(F[0], F[1], (DIMENSIONS==3)?F[2]:nullptr, Q);
}


void GRMHD::GRMHDSolver::source(const double* const Q, double* S) {
  pdesource_(S, Q);
  
  for(int l=0; l<NumberOfVariables; l++) {
    S[l] = 0.0;
  }
}

void GRMHD::GRMHDSolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {

  // Exact ID, not time integration
  initialdata_(x, &t, stateOutside);
}

void GRMHD::GRMHDSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
  
  for(int l=0; l<NumberOfVariables; l++) {
    BgradQ[l] = 0.0;
  }
}

// new FV scheme has no coefficient matrix
void GRMHD::GRMHDSolver::coefficientMatrix(const double* const Q,const int d,double* Bn) {
  static tarch::logging::Log _log("GRMHDSolver");
  logError("coefficientMatrix()", "Coefficient Matrix invoked");
  exit(-2);
  /*
  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
  
  for(int l=0; l<NumberOfVariables*NumberOfVariables; l++) {
    Bn[l] = 0.0;
  }
  */
}





