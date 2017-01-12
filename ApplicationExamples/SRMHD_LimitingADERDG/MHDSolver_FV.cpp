#include "MHDSolver_FV.h"

#include "PDE.h"

#include <memory>

void MHD::MHDSolver_FV::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

bool MHD::MHDSolver_FV::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) {
  return (t < 1e-10);
}

void MHD::MHDSolver_FV::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
    // Fortran call:
  adjustedsolutionvalues_(x, &w, &t, &dt, Q);
}

exahype::solvers::Solver::RefinementControl MHD::MHDSolver_FV::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void MHD::MHDSolver_FV::eigenvalues(const double* const Q,const int normalNonZeroIndex,double* lambda) {
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void MHD::MHDSolver_FV::flux(const double* const Q,double** F) {
  // Caveats: Fortran accepts a uniform array of size (nVar*nDim), however C passes an array of pointers.
  // This Fortran interface works only if F is a continous array and F[1]==F[nDim+1] etc!
  pdeflux_(F[0], Q);
}


void MHD::MHDSolver_FV::source(const double* const Q,double* S) {
  constexpr int nVar = 9;

  for(int i=0; i < nVar; i++) {
    S[i] = 0.0;
  }
}


void MHD::MHDSolver_FV::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
const double* const stateIn,double* stateOut) {
  constexpr int nVar = 9;
  
  for (int i=0; i<nVar; ++i) {
    stateOut[i] = stateIn[i];
  }

//  const double* normalVelocity = stateIn+1+normalNonZero;
//  const double sign = (faceIndex-2*normalNonZero)==0 ? -1.0 : 1.0;
//  const double outwardDirectedVelocity = sign * (*normalVelocity);
//
//  std::cout << "normalVelocity="<<normalVelocity<<",faceIndex="<<faceIndex<<",outwardDirectedVelocity="<<outwardDirectedVelocity << std::endl;
//
//  if (outwardDirectedVelocity>0) { // outflow; take inside value
//    // do nothing
//  } else { // inflow; take outside values
//    stateOut[0] = 1.0e-4; // density for InitialBlast
//    stateOut[4] = 5.0e-4; // pressure for InitialBlast
//  }
}


//void MHD::MHDSolver_FV::ncp(const double* const Q,const double* const gradQ,double* BgradQ) {
//  constexpr int nVar = 9;
//  std::memset(BgradQ, 0, nVar * sizeof(double));
//}
//
//
//void MHD::MHDSolver_FV::matrixb(const double* const Q,const int normalNonZero,double* Bn) {
//  constexpr int nVar = 9;
//  std::memset(Bn, 0, nVar * nVar * sizeof(double));
//}
