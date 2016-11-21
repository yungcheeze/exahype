/**
 * This is the fancy FORTRAN and C calling kernel generated stuff.
 *
 **/
#include "MHDSolver.h"

// includes fortran and c kernels
#include "kernels/aderdg/generic/Kernels.h"

// std::abs
#include <climits>
#include <cstdlib>


MHDSolver::MHDSolver::MHDSolver(double maximumMeshSize,exahype::solvers::Solver::TimeStepping timeStepping,std::vector<std::string>& cmdlineargs, exahype::Parser::ParserView constants):
  exahype::solvers::ADERDGSolver("MHDSolver",nVar /* numberOfUnknowns */,nParams /* numberOfParameters */,order + 1  /* nodesPerCoordinateAxis */,maximumMeshSize,timeStepping) {
  init(cmdlineargs, constants);
}



void MHDSolver::MHDSolver::spaceTimePredictor(double* lQhbnd,double* lFhbnd,double** tempSpaceTimeUnknowns,double** tempSpaceTimeFluxUnknowns,double* tempUnknowns,double* tempFluxUnknowns,double* tempStateSizedVectors,const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt) {
  kernels::aderdg::generic::c::spaceTimePredictorNonlinear<MHDSolver>(*this,lQhbnd,lFhbnd,tempSpaceTimeUnknowns,tempSpaceTimeFluxUnknowns,tempUnknowns,tempFluxUnknowns,tempStateSizedVectors,luh,dx,dt);
}



void MHDSolver::MHDSolver::solutionUpdate(double* luh,const double* const lduh,const double dt) {
  kernels::aderdg::generic::c::solutionUpdate(luh,lduh,dt,getNumberOfVariables(),getNumberOfParameters(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::volumeIntegral(double* lduh,const double* const lFhi,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  kernels::aderdg::generic::c::volumeIntegralNonlinear(lduh,lFhi,dx,getNumberOfVariables(),getNumberOfParameters(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::surfaceIntegral(double* lduh,const double* const lFhbnd,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  kernels::aderdg::generic::c::surfaceIntegralNonlinear(lduh,lFhbnd,dx,getNumberOfVariables(),getNodesPerCoordinateAxis());
}

int globalRiemannSolverCounter(0);

void MHDSolver::MHDSolver::riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,double* tempFaceUnknownsArray,double** tempStateSizedVectors,double** tempStateSizedSquareMatrices,const double dt,const int normalNonZeroIndex) {
  assertion2(normalNonZeroIndex>=0,dt,normalNonZeroIndex);
  assertion2(normalNonZeroIndex<DIMENSIONS,dt,normalNonZeroIndex);
  const int basissize = order + 1;
  globalRiemannSolverCounter++;

  // First call Fortran as we know it works and doesn't corrupt data.
  // output quantities
  double fort_FL[basissize*basissize*nVar];
  double fort_FR[basissize*basissize*nVar];
  // input quantities. May be copied for security. Not used do far.
  double fort_QL[basissize*basissize*nVar];
  double fort_QR[basissize*basissize*nVar];
  // temporary storage for computation
  double *fort_tempFaceUnknownsArray = NULL; // This is never used.
  double *fort_tempStateSizedVectors[nDim];
  double v1[nVar], v2[nVar], v3[nVar];
  fort_tempStateSizedVectors[0] = v1;
  fort_tempStateSizedVectors[1] = v2;
  if(nDim>2) fort_tempStateSizedVectors[2] = v3;
  double **fort_tempStateSizedSquareMatrices = NULL; // Never used

  kernels::aderdg::generic::fortran::riemannSolverNonlinear<MHDSolver>(*this,fort_FL,fort_FR,QL,QR,fort_tempFaceUnknownsArray,fort_tempStateSizedVectors,fort_tempStateSizedSquareMatrices,dt,normalNonZeroIndex);

  // then call C on the actual data
  kernels::aderdg::generic::c::riemannSolverNonlinear<MHDSolver>(*this,FL,FR,QL,QR,tempFaceUnknownsArray,tempStateSizedVectors,tempStateSizedSquareMatrices,dt,normalNonZeroIndex);

  // now compare data
  const double accepterror = 1e-10;
  for(int i=0; i<basissize*basissize*nVar; i++) {
    double errorFL = std::abs(fort_FL[i] - FL[i]);
    double errorFR = std::abs(fort_FR[i] - FR[i]);
    if(errorFL > accepterror) printf("%d. FL[%d]=%+e fort_FL[%d]=%+e errorFL=%e\n", globalRiemannSolverCounter, i, FL[i], i, fort_FL[i], errorFL);
    if(errorFR > accepterror) printf("%d. FR[%d]=%+e fort_FR[%d]=%+e errorFR=%e\n", globalRiemannSolverCounter, i, FR[i], i, fort_FR[i], errorFR);
  }
}



void MHDSolver::MHDSolver::boundaryConditions(double* fluxOut,double* stateOut,const double* const fluxIn,const double* const stateIn,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const int faceIndex,const int normalNonZero) {
  kernels::aderdg::generic::c::boundaryConditions<MHDSolver>(*this,fluxOut,stateOut,fluxIn,stateIn,cellCentre,cellSize,t,dt,faceIndex,normalNonZero);
}



double MHDSolver::MHDSolver::stableTimeStepSize(const double* const luh,double* tempEigenvalues,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  double d = kernels::aderdg::generic::c::stableTimeStepSize<MHDSolver>(*this,luh,tempEigenvalues,dx);
  return d;
}



void MHDSolver::MHDSolver::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  kernels::aderdg::generic::c::solutionAdjustment<MHDSolver>(*this,luh,center,dx,t,dt);
}



void MHDSolver::MHDSolver::faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex) {
  kernels::aderdg::generic::c::faceUnknownsProlongation(lQhbndFine,lFhbndFine,lQhbndCoarse,lFhbndCoarse,coarseGridLevel,fineGridLevel,subfaceIndex,getNumberOfVariables(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex) {
  kernels::aderdg::generic::c::faceUnknownsRestriction(lQhbndCoarse,lFhbndCoarse,lQhbndFine,lFhbndFine,coarseGridLevel,fineGridLevel,subfaceIndex,getNumberOfVariables(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::volumeUnknownsProlongation(double* luhFine,const double* luhCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  kernels::aderdg::generic::c::volumeUnknownsProlongation(luhFine,luhCoarse,coarseGridLevel,fineGridLevel,subcellIndex,getNumberOfVariables(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::volumeUnknownsRestriction(double* luhCoarse,const double* luhFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  kernels::aderdg::generic::c::volumeUnknownsRestriction(luhCoarse,luhFine,coarseGridLevel,fineGridLevel,subcellIndex,getNumberOfVariables(),getNodesPerCoordinateAxis());
}