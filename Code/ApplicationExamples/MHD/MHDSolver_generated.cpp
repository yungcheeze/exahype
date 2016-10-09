// ==============================================
// Please do not change the implementations below
// =============================---==============
#include "MHDSolver.h"
#include "kernels/aderdg/generic/Kernels.h"


MHDSolver::MHDSolver::MHDSolver(double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler, exahype::Parser::ParserView constants):
  exahype::solvers::ADERDGSolver("MHDSolver", 9 /* numberOfUnknowns */, 0/* numberOfParameters */, 4 /* nodesPerCoordinateAxis */, maximumMeshSize, timeStepping, std::move(profiler)) {
  init();
  MHDSolver::MHDSolver::constants = &constants;
}



void MHDSolver::MHDSolver::spaceTimePredictor(double* lQhbnd,double* lFhbnd,double** tempSpaceTimeUnknowns,double** tempSpaceTimeFluxUnknowns,double*  tempUnknowns,double*  tempFluxUnknowns,const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt) {
  kernels::aderdg::generic::c::spaceTimePredictorNonlinear<flux, source>( lQhbnd, lFhbnd, tempSpaceTimeUnknowns, tempSpaceTimeFluxUnknowns, tempUnknowns, tempFluxUnknowns, luh, dx, dt, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );
}



void MHDSolver::MHDSolver::solutionUpdate(double* luh, const double* const lduh, const double dt) {
  kernels::aderdg::generic::c::solutionUpdate( luh, lduh, dt, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );
}



void MHDSolver::MHDSolver::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {
  kernels::aderdg::generic::c::volumeIntegralNonlinear( lduh, lFhi, dx, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );
}



void MHDSolver::MHDSolver::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {
  kernels::aderdg::generic::c::surfaceIntegralNonlinear( lduh, lFhbnd, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void MHDSolver::MHDSolver::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, double* tempFaceUnknownsArray, double** tempStateSizedVectors, double** tempStateSizedSquareMatrices, const double dt, const int normalNonZeroIndex) {
  assertion2(normalNonZeroIndex>=0,dt,normalNonZeroIndex);
  assertion2(normalNonZeroIndex<DIMENSIONS,dt,normalNonZeroIndex);
  kernels::aderdg::generic::c::riemannSolverNonlinear<eigenvalues>( FL, FR, QL, QR, tempFaceUnknownsArray, tempStateSizedVectors, tempStateSizedSquareMatrices, dt, normalNonZeroIndex, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );
}



void MHDSolver::MHDSolver::boundaryConditions(double* fluxOut,double* stateOut,const double* const fluxIn,const double* const stateIn,const tarch::la::Vector<DIMENSIONS, double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const int faceIndex,const int normalNonZero) {
  kernels::aderdg::generic::c::boundaryConditions<boundaryValues>( fluxOut, stateOut, fluxIn, stateIn, cellCentre, cellSize, t, dt, faceIndex, normalNonZero, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



double MHDSolver::MHDSolver::stableTimeStepSize(const double* const luh,double* tempEigenvalues,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  double d = kernels::aderdg::generic::c::stableTimeStepSize<eigenvalues>( luh, tempEigenvalues, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );
  return d;
}



void MHDSolver::MHDSolver::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  kernels::aderdg::generic::c::solutionAdjustment<adjustedSolutionValues>( luh, center, dx, t, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void MHDSolver::MHDSolver::faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {
  kernels::aderdg::generic::c::faceUnknownsProlongation( lQhbndFine, lFhbndFine, lQhbndCoarse, lFhbndCoarse, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void MHDSolver::MHDSolver::faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {
  kernels::aderdg::generic::c::faceUnknownsRestriction( lQhbndCoarse, lFhbndCoarse, lQhbndFine, lFhbndFine, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void MHDSolver::MHDSolver::volumeUnknownsProlongation(double* luhFine, const double* luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  kernels::aderdg::generic::c::volumeUnknownsProlongation( luhFine, luhCoarse, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void MHDSolver::MHDSolver::volumeUnknownsRestriction(double* luhCoarse, const double* luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  kernels::aderdg::generic::c::volumeUnknownsRestriction( luhCoarse, luhFine, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



