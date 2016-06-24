// ==============================================
// Please do not change the implementations below
// =============================---==============
#include "SRHDSolver.h"
#include "kernels/aderdg/generic/Kernels.h"



void SRHD::SRHDSolver::spaceTimePredictor( double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {
  kernels::aderdg::generic::c::spaceTimePredictorNonlinear<flux>( lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd, luh, dx, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void SRHD::SRHDSolver::solutionUpdate(double* luh, const double* const lduh, const double dt) {
  kernels::aderdg::generic::c::solutionUpdate( luh, lduh, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void SRHD::SRHDSolver::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {
  kernels::aderdg::generic::c::volumeIntegralNonlinear( lduh, lFhi, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void SRHD::SRHDSolver::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {
  kernels::aderdg::generic::c::surfaceIntegralNonlinear( lduh, lFhbnd, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void SRHD::SRHDSolver::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {
  kernels::aderdg::generic::c::riemannSolverNonlinear<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



double SRHD::SRHDSolver::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {
   return kernels::aderdg::generic::c::stableTimeStepSize<eigenvalues>( luh, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void SRHD::SRHDSolver::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  kernels::aderdg::generic::c::solutionAdjustment<adjustedSolutionValues>( luh, center, dx, t, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void SRHD::SRHDSolver::faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {
  kernels::aderdg::generic::c::faceUnknownsProlongation( lQhbndFine, lFhbndFine, lQhbndCoarse, lFhbndCoarse, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void SRHD::SRHDSolver::faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {
  kernels::aderdg::generic::c::faceUnknownsRestriction( lQhbndCoarse, lFhbndCoarse, lQhbndFine, lFhbndFine, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void SRHD::SRHDSolver::volumeUnknownsProlongation(double* luhFine, const double* luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  kernels::aderdg::generic::c::volumeUnknownsProlongation( luhFine, luhCoarse, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void SRHD::SRHDSolver::volumeUnknownsRestriction(double* luhCoarse, const double* luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
  kernels::aderdg::generic::c::volumeUnknownsRestriction( luhCoarse, luhFine, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



