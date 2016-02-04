// ==============================================
// Please do not change the implementations below
// =============================---==============
#include "MyEulerSolver.h"
#include "kernels/aderdg/generic/Kernels.h"



#include "kernels/aderdg/generic/Kernels.h"


void Euler2d::MyEulerSolver::spaceTimePredictor( double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {
   kernels::aderdg::generic::spaceTimePredictor<flux>( lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd, luh, dx, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void Euler2d::MyEulerSolver::solutionUpdate(double* luh, const double* const lduh, const double dt) {
   kernels::aderdg::generic::solutionUpdate( luh, lduh, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void Euler2d::MyEulerSolver::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {
   kernels::aderdg::generic::volumeIntegral( lduh, lFhi,dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void Euler2d::MyEulerSolver::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {
   kernels::aderdg::generic::surfaceIntegral( lduh, lFhbnd, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void Euler2d::MyEulerSolver::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {
   kernels::aderdg::generic::riemannSolver<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



double Euler2d::MyEulerSolver::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {
   return kernels::aderdg::generic::stableTimeStepSize<eigenvalues>( luh, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



void Euler2d::MyEulerSolver::initialCondition(double* luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx) {
   kernels::aderdg::generic::initialCondition<initialValues>( luh, center, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );
}



