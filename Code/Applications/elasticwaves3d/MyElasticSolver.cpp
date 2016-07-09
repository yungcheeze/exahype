#include "MyElasticSolver.h"

#include <memory>

elasticwaves3d::MyElasticSolver::MyElasticSolver(const std::string& identifier, exahype::solvers::Solver::Type type, int kernelNumber, int numberOfVariables, int numberOfParameters, int nodesPerCoordinateAxis, tarch::la::Vector<DIMENSIONS,double> maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler)
  : exahype::solvers::Solver(
            identifier, type, kernelNumber, numberOfVariables, numberOfParameters, nodesPerCoordinateAxis, maximumMeshSize, timeStepping, std::move(profiler)) {
  // @todo Please implement/augment if required
}



int elasticwaves3d::MyElasticSolver::getMinimumTreeDepth() const {
  // @todo Please implement
#if defined(Asserts)
  return 2;
#else
  return 4;
#endif
}



void elasticwaves3d::MyElasticSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)

  double* f = F[0];
  double* g = F[1];
  double* h = F[2];

  // material parameters
  double irho = 1.0;        /* density */
  double ilambda = 1.0;     /* 1st Lame parameter */
  double imu = 1.0;         /* 2nd Lame parameter */
  

  // @todo Please implement
  // f
  f[0] = 1.0/irho*Q[3];
  f[1] = 1.0/irho*Q[6];
  f[2] = 1.0/irho*Q[7];
  f[3] = (2.0*imu+ilambda)*Q[0];
  f[4] = ilambda*Q[0];
  f[5] = ilambda*Q[0];
  f[6] = imu*Q[1];
  f[7] = imu*Q[2];
  f[8] = 0.0;
  // g
  // @todo Please implement
  g[0] = 1.0/irho*Q[6];
  g[1] = 1.0/irho*Q[4];
  g[2] = 1.0/irho*Q[8];
  g[3] = ilambda*Q[1];
  g[4] = (2.0*imu+ilambda)*Q[1];
  g[5] = ilambda*Q[1];
  g[6] = imu*Q[0];
  g[7] = 0.0;
  g[8] = imu*Q[2];
  // h
  // @todo Please implement
  h[0] = 1.0/irho*Q[7];
  h[1] = 1.0/irho*Q[8];
  h[2] = 1.0/irho*Q[5];
  h[3] = ilambda*Q[2];
  h[4] = ilambda*Q[2];
  h[5] = (2.0*imu+ilambda)*Q[2];
  h[6] = 0.0;
  h[7] = imu*Q[0];
  h[8] = imu*Q[1];
}



void elasticwaves3d::MyElasticSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)
  // @todo Please implement
  
  // material parameters
  double irho = 1.0;        /* density */
  double ilambda = 1.0;     /* 1st Lame parameter */
  double imu = 1.0;         /* 2nd Lame parameter */
  

  // extract p-wave speed (cp) and s-wave speed (cs)
  double cp = std::sqrt((2.0*imu + ilambda)/irho);
  double cs = std::sqrt(imu/irho);
  
  lambda[0] = -cp;
  lambda[1] = -cs;
  lambda[2] = -cs;
  lambda[3] = cp;
  lambda[4] = cs;
  lambda[5] = cs;
  lambda[6] = 0.0;
  lambda[7] = 0.0;
  lambda[8] = 0.0;
}



bool elasticwaves3d::MyElasticSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  // @todo Please implement
  //return false;
  if (tarch::la::equals(t, 0.0)) {
    return true;
  }
  return false;
}



void elasticwaves3d::MyElasticSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)
  // @todo Please implement

  // initialize particle velocities with a Gaussian
  Q[0] = std::exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)  + (x[2] - 0.5) * (x[2] - 0.5)) /
		  (0.05 * 0.05)) *1.0e-3;
  Q[1] = std::exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)  + (x[2] - 0.5) * (x[2] - 0.5)) /
                 (0.05 * 0.05)) *1.0e-3;
  Q[2] = std::exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)  + (x[2] - 0.5) * (x[2] - 0.5)) /
                 (0.05 * 0.05)) *1.0e-3;
  Q[3] = 0.0;
  Q[4] = 0.0;
  Q[5] = 0.0;
  Q[6] = 0.0;
  Q[7] = 0.0;
  Q[8] = 0.0;
}



exahype::solvers::Solver::RefinementControl elasticwaves3d::MyElasticSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}



