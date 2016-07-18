#include "MyElasticitySolver.h"

#include <memory>

Elasticity::MyElasticitySolver::MyElasticitySolver(int kernelNumber, int nodesPerCoordinateAxis, double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler):
  exahype::solvers::Solver("MyElasticitySolver", exahype::solvers::Solver::Type::ADER_DG, kernelNumber, 10, 0, nodesPerCoordinateAxis, maximumMeshSize, timeStepping, std::move(profiler)) {
  // @todo Please implement/augment if required
}



void Elasticity::MyElasticitySolver::flux(const double* const Q, double** F) {
  // Dimensions             = 2/3
  // Number of variables    = 10 (#unknowns + #parameters)

  double* f = F[0];
  double* g = F[1];

  // material parameters
  // field Q[9] models, in which material parameter regime we are
  double irho =     (Q[9] < 0 ? 1 : 2);  /* density */
  double ilambda =  (Q[9] < 0 ? 1 : 2);  /* 1st Lame parameter */
  double imu =      (Q[9] < 0 ? 1 : 2);  /* 2nd Lame parameter */
  
  // @todo Please implement
  // f
  f[0] = 1.0/irho*Q[2];
  f[1] = 1.0/irho*Q[4];
  f[2] = (2.0*imu +ilambda)*Q[0];
  f[3] = ilambda*Q[0];
  f[4] = imu*Q[1];
  f[5] = 0.0;
  f[6] = 0.0;
  f[7] = 0.0;
  f[8] = 0.0;
  
  // material parameters
  f[ 9] = 0.0;
  
  
  // g
  // @todo Please implement
  g[0] = 1.0/irho*Q[4];
  g[1] = 1.0/irho*Q[3];
  g[2] =  ilambda*Q[1];
  g[3] =  (2.0*imu +ilambda)*Q[1];
  g[4] =  imu *Q[0];
  g[5] = 0.0;
  g[6] = 0.0;
  g[7] = 0.0;
  g[8] = 0.0;
  
  // material parameters
  g[9] = 0.0;
}



void Elasticity::MyElasticitySolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // Dimensions             = 2/3
  // Number of variables    = 10 (#unknowns + #parameters)
  // @todo Please implement

  // material parameters
  // field Q[9] models, in which material parameter regime we are
  double irho =     (Q[9] < 0 ? 1 : 2);  /* density */
  double ilambda =  (Q[9] < 0 ? 1 : 2);  /* 1st Lame parameter */
  double imu =      (Q[9] < 0 ? 1 : 2);  /* 2nd Lame parameter */
  

  // extract p-wave speed (cp) and s-wave speed (cs)
  double cp = std::sqrt((2.0*imu + ilambda)/irho);
  double cs = std::sqrt(imu/irho);
  
  lambda[0] = -cp;
  lambda[1] = -cs;
  lambda[2] = cp;
  lambda[3] = cs;
  lambda[4] = 0.0;
  lambda[5] = 0.0;
  lambda[6] = 0.0;
  lambda[7] = 0.0;
  lambda[8] = 0.0;
  
  // material parameters
  lambda[9] = 0.0;
}



bool Elasticity::MyElasticitySolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  // @todo Please implement
  if (tarch::la::equals(t, 0.0)) {
    return true;
  }
  return false;
}



void Elasticity::MyElasticitySolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 2/3
  // Number of variables    = 10 (#unknowns + #parameters)
  // @todo Please implement
  // initialize particle velocities with a Gaussian
  Q[0] = std::exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)) /
		  (0.05 * 0.05)) *1.0e-3;
  Q[1] = std::exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)) /
                 (0.05 * 0.05)) *1.0e-3;
  
  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
  Q[5] = 0.0;
  Q[6] = 0.0;
  Q[7] = 0.0;
  Q[8] = 0.0;
  
  // material parameters
  // field Q[9] models, in which material parameter regime we are
  if ( (x[0] > 0.8) || (x[0] < 0.2) || (x[1] > 0.8) || (x[1] < 0.2)){
    Q[9] =  1.0;  //hard material
  }else{
    Q[9] = -1.0;  //soft material
  }
}



exahype::solvers::Solver::RefinementControl Elasticity::MyElasticitySolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}



