#include "MySWESolver.h"
#include "InitialData.h"
#include "MySWESolver_Variables.h"

#include "kernels/KernelUtils.h"

using namespace kernels;

const double grav= 9.81;

void SWE::MySWESolver::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
  printf("SWE was called with these parameters:\n");
  for(size_t i=0; i<cmdlineargs.size(); i++)
    printf("%i. %s\n", (int)i, cmdlineargs[i].c_str());

  static tarch::logging::Log _log("MySWESolver::init");
}

SWE::MySWESolver::AdjustSolutionValue SWE::MySWESolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  return tarch::la::equals(t,0.0) ? SWE::MySWESolver::AdjustSolutionValue::PointWisely : SWE::MySWESolver::AdjustSolutionValue::No;
}


void SWE::MySWESolver::adjustPointSolution(const double* const x,const double w,const double t,const double dt,double* Q) {
  assertion(tarch::la::equals(t, 0.0));

  initialData(x,Q);
}


void SWE::MySWESolver::eigenvalues(const double* const Q,const int normalNonZeroIndex,double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 4 (#unknowns + #parameters)
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);  

  if(vars.h()!=vars.h())
     exit(1);

  const double c= std::sqrt(grav*vars.h());
  const double ih = 1./vars.h();

  double u_n = Q[normalNonZeroIndex + 1] * ih;

  eigs.h() = u_n + c ;
  eigs.hu()= u_n -c;
  eigs.hv()= u_n ;

  //TODO: make b constant parameter
  eigs.b() = 0.0;

}


void SWE::MySWESolver::flux(const double* const Q,double** F) {
  ReadOnlyVariables vars(Q);

  const double ih = 1./vars.h();

  double* f = F[0];
  double* g = F[1];

  f[0]= vars.hu();
  f[1]= vars.hu()*vars.hu()*ih + 0.5*grav*vars.h()*vars.h();
  f[2]= vars.hu()*vars.hv()*ih;
  f[3]= 0.0;

  g[0]= vars.hv();
  g[1]= vars.hu()*vars.hv()*ih;
  g[2]= vars.hv()*vars.hv()*ih + 0.5*grav*vars.h()*vars.h();
  g[3]= 0.0;

}


void SWE::MySWESolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int dir,const double * const fluxIn,const double* const stateIn,double *fluxOut,double* stateOut)
 {
  // Dimensions             = 2
  // Number of variables    = 3 (#unknowns + #parameters)

  // Outflow 1
  //fluxOut[0]  = fluxIn[0];
  //fluxOut[1]  = fluxIn[1];
  //fluxOut[2]  = fluxIn[2];

  //Outflow 2 (not working)
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];

  double f[4];
  double g[4];
  double *F[DIMENSIONS];
  F[0] = f; 
  F[1] = g;
  
  F[dir]=fluxOut;
  flux(stateOut,F); 

  //TODO: make outflow work and implement wall boundaries!
}


exahype::solvers::Solver::RefinementControl SWE::MySWESolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo add refinement criterion
  return exahype::solvers::Solver::RefinementControl::Keep;
}

bool SWE::MySWESolver::useNonConservativeProduct() const 
{
   return true;
}

void SWE::MySWESolver::coefficientMatrix(const double* const Q,const int d,double* Bn) {
  // Dimensions             = 2
  // Number of variables    = 3 (#unknowns + #parameters)
  idx2 idx_Bn(NumberOfVariables,NumberOfVariables);  

  Bn[0] = 0.0;
  Bn[1] = 0.0;
  Bn[2] = 0.0;
  Bn[3] = 0.0;
  Bn[4] = 0.0;
  Bn[5] = 0.0;
  Bn[6] = 0.0;
  Bn[7] = 0.0;
  Bn[8] = 0.0;
  Bn[9] = 0.0;
  Bn[10]= 0.0;
  Bn[11]= 0.0;
  Bn[12]= 0.0;
  Bn[13]= 0.0;
  Bn[14]= 0.0;
  Bn[15]= 0.0;

  Bn[idx_Bn(3,d+1)]=grav*Q[0];
}
 
void SWE::MySWESolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  // Dimensions             = 2
  // Number of variables    = 3 (#unknowns + #parameters)
  idx2 idx_gradQ(DIMENSIONS,NumberOfVariables);

  BgradQ[0] = 0.0;
  BgradQ[1] = grav*Q[0]*gradQ[idx_gradQ(0,3)];
  BgradQ[2] = grav*Q[0]*gradQ[idx_gradQ(1,3)]; 
  BgradQ[3] = 0.0;

}



