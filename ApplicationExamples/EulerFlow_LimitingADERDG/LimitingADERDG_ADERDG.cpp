#include "LimitingADERDG_ADERDG.h"

#include "InitialData.h"
#include "LimitingADERDG_ADERDG_Variables.h"

#include <cstring>

void Euler::LimitingADERDG_ADERDG::init(std::vector<std::string>& cmdlineargs) {
  // This function is called inside the constructur.
  // @todo Please implement/augment if required.
}

void Euler::LimitingADERDG_ADERDG::flux(const double* const Q, double** F) {
  ReadOnlyVariables vars(Q);
  Fluxes f(F);

  tarch::la::Matrix<3,3,double> I;
  I = 1, 0, 0,
      0, 1, 0,
      0, 0, 1;

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  f.rho ( vars.j()                                 );
  f.j   ( irho * outerDot(vars.j(),vars.j()) + p*I );
  f.E   ( irho * (vars.E() + p) * vars.j()         );
}


void Euler::LimitingADERDG_ADERDG::algebraicSource(const double* const Q, double* S) {
  Variables s(S);
  s.rho()=0;
  s.E()=0;
  s.j(0,0,0);
}


void Euler::LimitingADERDG_ADERDG::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  double u_n = vars.j(normalNonZeroIndex) * irho;
  double c  = std::sqrt(GAMMA * p * irho);

  eigs.rho()=u_n - c;
  eigs.E()  =u_n + c;
  eigs.j(u_n,u_n,u_n);
}


void Euler::LimitingADERDG_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double * const fluxIn,const double* const stateIn,double *fluxOut,double* stateOut) {
  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];
  fluxOut[3] = fluxIn[3];
  fluxOut[4] = fluxIn[4];

  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];
}


bool Euler::LimitingADERDG_ADERDG::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t, double dt) {
  return tarch::la::equals(t, 0.0);
}

void Euler::LimitingADERDG_ADERDG::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  // @todo Please implement
  Euler::initialData(x,Q);
}

void Euler::LimitingADERDG_ADERDG::nonConservativeProduct(const double* const Q, const double* const gradQ, double* BgradQ) {
  std::memset(BgradQ, 0, NumberOfVariables * sizeof(double));
}

void Euler::LimitingADERDG_ADERDG::coefficientMatrix(const double* const Q, const int normalNonZero, double* Bn) {
  std::memset(Bn, 0, NumberOfVariables * NumberOfVariables * sizeof(double));
}

exahype::solvers::Solver::RefinementControl Euler::LimitingADERDG_ADERDG::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}

bool Euler::LimitingADERDG_ADERDG::physicalAdmissibilityDetection(const double* const QMin,const double* const QMax) {
  if (QMin[0] < 0.0) return false;
  if (QMin[4] < 0.0) return false;

  for (int i=0; i<5; ++i) {
    if (!std::isfinite(QMin[i])) return false;
    if (!std::isfinite(QMax[i])) return false;
  }

  return true;
}

bool Euler::LimitingADERDG_ADERDG::useAlgebraicSource() const {return true;}

bool Euler::LimitingADERDG_ADERDG::useNonConservativeProduct() const {return true;}

bool Euler::LimitingADERDG_ADERDG::useCoefficientMatrix() const {return true;}

bool Euler::LimitingADERDG_ADERDG::usePointSource() const {
  return false;
}

void Euler::LimitingADERDG_ADERDG::pointSource(const double* const x,const double t,const double dt, double* forceVector, double* x0) {
  // do nothing
}
