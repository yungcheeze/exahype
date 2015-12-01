/**
 * Declares flux functions for a elementwise constant DG method.
 */

#ifndef EXAHYPE_ADVECTIONDG0_PROBLEM_H_
#define EXAHYPE_ADVECTIONDG0_PROBLEM_H_

namespace exahype {
  namespace problem {
    double PDEInitialValue(const double x,const double y);
    double PDEInflow(const double x,const double y);

    double PDEVolumeFlux(double x, double y);
    double PDENormalFlux(const double x,const double y,const double nx,const double ny);

    double DGNormalFlux(const double x,const double y,const double nx,const double ny);
    void   DGRiemannSolver(const double x,const double y,const double nx,const double ny,double *selfFlux,double *neighbourFlux);
  }  // namespace problem
}  // namespace exahype

#endif /* EXAHYPE_ADVECTIONDG0_PROBLEM_H_ */
