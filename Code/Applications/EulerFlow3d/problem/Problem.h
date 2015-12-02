/**
 * Declares flux functions for a elementwise constant DG method.
 */

#ifndef EXAHYPE_EULERFLOW3D_PROBLEM_H_
#define EXAHYPE_EULERFLOW3D_PROBLEM_H_

namespace exahype {
  namespace problem {
    /**
     * @brief Computes the initial value of the problem at point (x,y)
     *
     * @param[in] x physical x coordinate
     * @param[in] y physical y coordinate
     * @param[in] numberOfVariables the number of consered quantities
     * @param[out] values array of size numberOfVariables
     */
    void PDEInitialValue2d(const double x,const double y,const int numberOfVariables,double * value);
    double PDEInflow(const double x,const double y);

    double PDEVolumeFlux(double x, double y);
    double PDENormalFlux(const double x,const double y,const double nx,const double ny);

    double DGNormalFlux(const double x,const double y,const double nx,const double ny);
    void   DGRiemannSolver(const double x,const double y,const double nx,const double ny,double *selfFlux,double *neighbourFlux);
  }  // namespace problem
}  // namespace exahype

#endif /* EXAHYPE_EULERFLOW3D_PROBLEM_H_ */
