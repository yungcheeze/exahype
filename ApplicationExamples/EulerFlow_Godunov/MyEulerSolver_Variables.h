

// @todo Has to be included by generator
#include "MyEulerSolver.h"

#include "tarch/la/Matrix.h"


/**
 * @todo: Diese Klasse sollte Variables heissen. Variables selber gibt es auch, aber das darf dann nur getter offerieren. Es kann dafuer mit const double const * arbieten
 * @todo Entweder RWVariables order ReadOnlyVariables
 */
class EulerFVM::MyEulerSolver::Variables {
  private:
    double* _Q;
  public:
    Variables(double* Q);

    double  rho() const;
    double& rho();

    double  E() const;
    double& E();

    double v(int index) const;
    tarch::la::Vector<3,double> v() const;

    /**
     * Jetzt kommen die Vektor Writer und hier wird's spassig
     */
    double& v(int index);
    void v(const tarch::la::Vector<3,double>& values);
    void v(double v0, double v1, double v2);
};


class EulerFVM::MyEulerSolver::Fluxes {
  private:
    double** _F;
  public:
    Fluxes(double** F);

    tarch::la::Vector<3,double>  rho() const;
    double& rho(int index);

    tarch::la::Vector<3,double> E() const;
    double& E(int index);

    tarch::la::Matrix<3,3,double> v() const;
    tarch::la::Vector<3,double> v(int index) const;

    /**
     * Jetzt kommen die Vektor Writer und hier wird's spassig
     */
    double& v(int index);
    void v(const tarch::la::Vector<3,double>& values);
    void v(double v0, double v1, double v2);
};

