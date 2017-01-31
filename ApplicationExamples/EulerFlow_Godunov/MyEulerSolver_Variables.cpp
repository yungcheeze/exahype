#include "MyEulerSolver_Variables.h"



EulerFVM::MyEulerSolver::Variables::Variables(double* Q):
  _Q(Q) {
}


double  EulerFVM::MyEulerSolver::Variables::rho() const {
  return _Q[0];
}


double& EulerFVM::MyEulerSolver::Variables::rho() {
  return _Q[0];
}


double  EulerFVM::MyEulerSolver::Variables::E() const {
  return _Q[0+4];
}


double& EulerFVM::MyEulerSolver::Variables::E() {
  return _Q[0+4];
}


double EulerFVM::MyEulerSolver::Variables::v(int index) const {
  return _Q[0+1+index];
}


double& EulerFVM::MyEulerSolver::Variables::v(int index) {
  return _Q[0+1+index];
}


tarch::la::Vector<3,double> EulerFVM::MyEulerSolver::Variables::v() const {
  return tarch::la::Vector<3,double>(_Q+1);
}


void EulerFVM::MyEulerSolver::Variables::v(const tarch::la::Vector<3,double>& values) {
  _Q[1] = values[0];
  _Q[2] = values[1];
  _Q[3] = values[2];
}


void EulerFVM::MyEulerSolver::Variables::v(double v0, double v1, double v2) {
  _Q[1] = v0;
  _Q[2] = v1;
  _Q[3] = v2;
}




/*

class EulerFVM::State {
private:
  double _rho;
  tarch::la::Vector<3,double> _u; // u has dim 3 for 2.5D Euler formulation
  double _E;
public:
  double rho() {
    return _rho;
  }
  const tarch::la::Vector<3,double>& u() {
    return _u;
  }
  double E() {
    return _E;
  }

  State(const double* const Q) {
    _rho  = Q[0];
    for (int i=0; i<3; ++i) {
      _u[i] = Q[i+1];
    }
    _E = Q[1+3];
  }
};

class EulerFVM::Flux {
private:
  double** _F;
  static constexpr int _variables = 5;
public:
  Flux(double** F) : _F(F) {}

  void writeRow(int index, const tarch::la::Vector<2,double>& rowVector) { // We need this method for 2D applications
    assertion2(index>-1,index,_variables);
    assertion1(index<_variables,_variables);
    _F[0][index] = rowVector[0];
    _F[1][index] = rowVector[1];
  }

  void writeRow(int index, const tarch::la::Vector<3,double>& rowVector) { // We need this method for 2.5D applications
    assertion2(index>-1,index,_variables);
    assertion1(index<_variables,_variables);
    _F[0][index] = rowVector[0];
    _F[1][index] = rowVector[1];
#if DIMENSIONS==3
    _F[2][index] = rowVector[2];
#endif
  }
};
*/

