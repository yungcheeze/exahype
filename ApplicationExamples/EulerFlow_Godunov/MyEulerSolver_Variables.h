#ifndef __MyEulerSolver_Variables_CLASS_HEADER__
#define __MyEulerSolver_Variables_CLASS_HEADER__

#include "MyEulerSolver.h"

#include "tarch/la/Matrix.h"

class EulerFVM::MyEulerSolver::Variables {
  private:
    double* _Q;
  public:
    Variables(double* Q) : _Q(Q) {}

    double rho() const { return _Q[0]; }

    double u(int index) const {
      assertion(index >= 0 && index<3);
      return _Q[1+index];
    }

    tarch::la::Vector<3,double> u() const {
      tarch::la::Vector<3,double> values(_Q[1],_Q[2],_Q[3]);
      return values;
    }

    double E() const { return _Q[4]; }



    double& rho() { return _Q[0]; }

    double& u(int index) { return _Q[1+index]; }

    void u(const tarch::la::Vector<3,double>& values) {
      *(_Q+1)=values[0];
      *(_Q+2)=values[1];
      *(_Q+3)=values[2];
    }

    void u(double u0,double u1,double u2) {
      *(_Q+1)=u0;
      *(_Q+2)=u1;
      *(_Q+3)=u2;
    }

    double& E() { return _Q[4]; }


};


class EulerFVM::MyEulerSolver::ReadOnlyVariables {
  private:
    const double* const _Q;
  public:
    ReadOnlyVariables(const double* const Q) : _Q(Q) {}

    double rho() const { return _Q[0]; }

    double u(int index) const {
      assertion(index >= 0 && index<3);
      return _Q[1+index];
    }

    tarch::la::Vector<3,double> u() const {
      tarch::la::Vector<3,double> values(_Q[1],_Q[2],_Q[3]);
      return values;
    }

    double E() const { return _Q[4]; }


};


class EulerFVM::MyEulerSolver::Fluxes {
  private:
    double** _F;
  public:
    Fluxes(double** F) : _F(F) {}

    double rho(int column) const {
      assertion(column >= 0 && column<DIMENSIONS);
      return _F[column][0];
    }

    tarch::la::Vector<DIMENSIONS,double> rho() const {
      tarch::la::Vector<DIMENSIONS,double> values(_F[0][0],_F[1][0]);
      return values;
    }

    double u(int row, int column) const {
      assertion(row >= 0 && row<3);
      assertion(column >= 0 && column<DIMENSIONS);
      return _F[column][1+row];
    }

    tarch::la::Vector<DIMENSIONS,double> u(int row) const {
      assertion(row >= 0 && row<3);
      tarch::la::Vector<DIMENSIONS,double> values(_F[0][1+row],_F[1][1+row]);
      return values;
    }

    tarch::la::Matrix<3,DIMENSIONS,double> u() const {
      tarch::la::Matrix<3,DIMENSIONS,double> values;
      values = _F[0][1],_F[1][1],
               _F[0][2],_F[1][2],
               _F[0][3],_F[1][3];
      return values;
    }

    double E(int column) const {
      assertion(column >= 0 && column<DIMENSIONS);
      return _F[column][4];
    }

    tarch::la::Vector<DIMENSIONS,double> E() const {
      tarch::la::Vector<DIMENSIONS,double> values(_F[0][4],_F[1][4]);
      return values;
    }



    double& rho(int column) {
      assertion(column >= 0 && column<DIMENSIONS);
      return _F[column][0];
    }

    void rho(const tarch::la::Vector<DIMENSIONS,double>& values) {
      _F[0][0]=values[0];
      _F[1][0]=values[1];
      #if DIMENSIONS==3
      _F[2][0]=values[2];
      #endif
    }
    #if DIMENSIONS==2
    /** Setter for 2.5D calculations. Third vector element is ignored.*/
    void rho(const tarch::la::Vector<3,double>& values) {
      _F[0][0]=values[0];
      _F[1][0]=values[1];
    }
    #endif

    /** Setter for 3D and 2.5D calculations. Third argument is ignored for the latter.*/
    void rho(double v0,double v1,double v2) {
      _F[0][0]=v0;
      _F[1][0]=v1;
      #if DIMENSIONS==3
      _F[2][0]=v2;
      #endif
    }
    #if DIMENSIONS==2
    void rho(double v0,double v1) {
      _F[0][0]=v0;
      _F[1][0]=v1;
    }
    #endif

    double& u(int row, int column) {
      assertion(row >= 0 && row<3);
      assertion(column >= 0 && column<DIMENSIONS);
      return _F[column][1+row];
    }

    void u(int row, const tarch::la::Vector<DIMENSIONS,double>& values) {
      assertion(row >= 0 && row<3);
      _F[0][1+row]=values[0];
      _F[1][1+row]=values[1];
      #if DIMENSIONS==2
      _F[2][1+row]=values[2];
      #endif
    }
    #if DIMENSIONS==2
    /** Setter for 2.5D calculations. Third vector element is ignored.*/
    void u(int row, const tarch::la::Vector<3,double>& values) {
      assertion(row >= 0 && row<3);
      _F[0][1+row]=values[0];
      _F[1][1+row]=values[1];
    }
    #endif

    void u(const tarch::la::Matrix<3,DIMENSIONS,double>& values) {
      _F[0][1]=values(0,0);
      _F[0][2]=values(1,0);
      _F[0][3]=values(2,0);
      _F[1][1]=values(0,1);
      _F[1][2]=values(1,1);
      _F[1][3]=values(2,1);
      #if DIMENSIONS==3
      _F[2][1]=values(0,2);
      _F[2][2]=values(1,2);
      _F[2][3]=values(2,2);
      #endif
    }
    #if DIMENSIONS==2
    /** Setter for 2.5D calculations. Third matrix column is ignored.*/
    void u(const tarch::la::Matrix<3,3,double>& values) {
      _F[0][1]=values(0,0);
      _F[0][2]=values(1,0);
      _F[0][3]=values(2,0);
      _F[1][1]=values(0,1);
      _F[1][2]=values(1,1);
      _F[1][3]=values(2,1);
    }
    #endif

    /** Setter for 3D and 2.5D calculations. Third argument is ignored for the latter.*/
    void u(int row, double v0,double v1,double v2) {
      assertion(row >= 0 && row<3);
      _F[0][1+row]=v0;
      _F[1][1+row]=v1;
      #if DIMENSIONS==3
      _F[2][1+row]=v2;
      #endif
    }
    #if DIMENSIONS==2
    /** Setter for 2D calculations.*/
    void u(int row, double v0,double v1) {
      assertion(row >= 0 && row<3);
      _F[0][1+row]=v0;
      _F[1][1+row]=v1;
    }
    #endif

    /** Setter for 3D and 2.5D calculations. Third column values are ignored for the latter.*/
    void u(double v00, double v01, double v02,
           double v10, double v11, double v12,
           double v20, double v21, double v22) {
      _F[0][1]=v00;
      _F[0][2]=v10;
      _F[0][3]=v20;
      _F[1][1]=v01;
      _F[1][2]=v11;
      _F[1][3]=v21;
      #if DIMENSIONS==3
      _F[2][1]=v02;
      _F[2][2]=v12;
      _F[2][3]=v22;
      #endif
    }
    #if DIMENSIONS==2
    void u(double v00, double v01,
           double v10, double v11,
           double v20, double v21) {
      _F[0][1]=v00;
      _F[0][2]=v10;
      _F[0][3]=v20;
      _F[1][1]=v01;
      _F[1][2]=v11;
      _F[1][3]=v21;
    }
    #endif

    double& E(int column) {
      assertion(column >= 0 && column<DIMENSIONS);
      return _F[column][4];
    }

    void E(const tarch::la::Vector<DIMENSIONS,double>& values) {
      _F[0][4]=values[0];
      _F[1][4]=values[1];
      #if DIMENSIONS==3
      _F[2][4]=values[2];
      #endif
    }
    #if DIMENSIONS==2
    /** Setter for 2.5D calculations. Third vector element is ignored.*/
    void E(const tarch::la::Vector<3,double>& values) {
      _F[0][4]=values[0];
      _F[1][4]=values[1];
    }
    #endif

    /** Setter for 3D and 2.5D calculations. Third argument is ignored for the latter.*/
    void E(double v0,double v1,double v2) {
      _F[0][4]=v0;
      _F[1][4]=v1;
      #if DIMENSIONS==3
      _F[2][4]=v2;
      #endif
    }
    #if DIMENSIONS==2
    void E(double v0,double v1) {
      _F[0][4]=v0;
      _F[1][4]=v1;
    }
    #endif


};

#endif