#ifndef _EXAHYPE_SOLVERS_SOLVER_H_
#define _EXAHYPE_SOLVERS_SOLVER_H_


#include <string>
#include <vector>


#include "peano/utils/Globals.h"
#include "tarch/la/Vector.h"


namespace exahype {
  namespace solvers {
    class Solver;

    extern std::vector<Solver*> RegisteredSolvers;
  }
}

/**
 * Describes one solver.
 */
class exahype::solvers::Solver {
public:
  struct Plot {
    int          variable;
    double       nextSnapshot;
    bool         repeat;
    std::string  filename;

    Plot( int variable_, double nextSnapshot_, bool repeat_, const std::string& filename_);
  };

  enum Type {
    ADER_DG
  };
protected:
  /**
   * Each solver has an identifier/name. It is used for debug purposes only.
   */
  const std::string _identifier;

  const Type        _type;

  /**
   * Each solver has a kernel number that says which kernel is to be
   * executed. Typically this is an ascending index starting from 0.
   */
  const int         _kernelNumber;

  const int         _numberOfVariables;

  const int         _nodesPerCoordinateAxis;
public:
  Solver(const std::string& identifier, Type type, int kernelNumber, int numberOfVariables, int nodesPerCoordinateAxis);

  virtual ~Solver() {}

  /**
   * Identify minimal mesh width at a certain point in the domain. This
   * minimal mesh width is used both as a constraint on the AMR as well
   * as to set up the initial grid. If you return 0, you indicate that
   * this PDE might not exist in the domain.
   */
  virtual int getMinimumTreeDepth() const = 0;

  Type getType() const;

  int getNumberOfVariables() const;

  /**
   * If you use a higher order method, then this operation returns the
   * polynomial degree plus one. If you use a Finite Volume method, it
   * returns the number of cells within a patch per coordinate axis.
   */
  int getNodesPerCoordinateAxis() const;

  /**
   * @brief Adds the solution update to the solution.
   *
   * @param[inout] luh  Cell-local solution DoF.
   * @param[in]    lduh Cell-local update DoF.
   * @param[in]    dx   Extent of the cell in each coordinate direction.
   * @param[dt]    dt   Time step size.
   */
  virtual void solutionUpdate(
      double * luh,
      const double * const lduh,
      const tarch::la::Vector<DIMENSIONS,double>&  dx,
      const double dt
  ) = 0;

  /**
   * @brief Computes the volume flux contribution to the cell update.
   *
   * @param[inout] lduh Cell-local update DoF.
   * @param[in]    dx   Extent of the cell in each coordinate direction.
   * @param[dt]    dt   Time step size.
   */
  virtual void volumeIntegral(
      double * lduh,
      const double * const lFhi,
      const tarch::la::Vector<DIMENSIONS,double>&  dx
  ) = 0;

  /**
   * @brief Computes the surface integral contributions
   * to the cell update.
   *
   * @param[inout] lduh   Cell-local update DoF.
   * @param[in]    lFhbnd Cell-local DoF of the boundary extrapolated fluxes.
   * @param[in]    dx     Extent of the cell in each coordinate direction.
   */
  virtual void surfaceIntegral(
      double * lduh,
      const double * const lFhbnd,
      const tarch::la::Vector<DIMENSIONS,double>&  dx
  ) = 0;

  /**
   * @brief Computes the normal fluxes (or fluctuations) at the interface of two cells.
   *
   * @param[inout] FL    Flux DoF belonging to the left cell.
   * @param[inout] FR    Flux DoF belonging the right cell.
   * @param[in]    QL    DoF of the boundary extrapolated predictor belonging to the left cell.
   * @param[in]    QR    DoF of the boundary extrapolated predictor belonging to the right cell.
   * @param[in]    hFace The (d-1)-measure of the interface: Length of an edge in 2D, area of a face in 3D, etc.
   * @param[in]    n     Unit vector normal to the interface.
   */
  virtual void riemannSolver(
      double * FL,
      double * FR,
      const double * const QL,
      const double * const QR,
      const double dt,
      const double hFace,
      const tarch::la::Vector<DIMENSIONS,double>&  n
  ) = 0;

  /**
   * @brief Computes cell-local predictor space-time, volume, and face DoF.
   *
   * Computes the cell-local space-time predictor lQi, the space-time volume flux lFi,
   * the predictor lQhi, the volume flux lFhi, the boundary
   * extrapolated predictor lQhbnd and normal flux lFhbnd.
   *
   * @param[inout] lQi    Space-time predictor DoF.
   * @param[inout] lFi    Space-time flux DoF.
   * @param[inout] lQhi   Predictor DoF
   * @param[inout] lFhi   Volume flux DoF.
   * @param[inout] lQhbnd Boundary extrapolated predictor DoF.
   * @param[inout] lFhbnd Boundary extrapolated normal fluxes (or fluctuations).
   * @param[ino]   luh    Solution DoF.
   * @param[in]    dx     Extent of the cell in each coordinate direction.
   * @param[in]    dt     Time step size.
   */
  virtual void spaceTimePredictor(
      double * lQi,
      double * lFi,
      double * lQhi,
      double * lFhi,
      double * lQhbnd,
      double * lFhbnd,
      const double * const luh,
      const tarch::la::Vector<DIMENSIONS,double>&  dx,
      const double dt
  ) = 0;

  /**
   * @brief Returns a stable time step size.
   *
   * @param[in] luh       Cell-local solution dof
   * @param[in] dx        Extent of the cell in each coordinate direction.
   */
  virtual double stableTimeStepSize(
      const double * const luh,
      const tarch::la::Vector<DIMENSIONS,double>&  dx
  ) = 0;

  /**
   * Sets the initial values.
   *
   * @param[inout] luh     Cell-local solution DoF.
   * @param[in]    center  Element center.
   * @param[in]    dx      Extent of the cell in each coordinate direction.
   */
  virtual void initialValues(
      double * luh,
      const tarch::la::Vector<DIMENSIONS,double>&  center,
      const tarch::la::Vector<DIMENSIONS,double>&  dx
  ) = 0;

  // @todo eher raus

  /**
   * Plot the solution.
   *
   * @param[inout] luh     Cell-local solution DoF.
   * @param[in]    center  Element center.
   * @param[in]    dx      Extent of the cell in each coordinate direction.
   */
  virtual void plot(
      const double * luh,
      const tarch::la::Vector<DIMENSIONS,double>&  center,
      const tarch::la::Vector<DIMENSIONS,double>&  dx
  ) = 0;
};

#endif

