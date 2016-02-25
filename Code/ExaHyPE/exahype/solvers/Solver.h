//$<EXAHYPE_HEADER_FILE_COPYRIGHT_NOTE>$

#ifndef _EXAHYPE_SOLVERS_SOLVER_H_
#define _EXAHYPE_SOLVERS_SOLVER_H_


#include <string>
#include <vector>

#include "peano/utils/Globals.h"

#include "tarch/la/Vector.h"

#include "exahype/records/ADERDGCellDescription.h"


#define EXAHYPE_FACE_LEFT   0
#define EXAHYPE_FACE_RIGHT  1
#define EXAHYPE_FACE_FRONT  2
#define EXAHYPE_FACE_BACK   3
#define EXAHYPE_FACE_BOTTOM 4
#define EXAHYPE_FACE_TOP    5

// todo 08/02/16:Dominic Etienne Charrier
// move somewhere else
// is evaluated at compile time
constexpr int power(int basis, int exp) {
  return (exp==0)? 1 : basis * power(basis, exp-1);
}

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
  enum Type {
    ADER_DG
  };

/*
  enum Type {
    SOLVE, SUBSOLVE
  };
*/

  enum TimeStepping {
    GlobalTimeStepping, // ANARCHIC
  };

protected:
  /**
   * Each solver has an identifier/name. It is used for debug purposes only.
   */
  const std::string  _identifier;

  const Type         _type;

  /**
   * Each solver has a kernel number that says which kernel is to be
   * executed. Typically this is an ascending index starting from 0.
   */
  const int         _kernelNumber;

  const int         _numberOfVariables;

  const int         _nodesPerCoordinateAxis;

  const int         _unknownsPerFace;

  const int         _unknownsPerCellBoundary;

  const int         _unknownsPerCell;

  const int         _fluxUnknownsPerCell;

  const int         _spaceTimeUnknownsPerCell;

  const int         _spaceTimeFluxUnknownsPerCell;

  const TimeStepping _timeStepping;

  /**
   * Minimum corrector time stamp.
   */
  double             _minCorrectorTimeStamp;

  /**
   * Minimum predictor time stamp. Always equal or larger
   * than the minimum corrector time stamp.
   */
  double             _minPredictorTimeStamp;

  /**
   * Corrector time step size.
   */
  double             _correctorTimeStepSize;

  /**
   * Predictor time step size.
   */
  double             _predictorTimeStepSize;

  /**
   * Predictor time step size.
   */
  double             _nextPredictorTimeStepSize;
public:
  Solver(const std::string& identifier, Type type, int kernelNumber, int numberOfVariables, int nodesPerCoordinateAxis, TimeStepping timeStepping);

  virtual ~Solver() {}

  /**
   * Identify minimal mesh width at a certain point in the domain. This
   * minimal mesh width is used both as a constraint on the AMR as well
   * as to set up the initial grid. If you return 0, you indicate that
   * this PDE might not exist in the domain.
   */
  virtual int getMinimumTreeDepth() const = 0;

  std::string getIdentifier() const;

  Type getType() const;

  TimeStepping getTimeStepping() const;

  int getNumberOfVariables() const;

  /**
   * This operation returns the number of space time
   * unknowns per cell.
   *
   * Note that this operation might only have a meaning for space-time type
   * discretisation methods.
   */
  int getSpaceTimeUnknownsPerCell() const;

  /**
   * This operation returns the number of space time
   * flux unknowns per cell.
   *
   * Note that this operation might only have a meaning for space-time type
   * discretisation methods.
   */
  int getSpaceTimeFluxUnknownsPerCell() const;

  /**
   * This operation returns the number of unknowns per cell located in
   * the interior of a cell.
   */
  int getUnknownsPerCell() const;

  /**
   * This operation returns the number of flux unknowns per cell
   * located in the interior of a cell.
   */
  int getFluxUnknownsPerCell() const;

  /**
   * This operation returns the number of unknowns that are located
   * on or in the vicinity of the boundary of a cell.
   */
  int getUnknownsPerCellBoundary() const;

  /**
   * This operation returns the number of unknowns that are located
   * on or in the vicinity of each face of a cell.
   */
  int getUnknownsPerFace() const;

  /**
   * If you use a higher order method, then this operation returns the
   * polynomial degree plus one. If you use a Finite Volume method, it
   * returns the number of cells within a patch per coordinate axis.
   */
  int getNodesPerCoordinateAxis() const;

//  /**
//   * @brief Prolongates coarse grid face unknowns
//   * \p levelDifference levels down to the fine grid unknowns.
//   */
//  virtual void prolongateCoarseGridFaceUnknowns(
//      double* fineGridUnknowns,
//      const double* const coarseGridUnknowns,
//      const int levelDifference
//  ) = 0;
//
//  /**
//   * @brief Restricts fine grid face unknowns
//   * \p levelDifference levels up to the coarse grid unknowns.
//   */
//  virtual void updateFineGridFaceUnknowns(
//      double* fineGridUnknowns,
//      const double* const fineGridUnknowns,
//      const int levelDifference
//  ) = 0;

  /**
   * @brief Adds the solution update to the solution.
   *
   * @param[inout] luh  Cell-local solution DoF.
   * @param[in]    lduh Cell-local update DoF.
   * @param[dt]    dt   Time step size.
   */
  virtual void solutionUpdate(
      double * luh,
      const double * const lduh,
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
      const tarch::la::Vector<DIMENSIONS,double>& dx
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
      const tarch::la::Vector<DIMENSIONS,double>& dx
  ) = 0;

  /**
   * @brief Computes the normal fluxes (or fluctuations) at the interface of two cells.
   *
   * @param[inout] FL             Flux DoF belonging to the left cell.
   * @param[inout] FR             Flux DoF belonging the right cell.
   * @param[in]    QL             DoF of the boundary extrapolated predictor belonging to the left cell.
   * @param[in]    QR             DoF of the boundary extrapolated predictor belonging to the right cell.
   * @param[in]    normalNonZero  Index of the nonzero normal vector component, i.e., 0 for e_x, 1 for e_y, and 2 for e_z.
   */
  virtual void riemannSolver(
      double * FL,
      double * FR,
      const double * const QL,
      const double * const QR,
      const double dt,
      const int normalNonZero
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
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double dt
  ) = 0;

  /**
   * @brief Returns a stable time step size.
   *
   * @param[in] luh       Cell-local solution DoF.
   * @param[in] dx        Extent of the cell in each coordinate direction.
   */
  virtual double stableTimeStepSize(
      const double * const luh,
      const tarch::la::Vector<DIMENSIONS,double>& dx
  ) = 0;

  /**
   * @brief Sets the initial values.
   *
   * @param[in] luh    Cell-local solution DoF.
   * @param[in] center Element center.
   * @param[in] dx     Extent of the cell in each coordinate direction.
   */
  virtual void initialCondition(
      double * luh,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx
  ) = 0;

  /**
   * This operation allows you to realise time-dependent conditions.
   * Please be aware that this operation is called per time step if
   * the corresponding predicate hasToUpdateSolution() yields true for the
   * region.
   */
  virtual void updateSolution(
    double *                                      luh,
    const tarch::la::Vector<DIMENSIONS,double>&   center,
    const tarch::la::Vector<DIMENSIONS,double>&   dx,
    double                                        t,
    double                                        dt
  ) = 0;

  virtual bool hasToUpdateSolution(
    const tarch::la::Vector<DIMENSIONS,double>&   center,
    const tarch::la::Vector<DIMENSIONS,double>&   dx
  ) = 0;

  /**
   * @todo Dominic, please add a description what this routine does.
   */
  void synchroniseTimeStepping(exahype::records::ADERDGCellDescription& p) const;

  void startNewTimeStep();

  double getPredictorTimeStamp() const;
};

#endif

