//$<EXAHYPE_HEADER_FILE_COPYRIGHT_NOTE>$
#ifndef SOLVE_H_
#define SOLVE_H_

#include <string>
#include <vector>

namespace exahype {
  namespace solvers {
    class Solve;
  }
}
class exahype::solvers::Solve {
public:
  const static int InvalidParentSolveIdentifier = -1;

  enum Type {
    SOLVE, SUBSOLVE
  };

  enum TimeStepping {
    GLOBAL //, ANARCHIC
  };

  Solve (
      int solverNumber,
      int parentSolve,
      exahype::solvers::Solve::Type type,
      exahype::solvers::Solve::TimeStepping timeStepping,
      bool useSameTimeStepSize,
      bool active
  );

  /**
   * Copy constructor.
   */
  Solve (const Solve& anotherSolve);

  /**
   * Each solve has a solver number that says which solver is to be
   * used. Typically this is an ascending index starting from 0.
   */
  const int getSolverNumber() const;

  int getParentSolve () const;

  void setParentSolve (int parentSolve);

  void setActive(bool active);

  bool isActive() const;

  bool isSameTimeStepSize() const;

  Type getType() const;

  void setType (Type type);

  TimeStepping getTimeStepping() const;

  double getPredictorMinTimeStamp () const;

  void setPredictorMinTimeStamp (double predictorMinTimeStamp);

  double getCorrectorMinTimeStamp () const;

  void setCorrectorMinTimeStamp (double correctorMinTimeStamp);

  /**
   * Returns the maximum stable global corrector time step size for the current sweep.
   * Has no use in the FV method.
   */
  double getCorrectorTimeStepSize() const;

  /**
   * Returns the maximum stable global predictor time step size for the current sweep.
   * For the FV method, this represents the time step size.
   */
  double getPredictorTimeStepSize() const;

  /**
   * Computes the maximum stable global predictor time step size for the next solve
   * as the minimum of the stored next predictor time step size of this Solve
   * and the argument.
   */
  void updateNextPredictorTimeStepSize(const double& nextPredictorTimeStepSize);

  /**
   * Prepares the corrector, predictor, and next predictor time
   * step size for the next sweep.
   */
  void startNewTimeStep();

  /**
   * Merges this solve with another solve
   */
  void merge(const Solve& anotherSolve);

  virtual ~Solve () {}

private:

  /**
   * Each solve has an identifier/name. It is used for debug purposes only.
   */
  const std::string  _identifier;

  /**
   *  Refers to a Solver in the solver registry.
   */
  const int          _solverNumber;

  /**
   * Refers to the Solve that spawned this solve if this a sub solve.
   *
   *
   * @todo 11/02/16:Dominic Etienne Charrier
   * At the moment the parentSolve must always be
   * of type SOLVE. SUBSOLVE can lead to confusion.
   * The solver registry must be designed carefully
   * if we want to register solves on the fly.
   */
  int                _parentSolve;

  Type               _type;

  TimeStepping       _timeStepping;

  bool               _sameTimeStepSize;

  /**
   * Indicates if this solve is active.
   */
  bool               _active;

  /**
   * Minimum predictor time stamp. Always equal or larger
   * than the minimum corrector time stamp.
   */
  double             _predictorMinTimeStamp;

  /**
   * Minimum corrector time stamp.
   */
  double             _correctorMinTimeStamp;

  /**
   * Predictor time step size.
   */
  double             _predictorTimeStepSize;

  /**
   * Predictor time step size.
   */
  double             _nextPredictorTimeStepSize;

  /**
   * Corrector time step size.
   */
  double             _correctorTimeStepSize;
};

#endif /* SOLVE_H_ */
