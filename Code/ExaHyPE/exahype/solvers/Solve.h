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
      exahype::solvers::Solve::Type type,
      exahype::solvers::Solve::TimeStepping timeStepping,
      bool correctorTimeLagging,
      bool active
  );

  /**
   * Copy constructor.
   */
  Solve (const Solve& anotherSolve);

  /**
   * Copy constructor.
   */
  Solve (
      std::string name,
      int solverNumber,
      exahype::solvers::Solve::Type type,
      exahype::solvers::Solve::TimeStepping timeStepping,
      bool correctorTimeLagging,
      bool active,
      double correctorTimeStamp,
      double predictorTimeStamp,
      double correctorTimeStepSize,
      double predictorTimeStepSize,
      double nextPredictorTimeStepSize);

  /**
   * Each solve has an identifier/name. It is used for debug purposes only.
   */
  const std::string getName () const;

  /**
   * Each solve has a solver number that says which solver is to be
   * used. Typically this is an ascending index starting from 0.
   */
  const int getSolverNumber () const;

  // @todo 13/02/16:Dominic Etienne Charrier
  // This is just a dummy entry.
  // I think we need a tree-like structure of Solves
  // tracking which solve is a sub solve of which
  // other solve.
  //
  // In an uncertainty identification scenario,
  // sub solves in one cell can request data from
  // their parents in neighbouring cells if they
  // want to perform the Riemann solve.
  // In this case, we probably also have to
  // register new cell descriptions in the neighbouring
  // cells.
  int getParentSolve () const;

  void setParentSolve (int parentSolve);

  void setActive (bool active);

  bool isActive () const;

  bool isCorrectorTimeLagging () const;

  Type getType () const;

  void setType (Type type);

  TimeStepping getTimeStepping () const;

  double getPredictorTimeStamp () const;

  void setPredictorTimeStamp (double predictorTimeStamp);

  double getCorrectorTimeStamp () const;

  void setCorrectorTimeStamp (double correctorTimeStamp);

  /**
   * Returns the maximum stable global corrector time step size for the current sweep.
   * Has no use in the FV method.
   */
  double getCorrectorTimeStepSize () const;

  void setCorrectorTimeStepSize (double correctorTimeStepSize);

  /**
   * Returns the maximum stable global predictor time step size for the current sweep.
   * For the FV method, this represents the time step size.
   */
  double getPredictorTimeStepSize () const;

  void setPredictorTimeStepSize (double predictorTimeStepSize);

  /**
   * Returns the maximum stable global predictor time step size for the next sweep.
   * For the FV method, this represents the next time step size.
   */
    double getNextPredictorTimeStepSize () const;

  /**
   * Computes the maximum stable global predictor time step size for the next sweep
   * as the minimum of the stored next predictor time step size of this Solve
   * and the argument.
   */
  void updateNextPredictorTimeStepSize (const double& nextPredictorTimeStepSize);

  /**
   * Prepares the corrector, predictor, and next predictor time
   * step size for the next sweep.
   */
  void startNewTimeStep ();

  /**
   * Merges this solve with another solve
   */
  void merge (const Solve& otherSolve);

  virtual ~Solve () {}

private:
  /**
   * Each solve has an identifier/name. It is used for debug purposes only.
   */
  std::string        _name;

  /**
   *  Refers to a Solver in the solver registry.
   */
  int                _solverNumber;

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

  bool               _correctorTimeLagging;

  /**
   * Indicates if this solve is active.
   */
  bool               _active;

  /**
   * Minimum corrector time stamp.
   */
  double             _correctorTimeStamp;

  /**
   * Minimum predictor time stamp. Always equal or larger
   * than the minimum corrector time stamp.
   */
  double             _predictorTimeStamp;

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
};

#endif /* SOLVE_H_ */
