/*
 * ADERDGAPosterioriSubcellLimiter.h
 *
 *  Created on: 6 Oct 2016
 *      Author: dominic
 */

#ifndef ADERDGAPOSTERIORISUBCELLLIMITER_H_
#define ADERDGAPOSTERIORISUBCELLLIMITER_H_

#include "CellWiseCoupling.h"

namespace exahype {
namespace solvers {

class ADERDGAPosterioriSubcellLimiter;

} /* namespace solvers */
} /* namespace exahype */

class exahype::solvers::ADERDGAPosterioriSubcellLimiter : public exahype::solvers::CellWiseCoupling {
private:
  /**
   * Flag is updated in coupleSolversBeforeSolutionUpdate(...)
   * and read in coupleSolversAfterSolutionUpdate(...).
   */
  bool _limiterIsActive;

  /**
   * Element index of the ADER-DG solver whose solution is to limit in
   * registry exahype::solvers::RegisteredSolvers.
   */
  const int  _aderdgSolverNumber;
  /**
   * Element index of the Finite Volumes solver used
   * for the subcell limiting in registry
   * exahype::solvers::RegisteredSolvers.
   */
  const int  _finiteVolumesSolverNumber;

public:
  ADERDGAPosterioriSubcellLimiter(int aderdgSolverNumber,int finiteVolumesSolverNumber);
  virtual ~ADERDGAPosterioriSubcellLimiter() {};

  // Disallow copy and assignment
  ADERDGAPosterioriSubcellLimiter(const ADERDGAPosterioriSubcellLimiter& other) = delete;
  ADERDGAPosterioriSubcellLimiter& operator=(const ADERDGAPosterioriSubcellLimiter& other) = delete;

  /**
   * Always returns true.
   */
  bool isActive(double timeStamp) override;

  /**
   * Couple solvers before the solution update of the solvers is performed.
   */
  void coupleSolversBeforeSolutionUpdate(const int cellDescriptionsIndex) override;

  /**
   * Couple solvers after the solution update of the solvers has been performed.
   */
  void coupleSolversAfterSolutionUpdate(const int cellDescriptionsIndex) override;
};

#endif /* ADERDGAPOSTERIORISUBCELLLIMITER_H_ */
